"""
pypyrus_map.query
=================
The query layer is the only place in PyPyrusMap that touches COBRApy objects.

Everything returned from this module is in terms of PyPyrusMap schema
dataclasses — no cobra.Metabolite, cobra.Reaction, or cobra.Solution objects
leak past this boundary.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional

from .schema import MetaboliteNode, ReactionNode, MetabolicEdge
from .exceptions import MetaboliteNotFoundError

if TYPE_CHECKING:
    # Type hints only — avoids hard import at module load for environments
    # where COBRApy might not be installed yet.
    import cobra


# ---------------------------------------------------------------------------
# Currency metabolites — suppressed by default in pathway tracing
# These are the "highway metabolites" that connect everything to everything.
# This list follows the conventions in Huss & Holme (2007) and KEGG RPAIR.
# ---------------------------------------------------------------------------

DEFAULT_CURRENCY_IDS: frozenset[str] = frozenset({
    # Adenine nucleotides
    "atp_c", "adp_c", "amp_c", "atp_e", "adp_e",
    # Nicotinamide
    "nad_c", "nadh_c", "nadp_c", "nadph_c",
    # Protons and water
    "h_c", "h_e", "h_m", "h2o_c", "h2o_e", "h2o_m",
    # Phosphate
    "pi_c", "pi_e", "ppi_c",
    # CoA
    "coa_c", "coa_m",
    # CO2
    "co2_c", "co2_e",
    # Oxygen
    "o2_c", "o2_e",
    # Ubiquinone
    "q8_c", "q8h2_c",
    # Ammonia
    "nh4_c", "nh4_e",
})


def resolve_metabolite(
    model: "cobra.Model",
    metabolite_id: str,
) -> "cobra.Metabolite":
    """
    Look up a metabolite by exact BiGG ID in a COBRApy model.

    Parameters
    ----------
    model:
        A loaded COBRApy model.
    metabolite_id:
        Exact BiGG metabolite ID including compartment suffix, e.g. 'glc__D_c'.

    Returns
    -------
    cobra.Metabolite

    Raises
    ------
    MetaboliteNotFoundError
        If the ID is not present. Includes candidate suggestions found by
        prefix matching to help the user self-correct.
    """
    try:
        return model.metabolites.get_by_id(metabolite_id)
    except KeyError:
        # Build a short candidate list for the error message.
        # Match on the base name (strip compartment suffix) for helpfulness.
        base = metabolite_id.rsplit("_", 1)[0].lower()
        candidates = [
            m.id for m in model.metabolites
            if m.id.lower().startswith(base)
        ][:8]
        raise MetaboliteNotFoundError(
            metabolite_id=metabolite_id,
            model_id=model.id,
            candidates=candidates,
        ) from None


def extract_metabolite_node(met: "cobra.Metabolite") -> MetaboliteNode:
    """Convert a cobra.Metabolite into a MetaboliteNode schema object."""
    return MetaboliteNode(
        id=met.id,
        name=met.name or met.id,
        compartment=met.compartment or "",
        formula=getattr(met, "formula", None),
    )


def extract_reaction_node(
    rxn: "cobra.Reaction",
    solution: Optional["cobra.Solution"],
) -> ReactionNode:
    """
    Convert a cobra.Reaction into a ReactionNode schema object.

    Parameters
    ----------
    rxn:
        The COBRApy reaction object.
    solution:
        Optional FBA solution. If provided, flux is read from it.
        If None, ReactionNode.flux remains None.
    """
    # Resolve gene names — prefer gene.name over gene.id (locus tag)
    gene_names = []
    for gene in rxn.genes:
        name = getattr(gene, "name", None)
        gene_names.append(name if name and name.strip() else gene.id)

    flux: Optional[float] = None
    if solution is not None:
        try:
            flux = float(solution.fluxes[rxn.id])
        except (KeyError, AttributeError):
            flux = None

    return ReactionNode(
        id=rxn.id,
        name=rxn.name or rxn.id,
        subsystem=rxn.subsystem or "",
        genes=gene_names,
        lower_bound=float(rxn.lower_bound),
        upper_bound=float(rxn.upper_bound),
        flux=flux,
    )


def get_neighborhood(
    model: "cobra.Model",
    metabolite_id: str,
    solution: Optional["cobra.Solution"] = None,
    suppress_currency: bool = True,
    custom_currency_ids: Optional[frozenset[str]] = None,
) -> tuple[MetaboliteNode, list[ReactionNode], list[MetaboliteNode], list[MetabolicEdge]]:
    """
    Extract the depth-1 neighborhood of a metabolite from a COBRApy model.

    This is the core query: for a given metabolite, find every reaction that
    produces or consumes it, then find all OTHER metabolites involved in those
    reactions (the neighbors).

    Parameters
    ----------
    model:
        A loaded COBRApy model.
    metabolite_id:
        Exact BiGG metabolite ID.
    solution:
        Optional FBA solution for flux annotation.
    suppress_currency:
        If True (default), currency metabolites (ATP, NADH, H2O, etc.) are
        excluded from the neighbor metabolite list. The focal metabolite itself
        is always included regardless of this flag.
    custom_currency_ids:
        If provided, replaces the default currency ID set entirely.

    Returns
    -------
    A tuple of:
        focal_node      - MetaboliteNode for the queried metabolite
        reaction_nodes  - All reactions involving this metabolite
        neighbor_nodes  - All OTHER metabolites in those reactions
        edges           - All directed MetabolicEdge objects
    """
    currency_ids = custom_currency_ids if custom_currency_ids is not None else DEFAULT_CURRENCY_IDS

    cobra_met = resolve_metabolite(model, metabolite_id)
    focal_node = extract_metabolite_node(cobra_met)

    reaction_nodes: list[ReactionNode] = []
    neighbor_nodes_map: dict[str, MetaboliteNode] = {}
    edges: list[MetabolicEdge] = []

    for rxn in cobra_met.reactions:
        rxn_node = extract_reaction_node(rxn, solution)
        reaction_nodes.append(rxn_node)

        # Determine the stoichiometric coefficient for the focal metabolite.
        focal_coeff = rxn.metabolites[cobra_met]  # negative = consumed, positive = produced

        # --- Edges for the focal metabolite ---
        if focal_coeff < 0:
            # Focal metabolite is a reactant → metabolite → reaction edge
            edges.append(MetabolicEdge(
                source_id=focal_node.id,
                target_id=rxn_node.id,
                stoichiometry=abs(focal_coeff),
                direction="consumes",
            ))
        else:
            # Focal metabolite is a product → reaction → metabolite edge
            edges.append(MetabolicEdge(
                source_id=rxn_node.id,
                target_id=focal_node.id,
                stoichiometry=abs(focal_coeff),
                direction="produces",
            ))

        # --- Neighbor metabolites in this reaction ---
        for other_met, coeff in rxn.metabolites.items():
            if other_met.id == focal_node.id:
                continue
            if suppress_currency and other_met.id in currency_ids:
                continue

            if other_met.id not in neighbor_nodes_map:
                neighbor_nodes_map[other_met.id] = extract_metabolite_node(other_met)

            neighbor_node = neighbor_nodes_map[other_met.id]

            # Build edge for this neighbor relative to the reaction.
            if coeff < 0:
                # Neighbor is a reactant of the reaction
                edges.append(MetabolicEdge(
                    source_id=neighbor_node.id,
                    target_id=rxn_node.id,
                    stoichiometry=abs(coeff),
                    direction="consumes",
                ))
            else:
                # Neighbor is a product of the reaction
                edges.append(MetabolicEdge(
                    source_id=rxn_node.id,
                    target_id=neighbor_node.id,
                    stoichiometry=abs(coeff),
                    direction="produces",
                ))

    neighbor_nodes = list(neighbor_nodes_map.values())
    return focal_node, reaction_nodes, neighbor_nodes, edges
