"""
pypyrus_map.schema
==================
Central data contracts for PyPyrusMap.

All three layers (query, graph, render) communicate exclusively through
these dataclasses. COBRApy objects never leak past the query layer.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal, Optional
import networkx as nx


# ---------------------------------------------------------------------------
# Node types
# ---------------------------------------------------------------------------

@dataclass
class MetaboliteNode:
    """Represents a metabolite vertex in the pathway graph."""

    id: str
    """BiGG metabolite ID, e.g. 'glc__D_c'."""

    name: str
    """Human-readable name, e.g. 'D-Glucose'."""

    compartment: str
    """Compartment symbol, e.g. 'c', 'e', 'm'."""

    formula: Optional[str] = None
    """Chemical formula if available in the model."""

    is_anchor: bool = False
    """True if this node was explicitly added by the user via session.add()."""

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, MetaboliteNode) and self.id == other.id


@dataclass
class ReactionNode:
    """Represents a reaction vertex in the bipartite pathway graph."""

    id: str
    """BiGG reaction ID, e.g. 'PGI'."""

    name: str
    """Human-readable reaction name."""

    subsystem: str
    """Metabolic subsystem / pathway, e.g. 'Glycolysis'."""

    genes: list[str]
    """Human-readable gene names (not locus tags) for this reaction."""

    lower_bound: float
    """Reaction lower bound from the model."""

    upper_bound: float
    """Reaction upper bound from the model."""

    flux: Optional[float] = None
    """Flux value from a COBRApy Solution, None if no solution was provided."""

    @property
    def is_reversible(self) -> bool:
        """True if the reaction can carry flux in both directions."""
        return self.lower_bound < 0 and self.upper_bound > 0

    @property
    def is_blocked(self) -> bool:
        """
        True if a solution was provided AND flux is effectively zero.
        Reactions with no solution attached are never marked blocked.
        """
        if self.flux is None:
            return False
        return abs(self.flux) < 1e-9

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, ReactionNode) and self.id == other.id


# ---------------------------------------------------------------------------
# Edge type
# ---------------------------------------------------------------------------

@dataclass
class MetabolicEdge:
    """
    Directed edge in the bipartite graph.

    In a bipartite metabolic graph, edges always connect a MetaboliteNode
    to a ReactionNode or vice versa:

        MetaboliteNode --> ReactionNode   (metabolite is consumed by reaction)
        ReactionNode   --> MetaboliteNode (reaction produces metabolite)
    """

    source_id: str
    """ID of the source node (metabolite or reaction)."""

    target_id: str
    """ID of the target node (metabolite or reaction)."""

    stoichiometry: float
    """Absolute stoichiometric coefficient for this metabolite in this reaction."""

    direction: Literal["consumes", "produces"]
    """
    'consumes' → metabolite flows into the reaction (metabolite → reaction edge).
    'produces' → reaction flows into the metabolite (reaction → metabolite edge).
    """


# ---------------------------------------------------------------------------
# Top-level graph container
# ---------------------------------------------------------------------------

@dataclass
class PathwayGraph:
    """
    The complete, structured output of a PyPyrusMap Map session.

    Returned by session.build() and consumed by the render layer.
    Can also be serialized for reproducible figures.
    """

    nx_graph: nx.DiGraph
    """
    The underlying NetworkX DiGraph.

    Node attributes:
      - 'type'   : 'metabolite' | 'reaction'
      - 'data'   : MetaboliteNode | ReactionNode instance

    Edge attributes:
      - 'data'   : MetabolicEdge instance
    """

    anchor_ids: list[str]
    """
    Metabolite IDs that were explicitly added by the user via session.add().
    Used by the renderer to apply distinct visual styling to focal nodes.
    """

    has_flux: bool = False
    """True if at least one reaction node carries a non-None flux value."""

    @property
    def metabolite_nodes(self) -> list[MetaboliteNode]:
        return [
            d["data"]
            for _, d in self.nx_graph.nodes(data=True)
            if d.get("type") == "metabolite"
        ]

    @property
    def reaction_nodes(self) -> list[ReactionNode]:
        return [
            d["data"]
            for _, d in self.nx_graph.nodes(data=True)
            if d.get("type") == "reaction"
        ]

    def summary(self) -> str:
        """Human-readable one-line summary of the graph."""
        m = len(self.metabolite_nodes)
        r = len(self.reaction_nodes)
        blocked = sum(1 for rn in self.reaction_nodes if rn.is_blocked)
        flux_str = f", {blocked} blocked" if self.has_flux else ", no flux data"
        return f"PathwayGraph: {m} metabolites, {r} reactions{flux_str}"
