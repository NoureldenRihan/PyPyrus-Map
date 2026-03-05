"""
pypyrus_map.graph
=================
Builds and manages the NetworkX DiGraph from schema objects.

This layer owns all graph construction, node/edge merging for chaining,
and the session's internal mutable state. No COBRApy objects, no rendering.
"""

from __future__ import annotations

import networkx as nx

from .schema import (
    MetaboliteNode,
    ReactionNode,
    MetabolicEdge,
    PathwayGraph,
)


# ---------------------------------------------------------------------------
# Node / edge attribute key constants
# ---------------------------------------------------------------------------

NODE_TYPE_KEY = "type"
NODE_DATA_KEY = "data"
EDGE_DATA_KEY = "data"

TYPE_METABOLITE = "metabolite"
TYPE_REACTION = "reaction"


# ---------------------------------------------------------------------------
# Low-level graph construction
# ---------------------------------------------------------------------------

def _add_metabolite_node(g: nx.DiGraph, node: MetaboliteNode) -> None:
    """
    Add or update a metabolite node in the graph.

    If the node already exists (same ID), we update is_anchor but preserve
    all other existing attributes — this is the key deduplication logic for
    chaining. A node that was previously a neighbor can become an anchor
    when the user explicitly adds it, but never the reverse.
    """
    if node.id in g:
        existing: MetaboliteNode = g.nodes[node.id][NODE_DATA_KEY]
        if node.is_anchor:
            existing.is_anchor = True
    else:
        g.add_node(
            node.id,
            **{NODE_TYPE_KEY: TYPE_METABOLITE, NODE_DATA_KEY: node},
        )


def _add_reaction_node(g: nx.DiGraph, node: ReactionNode) -> None:
    """
    Add or update a reaction node in the graph.

    If the reaction already exists (can happen when two metabolites share a
    reaction), we update flux from the newer extraction in case the solution
    was provided on a second add() call. All other attributes are preserved.
    """
    if node.id in g:
        existing: ReactionNode = g.nodes[node.id][NODE_DATA_KEY]
        if node.flux is not None:
            existing.flux = node.flux
    else:
        g.add_node(
            node.id,
            **{NODE_TYPE_KEY: TYPE_REACTION, NODE_DATA_KEY: node},
        )


def _add_edge(g: nx.DiGraph, edge: MetabolicEdge) -> None:
    """
    Add a directed edge to the graph.

    Duplicate edges (same source → target) are skipped silently — this is
    the correct behavior when two metabolite neighborhoods share a reaction,
    as the edge is already present from the first add().
    """
    if not g.has_edge(edge.source_id, edge.target_id):
        g.add_edge(
            edge.source_id,
            edge.target_id,
            **{EDGE_DATA_KEY: edge},
        )


# ---------------------------------------------------------------------------
# Session graph state
# ---------------------------------------------------------------------------

class GraphState:
    """
    Mutable internal state of a PyPyrusMap Map session.

    Holds the growing NetworkX DiGraph across multiple add() calls.
    Encapsulates all graph mutation logic so the session class stays clean.
    """

    def __init__(self) -> None:
        self._g: nx.DiGraph = nx.DiGraph()
        self._anchor_ids: list[str] = []  # ordered, preserves add() sequence

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def merge(
        self,
        focal_node: MetaboliteNode,
        reaction_nodes: list[ReactionNode],
        neighbor_nodes: list[MetaboliteNode],
        edges: list[MetabolicEdge],
        is_anchor: bool = True,
    ) -> None:
        """
        Merge a depth-1 neighborhood query result into the session graph.

        This is called once per session.add() call. The focal_node is always
        marked as an anchor. All nodes and edges are unioned with deduplication.

        Parameters
        ----------
        focal_node:
            The metabolite that was explicitly queried.
        reaction_nodes:
            Reactions extracted by the query layer.
        neighbor_nodes:
            Neighboring metabolites extracted by the query layer.
        edges:
            Directed edges from the query layer.
        is_anchor:
            Whether to mark the focal node as an anchor (always True for
            direct add() calls; False for internal test use only).
        """
        focal_node.is_anchor = is_anchor

        _add_metabolite_node(self._g, focal_node)
        if is_anchor and focal_node.id not in self._anchor_ids:
            self._anchor_ids.append(focal_node.id)

        for rxn_node in reaction_nodes:
            _add_reaction_node(self._g, rxn_node)

        for nbr_node in neighbor_nodes:
            _add_metabolite_node(self._g, nbr_node)

        for edge in edges:
            _add_edge(self._g, edge)

    def remove_metabolite(self, metabolite_id: str) -> bool:
        """
        Remove a metabolite node and all its edges from the graph.

        Also removes any reaction nodes that become orphaned (i.e., connected
        to no metabolite nodes after this removal).

        Returns True if the node was found and removed, False otherwise.
        """
        if metabolite_id not in self._g:
            return False

        # Collect reaction neighbors before removal
        reaction_neighbors = [
            n for n in list(self._g.predecessors(metabolite_id)) +
                        list(self._g.successors(metabolite_id))
            if self._g.nodes[n].get(NODE_TYPE_KEY) == TYPE_REACTION
        ]

        self._g.remove_node(metabolite_id)

        if metabolite_id in self._anchor_ids:
            self._anchor_ids.remove(metabolite_id)

        # Remove orphaned reactions (reactions with no remaining metabolite neighbors)
        for rxn_id in reaction_neighbors:
            if rxn_id in self._g:
                all_neighbors = (
                    list(self._g.predecessors(rxn_id)) +
                    list(self._g.successors(rxn_id))
                )
                if not all_neighbors:
                    self._g.remove_node(rxn_id)

        return True

    def clear(self) -> None:
        """Reset the session to an empty state."""
        self._g.clear()
        self._anchor_ids.clear()

    # ------------------------------------------------------------------
    # Read
    # ------------------------------------------------------------------

    @property
    def is_empty(self) -> bool:
        return len(self._g.nodes) == 0

    def build_pathway_graph(self) -> PathwayGraph:
        """
        Snapshot the current mutable state into an immutable PathwayGraph.

        Returns a PathwayGraph that can be safely passed to the render layer
        without risk of mutation.
        """
        has_flux = any(
            self._g.nodes[n][NODE_DATA_KEY].flux is not None
            for n in self._g.nodes
            if self._g.nodes[n].get(NODE_TYPE_KEY) == TYPE_REACTION
        )

        return PathwayGraph(
            nx_graph=self._g.copy(),  # shallow copy — node data objects are shared but graph structure is independent
            anchor_ids=list(self._anchor_ids),
            has_flux=has_flux,
        )
