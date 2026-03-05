"""
tests/test_graph.py
===================
Tests for the graph layer (GraphState, PathwayGraph) and session-level
behaviors: chaining, deduplication, removal, and error handling.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pypyrus_map.graph import GraphState, TYPE_METABOLITE, TYPE_REACTION, NODE_TYPE_KEY, NODE_DATA_KEY
from pypyrus_map.query import get_neighborhood
from pypyrus_map.session import PyPyrusMap
from pypyrus_map.exceptions import EmptySessionError, MetaboliteNotFoundError


# ---------------------------------------------------------------------------
# GraphState tests
# ---------------------------------------------------------------------------

class TestGraphState:

    def _get_ggpp_neighborhood(self, mock_model, mock_solution=None):
        return get_neighborhood(mock_model, "ggpp_c", solution=mock_solution)

    def _get_phyto_neighborhood(self, mock_model, mock_solution=None):
        return get_neighborhood(mock_model, "phyto_c", solution=mock_solution)

    def test_single_merge_populates_graph(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        assert not state.is_empty
        g = state.build_pathway_graph().nx_graph
        assert "ggpp_c" in g
        assert "ISPA" in g
        assert "CRTB" in g

    def test_anchor_flag_set_on_focal(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        pg = state.build_pathway_graph()
        ggpp_data = pg.nx_graph.nodes["ggpp_c"][NODE_DATA_KEY]
        assert ggpp_data.is_anchor is True

    def test_neighbor_not_anchor(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        pg = state.build_pathway_graph()
        phyto_data = pg.nx_graph.nodes["phyto_c"][NODE_DATA_KEY]
        assert phyto_data.is_anchor is False

    def test_chaining_no_node_duplication(self, mock_model):
        """Adding GGPP then Phytoene should not duplicate GGPP node."""
        state = GraphState()

        focal1, rxns1, nbrs1, edges1 = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal1, rxns1, nbrs1, edges1, is_anchor=True)

        focal2, rxns2, nbrs2, edges2 = self._get_phyto_neighborhood(mock_model)
        state.merge(focal2, rxns2, nbrs2, edges2, is_anchor=True)

        pg = state.build_pathway_graph()
        g = pg.nx_graph

        # GGPP appears once
        assert list(g.nodes).count("ggpp_c") == 1

        # CRTB reaction shared by both neighborhoods — appears once
        assert list(g.nodes).count("CRTB") == 1

    def test_chaining_promotes_neighbor_to_anchor(self, mock_model):
        """
        Phytoene starts as a neighbor of GGPP. When explicitly added,
        it should be promoted to anchor without losing existing edges.
        """
        state = GraphState()

        focal1, rxns1, nbrs1, edges1 = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal1, rxns1, nbrs1, edges1, is_anchor=True)

        focal2, rxns2, nbrs2, edges2 = self._get_phyto_neighborhood(mock_model)
        state.merge(focal2, rxns2, nbrs2, edges2, is_anchor=True)

        pg = state.build_pathway_graph()
        phyto_data = pg.nx_graph.nodes["phyto_c"][NODE_DATA_KEY]
        assert phyto_data.is_anchor is True

    def test_chaining_lycopene_reachable(self, mock_model):
        """After adding GGPP + Phytoene, Lycopene should be in the graph."""
        state = GraphState()

        for met_id in ["ggpp_c", "phyto_c"]:
            focal, rxns, nbrs, edges = get_neighborhood(mock_model, met_id)
            state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        pg = state.build_pathway_graph()
        assert "lyco_c" in pg.nx_graph

    def test_anchor_ids_ordered(self, mock_model):
        """anchor_ids should preserve add() order."""
        state = GraphState()

        for met_id in ["ggpp_c", "phyto_c"]:
            focal, rxns, nbrs, edges = get_neighborhood(mock_model, met_id)
            state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        pg = state.build_pathway_graph()
        assert pg.anchor_ids == ["ggpp_c", "phyto_c"]

    def test_remove_metabolite(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)

        state.remove_metabolite("ggpp_c")
        pg = state.build_pathway_graph()
        assert "ggpp_c" not in pg.nx_graph

    def test_remove_nonexistent_returns_false(self, mock_model):
        state = GraphState()
        result = state.remove_metabolite("not_in_graph_c")
        assert result is False

    def test_clear_empties_state(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)
        state.clear()
        assert state.is_empty

    def test_has_flux_true_when_solution_provided(self, mock_model, mock_solution):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model, mock_solution)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)
        pg = state.build_pathway_graph()
        assert pg.has_flux is True

    def test_has_flux_false_without_solution(self, mock_model):
        state = GraphState()
        focal, rxns, nbrs, edges = self._get_ggpp_neighborhood(mock_model)
        state.merge(focal, rxns, nbrs, edges, is_anchor=True)
        pg = state.build_pathway_graph()
        assert pg.has_flux is False


# ---------------------------------------------------------------------------
# PyPyrusMap Map session tests
# ---------------------------------------------------------------------------

class TestPyPyrusSession:

    def test_basic_add_and_summary(self, mock_model):
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c")
        assert "ggpp_c" in session.summary()

    def test_chaining_add_method(self, mock_model):
        session = PyPyrusMap(mock_model)
        result = session.add("ggpp_c").add("phyto_c")
        # Method chaining returns self
        assert result is session

    def test_invalid_metabolite_raises(self, mock_model):
        session = PyPyrusMap(mock_model)
        with pytest.raises(MetaboliteNotFoundError):
            session.add("fake_met_c")

    def test_render_raises_on_empty(self, mock_model):
        session = PyPyrusMap(mock_model)
        with pytest.raises(EmptySessionError):
            session.render()

    def test_build_raises_on_empty(self, mock_model):
        session = PyPyrusMap(mock_model)
        with pytest.raises(EmptySessionError):
            session.build()

    def test_build_returns_pathway_graph(self, mock_model):
        from pypyrus_map.schema import PathwayGraph
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c")
        pg = session.build()
        assert isinstance(pg, PathwayGraph)

    def test_anchor_ids_property(self, mock_model):
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c").add("phyto_c")
        assert session.anchor_ids == ["ggpp_c", "phyto_c"]

    def test_clear_resets_session(self, mock_model):
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c")
        session.clear()
        with pytest.raises(EmptySessionError):
            session.build()

    def test_remove_warns_on_missing(self, mock_model):
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c")
        with pytest.warns(UserWarning, match="not in the current session"):
            session.remove("not_real_c")

    def test_solution_stored_on_session(self, mock_model, mock_solution):
        session = PyPyrusMap(mock_model, solution=mock_solution)
        assert session.solution is mock_solution

    def test_flux_in_build_when_solution_provided(self, mock_model, mock_solution):
        session = PyPyrusMap(mock_model, solution=mock_solution)
        session.add("ggpp_c")
        pg = session.build()
        assert pg.has_flux is True

    def test_summary_string_format(self, mock_model):
        session = PyPyrusMap(mock_model)
        session.add("ggpp_c")
        s = session.summary()
        assert "Model" in s
        assert "Anchors" in s
        assert "Graph" in s
