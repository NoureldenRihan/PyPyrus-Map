"""
tests/test_query.py
===================
Tests for the query layer — metabolite resolution, neighborhood extraction,
currency suppression, and error handling.
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pypyrus_map.query import resolve_metabolite, get_neighborhood, DEFAULT_CURRENCY_IDS
from pypyrus_map.exceptions import MetaboliteNotFoundError


class TestResolveMetabolite:

    def test_exact_id_found(self, mock_model):
        met = resolve_metabolite(mock_model, "ggpp_c")
        assert met.id == "ggpp_c"

    def test_missing_id_raises_error(self, mock_model):
        with pytest.raises(MetaboliteNotFoundError) as exc_info:
            resolve_metabolite(mock_model, "xyz_not_real_c")
        assert "xyz_not_real_c" in str(exc_info.value)

    def test_error_contains_model_id(self, mock_model):
        with pytest.raises(MetaboliteNotFoundError) as exc_info:
            resolve_metabolite(mock_model, "missing_c")
        assert mock_model.id in str(exc_info.value)

    def test_error_includes_candidates_on_prefix_match(self, mock_model):
        # "ggpp" prefix should suggest "ggpp_c"
        with pytest.raises(MetaboliteNotFoundError) as exc_info:
            resolve_metabolite(mock_model, "ggpp_m")  # wrong compartment
        err = exc_info.value
        assert "ggpp_c" in err.candidates

    def test_case_sensitive_no_fuzzy(self, mock_model):
        with pytest.raises(MetaboliteNotFoundError):
            resolve_metabolite(mock_model, "GGPP_C")  # wrong case


class TestGetNeighborhood:

    def test_ggpp_has_two_reactions(self, mock_model):
        """GGPP is involved in ISPA (produced) and CRTB (consumed)."""
        focal, reactions, neighbors, edges = get_neighborhood(mock_model, "ggpp_c")
        reaction_ids = {r.id for r in reactions}
        assert "ISPA" in reaction_ids
        assert "CRTB" in reaction_ids

    def test_focal_node_id_matches(self, mock_model):
        focal, _, _, _ = get_neighborhood(mock_model, "ggpp_c")
        assert focal.id == "ggpp_c"

    def test_neighbors_exclude_focal(self, mock_model):
        focal, _, neighbors, _ = get_neighborhood(mock_model, "ggpp_c")
        neighbor_ids = {n.id for n in neighbors}
        assert focal.id not in neighbor_ids

    def test_currency_suppression_default(self, mock_model):
        """By default, atp_c and adp_c should be suppressed in HEX1's neighborhood."""
        focal, _, neighbors, _ = get_neighborhood(mock_model, "glc__D_c")
        neighbor_ids = {n.id for n in neighbors}
        assert "atp_c" not in neighbor_ids
        assert "adp_c" not in neighbor_ids
        assert "g6p_c" in neighbor_ids  # non-currency should appear

    def test_currency_suppression_disabled(self, mock_model):
        """With suppress_currency=False, ATP and ADP should appear."""
        focal, _, neighbors, _ = get_neighborhood(
            mock_model, "glc__D_c", suppress_currency=False
        )
        neighbor_ids = {n.id for n in neighbors}
        assert "atp_c" in neighbor_ids
        assert "adp_c" in neighbor_ids

    def test_edges_have_correct_directions(self, mock_model):
        """GGPP is produced by ISPA and consumed by CRTB."""
        _, _, _, edges = get_neighborhood(mock_model, "ggpp_c")
        edge_map = {(e.source_id, e.target_id): e for e in edges}

        # ISPA produces GGPP → edge: ISPA → ggpp_c
        assert ("ISPA", "ggpp_c") in edge_map
        assert edge_map[("ISPA", "ggpp_c")].direction == "produces"

        # CRTB consumes GGPP → edge: ggpp_c → CRTB
        assert ("ggpp_c", "CRTB") in edge_map
        assert edge_map[("ggpp_c", "CRTB")].direction == "consumes"

    def test_flux_attached_when_solution_provided(self, mock_model, mock_solution):
        _, reactions, _, _ = get_neighborhood(
            mock_model, "ggpp_c", solution=mock_solution
        )
        rxn_map = {r.id: r for r in reactions}
        assert rxn_map["ISPA"].flux == pytest.approx(2.41)
        assert rxn_map["CRTB"].flux == pytest.approx(1.20)

    def test_flux_is_none_without_solution(self, mock_model):
        _, reactions, _, _ = get_neighborhood(mock_model, "ggpp_c")
        for rxn in reactions:
            assert rxn.flux is None

    def test_blocked_reaction_detected(self, mock_model, mock_solution):
        """CRTI carries 0 flux — should be flagged as blocked."""
        _, reactions, _, _ = get_neighborhood(
            mock_model, "phyto_c", solution=mock_solution
        )
        rxn_map = {r.id: r for r in reactions}
        assert rxn_map["CRTI"].is_blocked is True
        assert rxn_map["CRTB"].is_blocked is False

    def test_gene_names_extracted(self, mock_model):
        _, reactions, _, _ = get_neighborhood(mock_model, "ggpp_c")
        rxn_map = {r.id: r for r in reactions}
        assert "ispA" in rxn_map["ISPA"].genes
        assert "crtB" in rxn_map["CRTB"].genes

    def test_invalid_metabolite_raises(self, mock_model):
        with pytest.raises(MetaboliteNotFoundError):
            get_neighborhood(mock_model, "not_real_c")
