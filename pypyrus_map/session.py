"""
pypyrus_map.session
===================
The PyPyrusMap Map session — the primary object the user interacts with.

All public API flows through this class. Internally it delegates to:
  - query.py  → COBRApy extraction
  - graph.py  → graph construction and union
  - render.py → Matplotlib visualization
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Literal

from .query import get_neighborhood, DEFAULT_CURRENCY_IDS
from .graph import GraphState
from .render import render_pathway_graph, save_figure
from .schema import PathwayGraph
from .exceptions import EmptySessionError, MetaboliteNotFoundError

if TYPE_CHECKING:
    import cobra
    import matplotlib.pyplot as plt


class PyPyrusMap:
    """
    A stateful session for exploring and visualizing metabolic neighborhoods.

    Usage
    -----
    >>> import cobra
    >>> from pypyrus_map import PyPyrusMap
    >>>
    >>> model = cobra.io.load_json_model("iML1515.json")
    >>> solution = model.optimize()
    >>>
    >>> session = PyPyrusMap(model, solution=solution)
    >>> session.add("ggpp_c")
    >>> session.add("phyto_c")   # chains seamlessly — GGPP node not duplicated
    >>> fig = session.render()
    >>> session.export("carotenoid_pathway.svg")

    Parameters
    ----------
    model:
        A loaded COBRApy model. The model is not modified by PyPyrusMap.
    solution:
        Optional COBRApy Solution object (from model.optimize()). If provided,
        flux values are overlaid on all reaction nodes and edges. Fixed at
        session construction — run a new session to compare different solutions.
    suppress_currency:
        If True (default), currency metabolites (ATP, NADH, H₂O, etc.) are
        excluded from neighbor nodes. They are present in reactions but do not
        clutter the graph. The focal metabolite is always shown regardless.
    custom_currency_ids:
        Optional frozenset of metabolite IDs to treat as currency. Replaces
        the built-in default list entirely.
    """

    def __init__(
        self,
        model: "cobra.Model",
        solution: Optional["cobra.Solution"] = None,
        suppress_currency: bool = True,
        custom_currency_ids: Optional[frozenset[str]] = None,
    ) -> None:
        self._model = model
        self._solution = solution
        self._suppress_currency = suppress_currency
        self._currency_ids = custom_currency_ids  # None → use default in query layer

        self._state = GraphState()

        # Track the last rendered figure for convenience
        self._last_fig: Optional["plt.Figure"] = None

    # ------------------------------------------------------------------
    # Core session API
    # ------------------------------------------------------------------

    def add(self, metabolite_id: str) -> "PyPyrusMap":
        """
        Add a metabolite and its depth-1 neighborhood to the session graph.

        If the metabolite was already added in a previous call (e.g. it appeared
        as a neighbor of another metabolite), its node is promoted to 'anchor'
        status and its neighborhood is expanded. Existing nodes and edges are
        never duplicated.

        Parameters
        ----------
        metabolite_id:
            Exact BiGG metabolite ID including compartment suffix, e.g. 'ggpp_c'.

        Returns
        -------
        self — for optional method chaining: session.add("a").add("b").render()

        Raises
        ------
        MetaboliteNotFoundError
            If the ID is not found in the model.
        """
        focal, reactions, neighbors, edges = get_neighborhood(
            model=self._model,
            metabolite_id=metabolite_id,
            solution=self._solution,
            suppress_currency=self._suppress_currency,
            custom_currency_ids=self._currency_ids,
        )
        self._state.merge(focal, reactions, neighbors, edges, is_anchor=True)
        return self

    def remove(self, metabolite_id: str) -> "PyPyrusMap":
        """
        Remove a metabolite node (and its orphaned reactions) from the session.

        Parameters
        ----------
        metabolite_id:
            Exact BiGG metabolite ID to remove.

        Returns
        -------
        self
        """
        found = self._state.remove_metabolite(metabolite_id)
        if not found:
            import warnings
            warnings.warn(
                f"Metabolite '{metabolite_id}' is not in the current session graph. "
                f"Nothing was removed.",
                UserWarning,
                stacklevel=2,
            )
        return self

    def clear(self) -> "PyPyrusMap":
        """
        Reset the session to an empty state.

        Returns
        -------
        self
        """
        self._state.clear()
        self._last_fig = None
        return self

    # ------------------------------------------------------------------
    # Rendering API
    # ------------------------------------------------------------------

    def render(
        self,
        layout: Literal["dot", "spring", "circular"] = "dot",
        orientation: Literal["portrait", "landscape"] = "landscape",
        figsize: Optional[tuple[float, float]] = None,
        title: Optional[str] = None,
        label_mode: Literal["id", "truncate", "full"] = "id",
        show_stoichiometry: bool = False,
        show_flux_labels: bool = False,
        legend_loc: str = "upper left",
        dpi: int = 120,
        font_family: str = "DejaVu Sans",
    ) -> "plt.Figure":
        """
        Render the current session graph as a Matplotlib figure.

        Parameters
        ----------
        layout:
            'dot' (hierarchical, recommended) | 'spring' | 'circular'
        orientation:
            'portrait'  → reactions flow top-to-bottom (default, best for chains).
            'landscape' → reactions flow left-to-right (best for wide graphs).
        figsize:
            Figure size in inches (width, height). Auto-sized from node count
            if not provided — canvas grows with the graph so nothing gets cramped.
        title:
            Auto-generated from anchor IDs if None.
        label_mode:
            'id'       → BiGG ID (default — compact, precise, no overlap).
            'truncate' → Human name, truncated to 18 characters.
            'full'     → Full human name (may overlap on dense graphs).
        show_stoichiometry:
            Add stoichiometric coefficients to edge midpoint labels.
        show_flux_labels:
            Show flux value (4 dp) at edge midpoints when solution is attached.
        legend_loc:
            Matplotlib legend location string. Default 'upper left'. Options:
            'upper right', 'lower left', 'lower right', 'lower center',
            'upper center', 'right', 'center left', 'center right'.
        dpi:
            Screen resolution. SVG/PDF export is always lossless.
        font_family:
            'Arial' matches Nature/Cell figure guidelines.

        Returns
        -------
        matplotlib.figure.Figure

        Raises
        ------
        EmptySessionError
        """
        if self._state.is_empty:
            raise EmptySessionError()

        pg = self._state.build_pathway_graph()
        fig = render_pathway_graph(
            pg=pg,
            layout=layout,
            orientation=orientation,
            figsize=figsize,
            title=title,
            label_mode=label_mode,
            show_stoichiometry=show_stoichiometry,
            show_flux_labels=show_flux_labels,
            legend_loc=legend_loc,
            dpi=dpi,
            font_family=font_family,
        )
        self._last_fig = fig
        return fig

    def export(
        self,
        path: str | Path,
        dpi: int = 300,
        re_render: bool = False,
        **render_kwargs,
    ) -> Path:
        """
        Export the session graph to a publication-quality file.

        Supported formats (detected from extension):
          .svg   — Vector, lossless, editable in Inkscape/Illustrator
          .pdf   — Vector, lossless, direct journal submission
          .png   — Raster at specified dpi (min 300 recommended for journals)
          .tiff  — Raster at specified dpi (TIFF required by some journals)

        Parameters
        ----------
        path:
            Output file path including extension.
        dpi:
            Resolution for raster formats. 300 = print quality. 600 = ultra-high.
        re_render:
            If True, re-renders the figure before saving even if a cached figure
            exists. Pass extra render kwargs alongside this.
        **render_kwargs:
            Passed to render() if re_render=True or no figure is cached.

        Returns
        -------
        Path to the saved file.

        Raises
        ------
        EmptySessionError
            If no metabolites have been added yet.
        """
        if self._last_fig is None or re_render:
            self.render(**render_kwargs)

        return save_figure(self._last_fig, path, dpi=dpi)

    # ------------------------------------------------------------------
    # Introspection API
    # ------------------------------------------------------------------

    def build(self) -> PathwayGraph:
        """
        Return the current graph as an immutable PathwayGraph snapshot.

        Useful for programmatic access to the graph data without rendering.
        The PathwayGraph exposes the underlying nx.DiGraph, all node/edge
        schema objects, and summary statistics.

        Raises
        ------
        EmptySessionError
        """
        if self._state.is_empty:
            raise EmptySessionError()
        return self._state.build_pathway_graph()

    def summary(self) -> str:
        """Return a human-readable summary of the current session state."""
        if self._state.is_empty:
            return "PyPyrusMap Map session: empty. Call session.add('metabolite_id') to begin."
        pg = self._state.build_pathway_graph()
        anchors = ", ".join(pg.anchor_ids)
        return (
            f"PyPyrusMap Map session\n"
            f"  Model   : {self._model.id}\n"
            f"  Anchors : {anchors}\n"
            f"  Graph   : {pg.summary()}\n"
            f"  Flux    : {'yes (solution attached)' if self._solution is not None else 'no solution provided'}"
        )

    def __repr__(self) -> str:
        return self.summary()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def model(self) -> "cobra.Model":
        """The COBRApy model attached to this session (read-only)."""
        return self._model

    @property
    def solution(self) -> Optional["cobra.Solution"]:
        """The FBA solution attached to this session, or None."""
        return self._solution

    @property
    def anchor_ids(self) -> list[str]:
        """Ordered list of metabolite IDs added via session.add()."""
        return list(self._state._anchor_ids)

    @property
    def currency_ids(self) -> frozenset[str]:
        """The currency metabolite ID set currently in use."""
        return self._currency_ids if self._currency_ids is not None else DEFAULT_CURRENCY_IDS
