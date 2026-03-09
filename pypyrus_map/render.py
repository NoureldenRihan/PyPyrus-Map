"""
pypyrus_map.render
==================
Publication-quality rendering of PathwayGraph objects.

Changelog
---------
v2  : Vectorised scatter/annotate rendering.
v3  : Round 1 — label modes, arrows, inline flux labels, auto-canvas, orientation.
v3.1: Defaults → landscape, show_flux_labels=False.
v4  : Round 2 — scatter + plain text nodes, len()-based sizing.
v5  : Round 2 polish —
      - Separate sizing constants for metabolites (1-line) vs reactions (2-line)
      - Width-only expansion: height fixed per node type, width follows label len()
      - Node sizes passed to Graphviz dot as width/height attributes so the
        layout engine spaces nodes to avoid overlap — fixes canvas cramping,
        margins, and arrow visibility in one shot
      - shrinkB = half_node_width_pts + 14pt buffer (guaranteed arrowhead clearance)
      - Spring layout fallback unaffected (acceptable per design decision)

Node sizing strategy
--------------------
Height is FIXED per node type (single/double line). Width scales with label.

    met_w  = max_chars * MET_CHAR_W + 2 * MET_PAD_W   (variable)
    met_h  = MET_LINE_H + 2 * MET_PAD_H               (fixed, 1 line)

    rxn_w  = max_chars * RXN_CHAR_W + 2 * RXN_PAD_W   (variable)
    rxn_h  = 2 * RXN_LINE_H + 2 * RXN_PAD_H           (fixed, 2 lines)

    scatter s = w * h   (area in display points²)

Graphviz sizing
---------------
We convert pt → inches (÷ 72) and pass width/height + fixedsize=true
as node attributes. Graphviz dot then spaces nodes to prevent overlap,
so the rendered canvas naturally fits the content.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional, Literal

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import networkx as nx

from .schema import PathwayGraph, MetaboliteNode, ReactionNode
from .graph import NODE_TYPE_KEY, NODE_DATA_KEY, TYPE_METABOLITE, TYPE_REACTION


# ---------------------------------------------------------------------------
# Sizing constants — metabolites and reactions tuned independently
# ---------------------------------------------------------------------------

# Metabolites — single line, width-only expansion
MET_CHAR_W =  9.0   # pt per character
MET_PAD_W  = 10.0   # horizontal padding each side (pt)
MET_LINE_H = 14.0   # fixed single-line height (pt)
MET_PAD_H  =  7.0   # vertical padding each side (pt)

# Reactions — two lines (ID + genes), width-only expansion, smaller amplifier
RXN_CHAR_W =  7.5   # pt per character (lower than met to keep boxes tighter)
RXN_PAD_W  =  9.0   # horizontal padding each side (pt)
RXN_LINE_H = 11.0   # height per line (pt)
RXN_PAD_H  =  6.0   # vertical padding each side (pt)
RXN_LINES  =  2     # always 2 lines: reaction ID + gene names

# Arrow clearance buffer added on top of computed half-width
_ARROW_BUFFER = 5.0    # pt safety buffer added on top of computed visual half-width

# Graphviz spacing — extra margin around each node (inches) so nodes breathe
_GV_NODE_MARGIN = 0.15   # inches added to each side when telling dot node size


# ---------------------------------------------------------------------------
# Label helpers
# ---------------------------------------------------------------------------

_TRUNCATE_LEN = 20


def _metabolite_label(
    mn: MetaboliteNode,
    label_mode: Literal["id", "truncate", "full"],
) -> str:
    if label_mode == "id":
        return mn.id
    if label_mode == "full":
        return mn.name if mn.name else mn.id
    name = mn.name if mn.name else mn.id
    return name if len(name) <= _TRUNCATE_LEN else name[:_TRUNCATE_LEN - 1] + "…"


def _reaction_label(rn: ReactionNode) -> str:
    gene_label = ", ".join(rn.genes[:2])
    if len(rn.genes) > 2:
        gene_label += f" +{len(rn.genes) - 2}"
    return f"{rn.id}\n{gene_label}" if gene_label else rn.id


def _node_label(node_id: str, g: nx.DiGraph, label_mode: str) -> str:
    data      = g.nodes[node_id][NODE_DATA_KEY]
    node_type = g.nodes[node_id][NODE_TYPE_KEY]
    if node_type == TYPE_METABOLITE:
        return _metabolite_label(data, label_mode)
    return _reaction_label(data)


# ---------------------------------------------------------------------------
# Node size computation — pure Python, no Matplotlib
# ---------------------------------------------------------------------------

def _met_dimensions(label: str) -> tuple[float, float]:
    """Return (width_pt, height_pt) for a metabolite node."""
    max_chars = max(len(l) for l in label.split("\n"))
    w = max_chars * MET_CHAR_W + 2 * MET_PAD_W
    h = MET_LINE_H + 2 * MET_PAD_H          # height fixed — 1 line always
    return w, h


def _rxn_dimensions(label: str) -> tuple[float, float]:
    """Return (width_pt, height_pt) for a reaction node."""
    max_chars = max(len(l) for l in label.split("\n"))
    w = max_chars * RXN_CHAR_W + 2 * RXN_PAD_W
    h = RXN_LINES * RXN_LINE_H + 2 * RXN_PAD_H   # height fixed — 2 lines always
    return w, h


def _node_dimensions(node_id: str, g: nx.DiGraph, label_mode: str) -> tuple[float, float]:
    """Return (width_pt, height_pt) for any node."""
    node_type = g.nodes[node_id][NODE_TYPE_KEY]
    label     = _node_label(node_id, g, label_mode)
    if node_type == TYPE_METABOLITE:
        return _met_dimensions(label)
    return _rxn_dimensions(label)


def _scatter_s(w_pt: float, h_pt: float) -> float:
    """Convert box dimensions (pt) to scatter marker area (pt²)."""
    return w_pt * h_pt


import math as _math


def _shrink_b(w_pt: float, h_pt: float) -> float:
    """
    Compute shrinkB (pt) so arrowheads land just outside the scatter marker.

    scatter renders at area s = w * h, so the actual visual marker
    half-width is sqrt(w*h)/2 — NOT w/2. The old w/2 formula was
    2-2.5x too large, making arrows vanish on short connections.
    """
    return _math.sqrt(w_pt * h_pt) / 2 + _ARROW_BUFFER


# ---------------------------------------------------------------------------
# Layout — with Graphviz node-size awareness for dot
# ---------------------------------------------------------------------------

def _compute_layout(
    g: nx.DiGraph,
    layout: Literal["dot", "spring", "circular"],
    orientation: Literal["portrait", "landscape"],
    label_mode: str,
) -> dict[str, tuple[float, float]]:
    """
    Compute 2D node positions.

    For 'dot' layout: tries pydot (pure Python, works on all platforms).
    Falls back to spring layout if Graphviz is not installed.

    For spring/circular: standard NetworkX layout.
    """
    rankdir = "TB" if orientation == "portrait" else "LR"

    if layout == "dot":
        try:
            pos = _dot_layout_pydot(g, rankdir, label_mode)
            if pos:
                return pos
        except Exception:
            pass

        warnings.warn(
            "Graphviz dot layout unavailable — falling back to spring layout.\n"
            "For publication-quality figures:\n"
            "  1. Install Graphviz from https://graphviz.org/download/\n"
            "     (tick 'Add to PATH' on Windows)\n"
            "  2. pip install pydot",
            UserWarning, stacklevel=3,
        )

    if layout in ("dot", "spring"):
        return nx.spring_layout(g, seed=42, k=3.5)
    if layout == "circular":
        return nx.circular_layout(g)
    return nx.spring_layout(g, seed=42, k=3.5)


def _dot_layout_pydot(
    g: nx.DiGraph,
    rankdir: str,
    label_mode: str,
) -> dict[str, tuple[float, float]]:
    """
    Dot layout via pydot — pure Python, works on all platforms.
    Requires: pip install pydot  +  Graphviz in system PATH.
    Node sizes passed as width/height so dot spaces nodes correctly.
    """
    import pydot

    PT_TO_IN = 1.0 / 72.0

    graph = pydot.Dot(
        graph_type="digraph",
        rankdir=rankdir,
        splines="spline",
        nodesep="0.6",
        ranksep="0.8",
    )

    for node_id in g.nodes:
        w_pt, h_pt = _node_dimensions(node_id, g, label_mode)
        w_in = w_pt * PT_TO_IN + _GV_NODE_MARGIN
        h_in = h_pt * PT_TO_IN + _GV_NODE_MARGIN
        graph.add_node(pydot.Node(
            f'"{node_id}"',
            shape="box",
            fixedsize="true",
            width=f"{w_in:.4f}",
            height=f"{h_in:.4f}",
        ))

    for src, tgt in g.edges():
        graph.add_edge(pydot.Edge(f'"{src}"', f'"{tgt}"'))

    dot_output = graph.create(prog="dot", format="dot").decode("utf-8")
    parsed = pydot.graph_from_dot_data(dot_output)
    if not parsed:
        return {}

    pos = {}
    for node in parsed[0].get_nodes():
        node_id = node.get_name().strip('"')
        if node_id in ("node", "graph", "edge", ""):
            continue
        pos_attr = node.get_pos()
        if pos_attr:
            coords = pos_attr.strip('"')
            x, y = coords.split(",")
            pos[node_id] = (float(x), float(y))
    return pos


# ---------------------------------------------------------------------------
# Colours — Okabe-Ito palette (colorblind-safe, high contrast)
# ---------------------------------------------------------------------------

_PALETTE = {
    "anchor_metabolite":   "#005AB5",
    "neighbor_metabolite": "#1AB2FF",
    "reaction":            "#F0A500",
    "edge_consumes":       "#D62598",
    "edge_produces":       "#00A875",
    "background":          "#F8F9FA",
    "text_dark":           "#1A1A1A",
    "text_light":          "#FFFFFF",
}


# ---------------------------------------------------------------------------
# Auto figure sizing
# ---------------------------------------------------------------------------

_MAX_DIM = 40


def _auto_figsize(
    g: nx.DiGraph,
    orientation: Literal["portrait", "landscape"],
    user_figsize: Optional[tuple[float, float]],
) -> tuple[float, float]:
    if user_figsize is not None:
        w, h = user_figsize
        return min(w, _MAX_DIM), min(h, _MAX_DIM)
    n       = g.number_of_nodes()
    primary = min(max(14.0, n * 0.8), _MAX_DIM)
    cross   = min(20.0, _MAX_DIM)
    return (cross, primary) if orientation == "portrait" else (primary, cross)


# ---------------------------------------------------------------------------
# Core render entry point
# ---------------------------------------------------------------------------

def render_pathway_graph(
    pg: PathwayGraph,
    layout: Literal["dot", "spring", "circular"] = "dot",
    orientation: Literal["portrait", "landscape"] = "landscape",
    figsize: Optional[tuple[float, float]] = None,
    title: Optional[str] = None,
    label_mode: Literal["id", "truncate", "full"] = "id",
    show_stoichiometry: bool = False,
    show_flux_labels: bool = False,
    legend_loc: str = "upper left",
    dpi: int = 110,
    font_family: str = "DejaVu Sans",
) -> plt.Figure:
    """
    Render a PathwayGraph as a publication-quality Matplotlib figure.

    Parameters
    ----------
    pg : PathwayGraph
    layout : 'dot' | 'spring' | 'circular'
    orientation : 'portrait' | 'landscape'  (default: landscape)
    figsize : (w, h) inches — auto-sized from node count if None
    title : auto-generated from anchor IDs if None
    label_mode : 'id' (default) | 'truncate' | 'full'
    show_stoichiometry : stoichiometric coefficients on edge midpoints
    show_flux_labels : flux value (4 dp) on edge midpoints (default: off)
    legend_loc : any Matplotlib legend location string, e.g.:
        'upper left' (default), 'upper right', 'lower left', 'lower right',
        'lower center', 'upper center', 'right', 'center left', 'center right'
    dpi : screen resolution — SVG/PDF always lossless
    font_family : 'Arial' matches Nature/Cell/Science guidelines
    """
    plt.close("all")

    matplotlib.rcParams["font.family"]  = font_family
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["svg.fonttype"] = "none"

    g = pg.nx_graph
    if g.number_of_nodes() == 0:
        raise ValueError("Graph is empty — call session.add() before rendering.")

    pos  = _compute_layout(g, layout, orientation, label_mode)
    w, h = _auto_figsize(g, orientation, figsize)

    fig, ax = plt.subplots(figsize=(w, h), dpi=dpi)
    fig.patch.set_facecolor(_PALETTE["background"])
    ax.set_facecolor(_PALETTE["background"])
    ax.axis("off")

    met_nodes = [n for n in g.nodes if g.nodes[n][NODE_TYPE_KEY] == TYPE_METABOLITE]
    rxn_nodes = [n for n in g.nodes if g.nodes[n][NODE_TYPE_KEY] == TYPE_REACTION]

    # Draw order: edges → markers → labels → legend
    _draw_edges(ax, g, pos, pg, label_mode, show_stoichiometry, show_flux_labels)
    _draw_metabolite_markers(ax, g, pos, met_nodes, label_mode)
    _draw_reaction_markers(ax, g, pos, rxn_nodes)
    _draw_metabolite_labels(ax, g, pos, met_nodes, label_mode)
    _draw_reaction_labels(ax, g, pos, rxn_nodes)
    _draw_legend(ax, legend_loc)

    anchor_ids    = [g.nodes[a][NODE_DATA_KEY].id for a in pg.anchor_ids if a in g]
    default_title = "Pathway — " + ", ".join(anchor_ids) if anchor_ids else "Pathway Graph"
    ax.set_title(
        title or default_title,
        fontsize=13, fontweight="bold", pad=14, color=_PALETTE["text_dark"],
    )
    fig.text(
        0.99, 0.005,
        "Generated by PyPyrusMap · Nourelden Rihan",
        ha="right", va="bottom", fontsize=7, color="#AAAAAA", style="italic",
    )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.02)
    return fig


# ---------------------------------------------------------------------------
# Node markers
# ---------------------------------------------------------------------------

def _draw_metabolite_markers(
    ax: plt.Axes,
    g: nx.DiGraph,
    pos: dict,
    met_nodes: list[str],
    label_mode: str,
) -> None:
    if not met_nodes:
        return

    for is_anchor_pass in (False, True):
        xs, ys, sizes = [], [], []
        for nid in met_nodes:
            mn: MetaboliteNode = g.nodes[nid][NODE_DATA_KEY]
            if mn.is_anchor != is_anchor_pass:
                continue
            label  = _metabolite_label(mn, label_mode)
            w, h   = _met_dimensions(label)
            xs.append(pos[nid][0])
            ys.append(pos[nid][1])
            sizes.append(_scatter_s(w, h))

        if not xs:
            continue

        ax.scatter(
            xs, ys,
            s=sizes,
            marker="s",
            c=_PALETTE["anchor_metabolite"] if is_anchor_pass else _PALETTE["neighbor_metabolite"],
            edgecolors="#1A1A1A",
            linewidths=2.2 if is_anchor_pass else 1.2,
            zorder=3,
        )


def _draw_reaction_markers(
    ax: plt.Axes,
    g: nx.DiGraph,
    pos: dict,
    rxn_nodes: list[str],
) -> None:
    if not rxn_nodes:
        return

    xs, ys, sizes = [], [], []
    for nid in rxn_nodes:
        rn: ReactionNode = g.nodes[nid][NODE_DATA_KEY]
        w, h = _rxn_dimensions(_reaction_label(rn))
        xs.append(pos[nid][0])
        ys.append(pos[nid][1])
        sizes.append(_scatter_s(w, h))

    ax.scatter(
        xs, ys,
        s=sizes,
        marker="D",
        c=_PALETTE["reaction"],
        edgecolors="#1A1A1A",
        linewidths=1.2,
        zorder=3,
    )


# ---------------------------------------------------------------------------
# Node labels — plain ax.text, no bbox kwarg
# ---------------------------------------------------------------------------

def _draw_metabolite_labels(
    ax: plt.Axes,
    g: nx.DiGraph,
    pos: dict,
    met_nodes: list[str],
    label_mode: str,
) -> None:
    for nid in met_nodes:
        mn: MetaboliteNode = g.nodes[nid][NODE_DATA_KEY]
        x, y  = pos[nid]
        label = _metabolite_label(mn, label_mode)
        ax.text(
            x, y, label,
            ha="center", va="center",
            fontsize=9,
            fontweight="bold" if mn.is_anchor else "normal",
            color=_PALETTE["text_light"],
            zorder=5,
            path_effects=[pe.withStroke(linewidth=1.5, foreground="#1A1A1A")],
        )


def _draw_reaction_labels(
    ax: plt.Axes,
    g: nx.DiGraph,
    pos: dict,
    rxn_nodes: list[str],
) -> None:
    for nid in rxn_nodes:
        rn: ReactionNode = g.nodes[nid][NODE_DATA_KEY]
        x, y = pos[nid]
        ax.text(
            x, y, _reaction_label(rn),
            ha="center", va="center",
            fontsize=8,
            color=_PALETTE["text_dark"],
            zorder=5,
            multialignment="center",
        )


# ---------------------------------------------------------------------------
# Edge drawing
# ---------------------------------------------------------------------------

def _draw_edges(
    ax: plt.Axes,
    g: nx.DiGraph,
    pos: dict,
    pg: PathwayGraph,
    label_mode: str,
    show_stoichiometry: bool,
    show_flux_labels: bool,
) -> None:
    for src, tgt, edata in g.edges(data=True):
        edge     = edata.get(NODE_DATA_KEY)
        tgt_type = g.nodes[tgt][NODE_TYPE_KEY]
        rxn_id   = tgt if tgt_type == TYPE_REACTION else src
        rn: ReactionNode = g.nodes[rxn_id][NODE_DATA_KEY]

        # Color always direction-based — blocked never overrides
        color = (
            _PALETTE["edge_produces"]
            if (edge is not None and edge.direction == "produces")
            else _PALETTE["edge_consumes"]
        )

        linestyle   = "-"
        alpha       = 0.90
        arrow_style = "<|-|>" if rn.is_reversible else "-|>"
        rad         = "0.14"  if g.has_edge(tgt, src) else "0.0"

        # shrinkB from target node actual visual size — sqrt(w*h)/2 + buffer
        tgt_label = _node_label(tgt, g, label_mode)
        tgt_type_ = g.nodes[tgt][NODE_TYPE_KEY]
        if tgt_type_ == TYPE_METABOLITE:
            tgt_w, tgt_h = _met_dimensions(tgt_label)
        else:
            tgt_w, tgt_h = _rxn_dimensions(tgt_label)

        # shrinkA from source node — same formula, fixes back-arrowhead on reversible
        src_label = _node_label(src, g, label_mode)
        src_type_ = g.nodes[src][NODE_TYPE_KEY]
        if src_type_ == TYPE_METABOLITE:
            src_w, src_h = _met_dimensions(src_label)
        else:
            src_w, src_h = _rxn_dimensions(src_label)

        shrink_b = _shrink_b(tgt_w, tgt_h)
        shrink_a = _shrink_b(src_w, src_h)   # same formula, applied to source

        base_props = dict(
            arrowstyle=arrow_style,
            linestyle=linestyle,
            alpha=alpha,
            mutation_scale=15,
            shrinkA=shrink_a,
            shrinkB=shrink_b,
            connectionstyle=f"arc3,rad={rad}",
        )

        # Blocked arrows: draw thick black outline first, then color on top
        if pg.has_flux and rn.is_blocked:
            ax.annotate(
                "",
                xy=pos[tgt], xytext=pos[src],
                arrowprops=dict(**base_props, color="#1A1A1A", lw=5.5),
                zorder=2,
            )

        ax.annotate(
            "",
            xy=pos[tgt], xytext=pos[src],
            arrowprops=dict(**base_props, color=color, lw=2.0),
            zorder=3,
        )

        # Mid-edge labels
        mid: list[str] = []
        if show_flux_labels and pg.has_flux and rn.flux is not None:
            mid.append(f"{rn.flux:.4f}")
        if show_stoichiometry and edge is not None and edge.stoichiometry != 1.0:
            mid.append(f"×{edge.stoichiometry:.0f}")

        if mid:
            sx, sy = pos[src]
            tx, ty = pos[tgt]
            ax.text(
                (sx + tx) / 2, (sy + ty) / 2,
                "  ".join(mid),
                fontsize=7, color=color,
                ha="center", va="center", zorder=6,
                path_effects=[pe.withStroke(linewidth=4, foreground="white")],
            )


# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------

def _draw_legend(ax: plt.Axes, loc: str = "upper left") -> None:
    handles = [
        mpatches.Patch(fc=_PALETTE["anchor_metabolite"],
                       ec="#1A1A1A", lw=2.0, label="Focal metabolite"),
        mpatches.Patch(fc=_PALETTE["neighbor_metabolite"],
                       ec="#1A1A1A", lw=1.2, label="Neighbor metabolite"),
        mpatches.Patch(fc=_PALETTE["reaction"],
                       ec="#1A1A1A", lw=1.2, label="Reaction"),
        plt.Line2D([0], [0], color=_PALETTE["edge_produces"],
                   lw=2.5, label="Produces →"),
        plt.Line2D([0], [0], color=_PALETTE["edge_consumes"],
                   lw=2.5, label="Consumes →"),
        plt.Line2D([0], [0], color=_PALETTE["edge_produces"],
                   lw=2.5, label="Produces (blocked)",
                   path_effects=[pe.withStroke(linewidth=5.5, foreground="#1A1A1A")]),
        plt.Line2D([0], [0], color=_PALETTE["edge_consumes"],
                   lw=2.5, label="Consumes (blocked)",
                   path_effects=[pe.withStroke(linewidth=5.5, foreground="#1A1A1A")]),
    ]
    ax.legend(
        handles=handles, loc=loc, fontsize=8,
        framealpha=0.95, edgecolor="#CCCCCC",
        facecolor="white", borderpad=1.0, handlelength=2.2,
    )


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def save_figure(fig: plt.Figure, path: str | Path, dpi: int = 300) -> Path:
    """Save to disk. SVG/PDF lossless vector. PNG/TIFF use specified dpi."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ext = path.suffix.lower()
    kw: dict = {"bbox_inches": "tight", "facecolor": fig.get_facecolor()}
    if ext not in {".svg", ".pdf"}:
        kw["dpi"] = dpi
    fig.savefig(path, **kw)
    return path
