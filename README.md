# PyPyrus Map

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18870982.svg)](https://doi.org/10.5281/zenodo.18870982)
[![PyPI](https://img.shields.io/pypi/v/pypyrus-map)](https://pypi.org/project/pypyrus-map/)

**Metabolic neighborhood visualization for COBRApy genome-scale models.**

PyPyrus Map lets you instantly visualize the full reaction neighborhood of any metabolite in a COBRApy-compatible genome-scale model — including heterologous and non-native pathways that curated databases like Escher cannot represent. Output is publication-quality vector figures (SVG, PDF) suitable for direct journal submission.

---

## Why PyPyrus Map exists

PyPyrus Map was built out of a concrete research problem during [Project Menhed](https://github.com/NoureldenRihan), an *in silico* metabolic engineering effort to optimize lycopene production across multiple chassis organisms. When working with heterologous carotenoid biosynthesis pathways, standard tools like Escher were limited to pre-drawn maps — they could not show the neighborhood of a non-native intermediate like GGPP or Phytoene, making it difficult to quickly identify overexpression and knockout targets.

The need was simple: given a metabolite of interest, show every reaction that produces or consumes it, who the other substrates and products are, and — when an FBA solution is available — what the flux through each reaction looks like. PyPyrus Map was written to answer that question in one function call.

---

## Installation

**Step 1: Install PyPyrus Map:**

```bash
pip install pypyrus-map
```

**Step 2: Install Graphviz** (required for `layout="dot"`):
```bash
# Ubuntu / Debian 
apt-get install graphviz

# macOS
brew install graphviz

# Windows
# Download and run the installer from https://graphviz.org/download/
# During installation, tick "Add Graphviz to system PATH"
# Then restart your terminal
```

---

## Quick start

```python
import cobra
from pypyrus_map import PyPyrusMap

# Load any COBRApy-compatible genome-scale model
model = cobra.io.load_model("iML1515")
solution = model.optimize()

# Create a session and explore a metabolite neighborhood
session = PyPyrusMap(model, solution=solution)
session.add("ipdp_c")

fig = session.render(title="IPP Neighborhood")
fig.show()

# Export publication-quality vector figures
session.export("ipdp_neighborhood.svg")
session.export("ipdp_neighborhood.pdf")
```

---

## Chaining metabolites

PyPyrus Map is designed for pathway chains. Call `.add()` multiple times to build up a connected pathway graph — shared nodes are never duplicated.

```python
session = PyPyrusMap(model, solution=solution)

# Build a carotenoid biosynthesis chain
session.add("ipdp_c")   # IPP neighborhood
session.add("frdp_c")   # FPP neighborhood 

fig = session.render(
    title="Isoprenoid Pathway — IPP → FPP",
)
session.export("isoprenoid_chain.svg")
```

When you call `session.add("frdp_c")`, PyPyrus Map queries FPP's full depth-1 neighborhood and merges it into the existing graph. Any node already present is promoted to anchor status rather than duplicated.

---

## Visualization options

```python
fig = session.render(
    layout="dot",              # 'dot' (default) | 'spring' | 'circular'
    orientation="landscape",   # 'landscape' (default) | 'portrait'
    figsize=None,              # auto-sized from node count; override with (width, height)
    title=None,                # auto-generated from anchor IDs if None
    label_mode="id",           # 'id' (default) | 'truncate' | 'full'
    show_flux_labels=False,    # show flux values on edge midpoints
    show_stoichiometry=False,  # show stoichiometric coefficients on edges
    legend_loc="upper left",   # any Matplotlib legend location string
    dpi=110,                   # screen resolution; SVG/PDF always lossless
    font_family="DejaVu Sans", # 'Arial' matches Nature/Cell/Science guidelines
)
```

### Label modes

| Mode | Shows | Best for |
|---|---|---|
| `"id"` | BiGG ID (e.g. `ggpp_c`) | Default — compact, unambiguous |
| `"truncate"` | Human name, max 20 chars | Readable without overlap |
| `"full"` | Full human name | Presentations, wide figures |

### Visual encoding

| Element | Meaning |
|---|---|
| **Deep blue square** | Focal metabolite (anchor) |
| **Sky blue square** | Neighbor metabolite |
| **Amber diamond** | Reaction |
| **Teal-green arrow** | Produces (reaction → metabolite) |
| **Magenta-pink arrow** | Consumes (metabolite → reaction) |
| **Double-headed arrow** | Reversible reaction |
| **Colored arrow with black outline** | Edge through a blocked reaction |

---

## Flux overlay

When a `cobra.Solution` is passed at session construction, flux values are attached to each reaction node. Blocked reactions (flux ≈ 0) are visually distinguished by a thick black border on the reaction diamond and a black outline on their connecting arrows.

```python
solution = model.optimize()
session = PyPyrusMap(model, solution=solution)
session.add("ggpp_c")

fig = session.render(show_flux_labels=True)
```

---

## Currency metabolite suppression

High-connectivity cofactors (ATP, NADH, H₂O, CoA, Pi, etc.) connect nearly every reaction in the model to every other. By default, PyPyrus Map suppresses these from neighbor nodes to keep figures clean. The focal metabolite is always shown regardless.

```python
# Default: currency suppressed
session = PyPyrusMap(model)

# Show everything including currency metabolites
session = PyPyrusMap(model, suppress_currency=False)

# Custom currency list
my_currency = frozenset({"h_c", "h_e", "atp_c", "adp_c"})
session = PyPyrusMap(model, custom_currency_ids=my_currency)
```

The default currency ID list is available as `pypyrus_map.DEFAULT_CURRENCY_IDS` and follows the conventions of Huss & Holme (2007) and KEGG RPAIR.

---

## Error handling

PyPyrus Map uses strict BiGG IDs — no fuzzy matching. If an ID is not found, a clear error is raised with prefix-matched candidate suggestions:

```python
from pypyrus_map import MetaboliteNotFoundError

try:
    session.add("ggpp_m")   # wrong compartment
except MetaboliteNotFoundError as e:
    print(e)
# MetaboliteNotFoundError: 'ggpp_m' not found in model.
# Did you mean: ggpp_c?
```

---

## Session API

### `PyPyrusMap(model, solution=None, suppress_currency=True, custom_currency_ids=None)`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `model` | `cobra.Model` | required | Any COBRApy-compatible genome-scale model |
| `solution` | `cobra.Solution` | `None` | FBA solution — enables flux overlay |
| `suppress_currency` | `bool` | `True` | Suppress high-connectivity cofactors from neighbor nodes |
| `custom_currency_ids` | `frozenset[str]` | `None` | Replace the default currency list entirely |

### Session methods

| Method | Returns | Description |
|---|---|---|
| `session.add(metabolite_id)` | `self` | Add a metabolite neighborhood. Chainable. |
| `session.remove(metabolite_id)` | `self` | Remove a metabolite and its orphaned reactions. Chainable. |
| `session.clear()` | `self` | Reset the session graph entirely. |
| `session.render(**kwargs)` | `Figure` | Render and return a Matplotlib figure. |
| `session.export(path, dpi=300)` | `Path` | Save figure to SVG, PDF, PNG, or TIFF. |
| `session.build()` | `PathwayGraph` | Return the raw graph object for custom processing. |
| `session.summary()` | `str` | Print a text summary of current session state. |

---

## Development

```bash
git clone https://github.com/NoureldenRihan/PyPyrus-Map.git
cd PyPyrus-Map
pip install -e ".[dev]"
pytest tests/ -v
```

---

## Citing PyPyrus Map

If PyPyrus Map contributes to a published work, please cite:

> Rihan, N. *PyPyrus Map: Metabolic neighborhood visualization for COBRApy genome-scale models.* (2026). DOI: 10.5281/zenodo.18870982

```bibtex
@software{rihan_2026_pypyrus_map,
  author    = {Rihan, Nourelden},
  title     = {PyPyrus Map: Metabolic neighborhood visualization for COBRApy genome-scale models},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.18870982},
  url       = {https://doi.org/10.5281/zenodo.18870982}
}
```

---

## Acknowledgements

PyPyrus Map was developed with the assistance of Claude (Anthropic), a large language model used for code generation and implementation. All research direction, architectural decisions, feature design, and validation were conceived and directed by the author.

PyPyrus Map was built as part of Project Menhed, an *in silico* metabolic engineering initiative targeting lycopene biosynthesis optimization.

---

## License

MIT License — Copyright (c) 2026 Nourelden Rihan. See [LICENSE](LICENSE) for details.
