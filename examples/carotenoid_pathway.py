"""
PyPyrus Map — Example: Isoprenoid Pathway Neighborhood
=======================================================
Demonstrates PyPyrus Map using two metabolites from the native
isoprenoid biosynthesis pathway in iML1515 (E. coli).

ipdp_c  — Isopentenyl diphosphate (IPP)
frdp_c  — Farnesyl diphosphate (FPP)

Both are present in iML1515 and do not require any heterologous
pathway additions. This example will work on any fresh iML1515 install.

Run
---
    pip install pypyrus-map
    python examples/isoprenoid_pathway.py
"""

import cobra
from pypyrus_map import PyPyrusMap
from pypyrus_map.exceptions import MetaboliteNotFoundError

# ── 1. Load model ────────────────────────────────────────────────────────────

print("Loading iML1515...")
try:
    model = cobra.io.load_json_model("iML1515.json")
except FileNotFoundError:
    raise SystemExit(
        "iML1515.json not found in current directory.\n"
        "Download it from: https://github.com/SBRG/iML1515\n"
        "  or via: cobra.io.load_model('iML1515')"
    )

print(f"Model loaded: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

# ── 2. Run FBA ───────────────────────────────────────────────────────────────

solution = model.optimize()
print(f"FBA objective: {solution.objective_value:.4f}")

# ── 3. Single metabolite — IPP neighborhood ──────────────────────────────────

print("\n--- IPP neighborhood (ipdp_c) ---")
session = PyPyrusMap(model, solution=solution)
session.add("ipdp_c")
print(session.summary())

fig = session.render(
    title="IPP (Isopentenyl Diphosphate) Neighborhood",
    orientation="landscape",
    label_mode="id",
)
session.export("ipdp_neighborhood.svg")
session.export("ipdp_neighborhood.pdf")
print("Saved: ipdp_neighborhood.svg / .pdf")

# ── 4. Chained — IPP + FPP pathway ──────────────────────────────────────────

print("\n--- IPP + FPP chained neighborhood ---")
session2 = PyPyrusMap(model, solution=solution)
session2.add("ipdp_c").add("frdp_c")
print(session2.summary())

fig2 = session2.render(
    title="Isoprenoid Pathway — IPP → FPP",
    orientation="portrait",
    label_mode="id",
)
session2.export("isoprenoid_chain.svg")
session2.export("isoprenoid_chain.pdf")
print("Saved: isoprenoid_chain.svg / .pdf")

# ── 5. Full names mode ───────────────────────────────────────────────────────

fig3 = session2.render(
    title="Isoprenoid Pathway — Full Names",
    orientation="portrait",
    label_mode="truncate",
)
session2.export("isoprenoid_chain_names.svg")
print("Saved: isoprenoid_chain_names.svg")

# ── 6. Error handling demo ───────────────────────────────────────────────────

print("\n--- Error handling ---")
try:
    session2.add("ipdp_m")   # wrong compartment — should be ipdp_c
except MetaboliteNotFoundError as e:
    print(f"Caught expected error: {e}")

print("\nDone.")

