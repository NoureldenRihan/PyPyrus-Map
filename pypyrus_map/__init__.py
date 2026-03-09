"""
PyPyrusMap
===========
A reusable Python tool for exploring and visualizing metabolic neighborhoods
in COBRApy genome-scale models.

Developed by Nourelden Rihan.

Quick start
-----------
>>> import cobra
>>> from pypyrus_map import PyPyrusMap
>>>
>>> model = cobra.io.load_json_model("iML1515.json")
>>> solution = model.optimize()
>>>
>>> session = PyPyrusMap(model, solution=solution)
>>> session.add("ggpp_c").add("phyto_c")
>>> fig = session.render(layout="dot", title="Carotenoid Pathway")
>>> session.export("carotenoid_pathway.svg")
>>> session.export("carotenoid_pathway.pdf")
"""

from .session import PyPyrusMap
from .schema import PathwayGraph, MetaboliteNode, ReactionNode, MetabolicEdge
from .exceptions import (
    PyPyrusMapError,
    MetaboliteNotFoundError,
    ReactionNotFoundError,
    NoPathFoundError,
    EmptySessionError,
)
from .query import DEFAULT_CURRENCY_IDS

__all__ = [
    # Primary API
    "PyPyrusMap",
    # Schema types (for type annotations in downstream code)
    "PathwayGraph",
    "MetaboliteNode",
    "ReactionNode",
    "MetabolicEdge",
    # Exceptions
    "PyPyrusMapError",
    "MetaboliteNotFoundError",
    "ReactionNotFoundError",
    "NoPathFoundError",
    "EmptySessionError",
    # Utilities
    "DEFAULT_CURRENCY_IDS",
]

__version__ = "0.1.1"
__author__ = "Nourelden Rihan"
__license__ = "MIT"

