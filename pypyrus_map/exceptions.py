"""
pypyrus_map.exceptions
======================
All custom exceptions raised by PyPyrusMap.
"""


class PyPyrusMapError(Exception):
    """Base exception for all PyPyrusMap errors."""


class MetaboliteNotFoundError(PyPyrusMapError):
    """
    Raised when a metabolite ID is not found in the provided COBRApy model.

    Carries the attempted ID and a list of close candidates (if any were
    found via prefix matching) to help the user self-correct.
    """

    def __init__(self, metabolite_id: str, model_id: str, candidates: list[str] = None):
        self.metabolite_id = metabolite_id
        self.model_id = model_id
        self.candidates = candidates or []

        hint = ""
        if self.candidates:
            formatted = ", ".join(f"'{c}'" for c in self.candidates[:5])
            hint = f"\n  Did you mean one of: {formatted}?"

        super().__init__(
            f"Metabolite '{metabolite_id}' not found in model '{model_id}'.{hint}\n"
            f"  Tip: Use exact BiGG IDs including compartment suffix, e.g. 'glc__D_c'."
        )


class ReactionNotFoundError(PyPyrusMapError):
    """Raised when a reaction ID is not found in the provided COBRApy model."""

    def __init__(self, reaction_id: str, model_id: str):
        self.reaction_id = reaction_id
        self.model_id = model_id
        super().__init__(
            f"Reaction '{reaction_id}' not found in model '{model_id}'."
        )


class NoPathFoundError(PyPyrusMapError):
    """
    Raised by trace_pathway() when no connected path exists between
    the source and target metabolites within the current graph.
    """

    def __init__(self, source_id: str, target_id: str):
        self.source_id = source_id
        self.target_id = target_id
        super().__init__(
            f"No pathway found between '{source_id}' and '{target_id}'.\n"
            f"  The metabolites may not be connected in this model, or the path "
            f"may require depth > 1. Try adding intermediate metabolites manually "
            f"via session.add() to chain the pathway."
        )


class EmptySessionError(PyPyrusMapError):
    """Raised when render() or build() is called on an empty session."""

    def __init__(self):
        super().__init__(
            "Cannot render an empty session. "
            "Call session.add('metabolite_id') at least once before rendering."
        )
