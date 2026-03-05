"""
tests/conftest.py
=================
Shared fixtures for the PyPyrusMap test suite.

Uses a minimal hand-built mock COBRApy model so that tests run without
needing iML1515 downloaded, and without network access.

The mock model represents a simplified carotenoid pathway:
    IPP + DMAPP → GGPP (ispA)
    GGPP → Phytoene (crtB)
    Phytoene → Lycopene (crtI)

Plus a glycolysis reaction for testing currency suppression:
    Glucose + ATP → G6P + ADP (HEX1)
"""

import pytest


# ---------------------------------------------------------------------------
# Mock COBRApy objects
# ---------------------------------------------------------------------------

class MockGene:
    def __init__(self, id_, name=""):
        self.id = id_
        self.name = name


class MockMetabolite:
    def __init__(self, id_, name="", compartment="c", formula=None):
        self.id = id_
        self.name = name
        self.compartment = compartment
        self.formula = formula


class MockReaction:
    def __init__(self, id_, name="", subsystem="", metabolites=None,
                 lower_bound=0.0, upper_bound=1000.0, genes=None):
        self.id = id_
        self.name = name
        self.subsystem = subsystem
        self.metabolites = metabolites or {}  # {MockMetabolite: float}
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.genes = genes or []


class MockSolution:
    def __init__(self, fluxes: dict):
        self.fluxes = fluxes
        self.status = "optimal"
        self.objective_value = 1.0


class MockMetaboliteList:
    def __init__(self, metabolites):
        self._mets = {m.id: m for m in metabolites}

    def get_by_id(self, id_):
        if id_ not in self._mets:
            raise KeyError(id_)
        return self._mets[id_]

    def __iter__(self):
        return iter(self._mets.values())


class MockModel:
    def __init__(self, id_, metabolites, reactions):
        self.id = id_
        self.metabolites = MockMetaboliteList(metabolites)
        self.reactions = reactions


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def carotenoid_mets():
    """Metabolites for the mock carotenoid pathway."""
    return {
        "ipp_c":    MockMetabolite("ipp_c",    "Isopentenyl diphosphate", "c"),
        "dmapp_c":  MockMetabolite("dmapp_c",  "Dimethylallyl diphosphate", "c"),
        "ggpp_c":   MockMetabolite("ggpp_c",   "Geranylgeranyl diphosphate", "c"),
        "phyto_c":  MockMetabolite("phyto_c",  "Phytoene", "c"),
        "lyco_c":   MockMetabolite("lyco_c",   "Lycopene", "c"),
        "ppi_c":    MockMetabolite("ppi_c",    "Pyrophosphate", "c"),
        "glc__D_c": MockMetabolite("glc__D_c", "D-Glucose", "c"),
        "g6p_c":    MockMetabolite("g6p_c",    "D-Glucose 6-phosphate", "c"),
        "atp_c":    MockMetabolite("atp_c",    "ATP", "c"),
        "adp_c":    MockMetabolite("adp_c",    "ADP", "c"),
    }


@pytest.fixture
def carotenoid_reactions(carotenoid_mets):
    """Reactions for the mock carotenoid pathway."""
    m = carotenoid_mets
    return [
        MockReaction(
            "ISPA",
            name="Geranylgeranyl diphosphate synthase",
            subsystem="Terpenoid biosynthesis",
            metabolites={
                m["ipp_c"]:   -3.0,
                m["dmapp_c"]: -1.0,
                m["ggpp_c"]:   1.0,
                m["ppi_c"]:    3.0,
            },
            lower_bound=0.0, upper_bound=1000.0,
            genes=[MockGene("ispA", "ispA")],
        ),
        MockReaction(
            "CRTB",
            name="Phytoene synthase",
            subsystem="Carotenoid biosynthesis",
            metabolites={
                m["ggpp_c"]:  -2.0,
                m["phyto_c"]:  1.0,
                m["ppi_c"]:    1.0,
            },
            lower_bound=0.0, upper_bound=1000.0,
            genes=[MockGene("crtB", "crtB")],
        ),
        MockReaction(
            "CRTI",
            name="Phytoene desaturase",
            subsystem="Carotenoid biosynthesis",
            metabolites={
                m["phyto_c"]: -1.0,
                m["lyco_c"]:   1.0,
            },
            lower_bound=0.0, upper_bound=1000.0,
            genes=[MockGene("crtI", "crtI")],
        ),
        MockReaction(
            "HEX1",
            name="Hexokinase",
            subsystem="Glycolysis",
            metabolites={
                m["glc__D_c"]: -1.0,
                m["atp_c"]:    -1.0,
                m["g6p_c"]:     1.0,
                m["adp_c"]:     1.0,
            },
            lower_bound=0.0, upper_bound=1000.0,
            genes=[MockGene("glk", "glk")],
        ),
    ]


@pytest.fixture
def mock_model(carotenoid_mets, carotenoid_reactions):
    """A minimal COBRApy-compatible mock model."""
    return MockModel(
        id_="mock_carotenoid_model",
        metabolites=list(carotenoid_mets.values()),
        reactions=carotenoid_reactions,
    )


@pytest.fixture
def mock_solution():
    """A mock FBA solution with realistic flux values."""
    return MockSolution(fluxes={
        "ISPA": 2.41,
        "CRTB": 1.20,
        "CRTI": 0.0,   # blocked — no desaturase activity
        "HEX1": 8.5,
    })
