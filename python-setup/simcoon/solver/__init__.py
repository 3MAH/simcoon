"""In-memory material-point solver API.

Example
-------
>>> from simcoon import solver
>>> import numpy as np
>>> step = solver.StepMeca(control=['strain'] + ['stress']*5,
...                        value=np.array([0.01, 0, 0, 0, 0, 0]), ninc=100)
>>> res = solver.solve(step, "ELISO", [70000., 0.3, 1.E-5], 1)
>>> res["Stress"][0]     # sigma_11 history, fedoo-style (6, N) layout

simcoon reads no legacy text format (path.txt, material.dat, N<kind>.dat) any more:
convert such a directory once with ``scripts/legacy_to_json.py`` and load the JSON with
:func:`load_simulation_json`.
"""

from .maps import (
    BLOCK_TYPES,
    CONTROL_TYPES,
    CORATE_TYPES,
    STEP_MODES,
    TANGENT_MODES,
    THERMAL_CONTROL,
    tangent_algorithmic,
    tangent_closest_point,
    tangent_continuum,
    tangent_default,
    tangent_none,
)
from .blocks import Block, StepMeca, StepThermomeca
from .results import SolverResults
from .core import solve
from .io import (
    load_material_json,
    load_path_json,
    load_simulation_json,
    save_material_json,
    save_path_json,
)
# Phase/geometry dataclasses and their JSON I/O for micromechanics: a composite is
# described, saved and reloaded here, and handed to the extension by solve() / L_eff.
from . import micromechanics

__all__ = [
    "micromechanics",
    "Block", "StepMeca", "StepThermomeca", "SolverResults", "solve",
    "BLOCK_TYPES", "CONTROL_TYPES", "CORATE_TYPES", "STEP_MODES",
    "TANGENT_MODES", "THERMAL_CONTROL",
    "tangent_none", "tangent_continuum", "tangent_algorithmic",
    "tangent_closest_point", "tangent_default",
    "load_material_json", "load_path_json", "load_simulation_json",
    "save_material_json", "save_path_json",
]
