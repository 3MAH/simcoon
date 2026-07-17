"""In-memory material-point solver API.

Example
-------
>>> from simcoon import solver
>>> import numpy as np
>>> step = solver.StepMeca(control=['strain'] + ['stress']*5,
...                        value=np.array([0.01, 0, 0, 0, 0, 0]), ninc=100)
>>> res = solver.solve(step, "ELISO", [70000., 0.3, 1.E-5], 1)
>>> res["Stress"][0]     # sigma_11 history, fedoo-style (6, N) layout

Legacy text files (path.txt / material.dat) are parsed into loading objects
with :func:`from_file` and :func:`material_from_file`:

>>> blocks, T_init = solver.from_file("data", "path.txt")
>>> res = solver.solve(blocks, T_init=T_init,
...                    **solver.material_from_file("data", "material.dat"))

(The pre-2.0 file-driven runner remains available as the low-level binding
``simcoon._core.solver`` -- it reads/writes files and returns nothing.)
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
from .files import from_file, material_from_file
from .io import (
    load_material_json,
    load_path_json,
    load_simulation_json,
    save_material_json,
    save_path_json,
)

__all__ = [
    "Block", "StepMeca", "StepThermomeca", "SolverResults", "solve",
    "from_file", "material_from_file",
    "BLOCK_TYPES", "CONTROL_TYPES", "CORATE_TYPES", "STEP_MODES",
    "TANGENT_MODES", "THERMAL_CONTROL",
    "tangent_none", "tangent_continuum", "tangent_algorithmic",
    "tangent_closest_point", "tangent_default",
    "load_material_json", "load_path_json", "load_simulation_json",
    "save_material_json", "save_path_json",
]
