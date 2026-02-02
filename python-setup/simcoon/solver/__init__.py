"""
Simcoon Python Solver Module.

This module provides a Python-based 0D material point solver that mirrors
the C++ solver architecture while enabling flexible material simulations
with Python control flow.

Classes
-------
StateVariables
    Base state variables class
StateVariablesM
    Mechanical state variables
StateVariablesT
    Thermomechanical state variables
Step
    Base loading step class
StepMeca
    Mechanical loading step
StepThermomeca
    Thermomechanical loading step
Block
    Loading block containing steps
Solver
    Main solver class with Newton-Raphson iteration

Micromechanics Classes (standalone, no _core required)
------------------------------------------------------
MaterialOrientation
    Material orientation via Euler angles
GeometryOrientation
    Geometry/phase orientation via Euler angles
Phase
    Generic phase for micromechanics homogenization
Layer
    Layer phase for laminate homogenization
Ellipsoid
    Ellipsoidal inclusion for Eshelby-based homogenization
Cylinder
    Cylindrical inclusion for micromechanics
Section
    Section/yarn for textile composite homogenization

Functions
---------
Lt_2_K
    Build 6x6 Jacobian for mechanical solver
Lth_2_K
    Build 7x7 Jacobian for thermomechanical solver

Examples
--------
>>> import numpy as np
>>> from simcoon.solver import Solver, Block, StepMeca, StateVariablesM
>>>
>>> # Material properties for ELISO (E, nu)
>>> props = np.array([210000.0, 0.3])
>>>
>>> # Uniaxial tension step
>>> step = StepMeca(
...     DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
...     control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
...     Dn_init=10
... )
>>>
>>> block = Block(
...     steps=[step],
...     umat_name="ELISO",
...     props=props,
...     nstatev=1
... )
>>>
>>> solver = Solver(blocks=[block])
>>> history = solver.solve()
"""

# Micromechanics classes and I/O functions (standalone - no _core required)
# These can be imported even without building simcoon._core
from .micromechanics import (
    # Data classes
    MaterialOrientation,
    GeometryOrientation,
    Phase,
    Layer,
    Ellipsoid,
    Cylinder,
    Section,
    # JSON I/O
    load_phases_json,
    save_phases_json,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
    load_cylinders_json,
    save_cylinders_json,
    load_sections_json,
    save_sections_json,
)

# Solver classes (require simcoon._core)
from .solver import (
    # Control type mappings
    CONTROL_TYPES,
    CORATE_TYPES,
    # History storage
    HistoryPoint,
    # State variable classes
    StateVariables,
    StateVariablesM,
    StateVariablesT,
    # Step classes
    Step,
    StepMeca,
    StepThermomeca,
    # Block class
    Block,
    # Solver class
    Solver,
    # Helper functions
    Lt_2_K,
    Lth_2_K,
)

# I/O functions (material/path JSON requires _core via lazy import)
from .io import (
    # JSON loading - Material/Path
    load_material_json,
    save_material_json,
    load_path_json,
    save_path_json,
    load_simulation_json,
)

__all__ = [
    # Control type mappings
    'CONTROL_TYPES',
    'CORATE_TYPES',
    # History storage
    'HistoryPoint',
    # State variable classes
    'StateVariables',
    'StateVariablesM',
    'StateVariablesT',
    # Step classes
    'Step',
    'StepMeca',
    'StepThermomeca',
    # Block class
    'Block',
    # Solver class
    'Solver',
    # Helper functions
    'Lt_2_K',
    'Lth_2_K',
    # JSON I/O - Material/Path
    'load_material_json',
    'save_material_json',
    'load_path_json',
    'save_path_json',
    'load_simulation_json',
    # Micromechanics data classes
    'MaterialOrientation',
    'GeometryOrientation',
    'Phase',
    'Layer',
    'Ellipsoid',
    'Cylinder',
    'Section',
    # Micromechanics JSON I/O
    'load_phases_json',
    'save_phases_json',
    'load_layers_json',
    'save_layers_json',
    'load_ellipsoids_json',
    'save_ellipsoids_json',
    'load_cylinders_json',
    'save_cylinders_json',
    'load_sections_json',
    'save_sections_json',
]
