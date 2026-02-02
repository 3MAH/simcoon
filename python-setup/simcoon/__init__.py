"""
Simcoon - A library for the simulation of heterogeneous materials.

This package provides tools for:
- Constitutive modeling (UMATs)
- Homogenization and micromechanics
- Material parameter identification
- Continuum mechanics operations

Modules
-------
solver
    Python 0D material point solver with Block/Step/Solver classes
identification
    Material parameter identification and calibration tools
properties
    Elastic properties computation and analysis
odf
    Orientation Distribution Function tools
pdf
    Probability Distribution Function tools
parameter
    Parameter management for DOE and optimization
constant
    Test-dependent constants management
data
    Experimental/numerical data handling

Quick Start
-----------
>>> import numpy as np
>>> from simcoon.solver import Solver, Block, StepMeca
>>>
>>> # Define material properties (E, nu, alpha)
>>> props = np.array([210000.0, 0.3, 1e-5])
>>>
>>> # Create loading step
>>> step = StepMeca(
...     DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
...     control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'],
...     Dn_init=50
... )
>>>
>>> # Create block and solver
>>> block = Block(steps=[step], umat_name='ELISO', props=props, nstatev=1)
>>> solver = Solver(blocks=[block])
>>> history = solver.solve()
"""

from simcoon._core import *
from simcoon.__version__ import __version__

# Alias for backward compatibility
from simcoon import _core as simmit

# Import submodules
from simcoon import solver
from simcoon import identification
from simcoon import properties
from simcoon import odf
from simcoon import pdf
from simcoon import parameter
from simcoon import constant
from simcoon import data
