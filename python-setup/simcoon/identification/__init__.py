"""
Simcoon Identification Module.

Material parameter identification and calibration tools using scipy.optimize.

Classes
-------
IdentificationProblem
    Main class for defining identification/calibration problems
OptimizationResult
    Result container for optimization runs

Functions
---------
levenberg_marquardt
    Levenberg-Marquardt optimization wrapper
differential_evolution
    Global optimization via differential evolution
hybrid_optimization
    Combined global + local optimization

Cost Functions
--------------
mse, mae, r2, weighted_mse
    Standard cost function metrics

Examples
--------
>>> import numpy as np
>>> from simcoon.identification import IdentificationProblem, levenberg_marquardt
>>> from simcoon.solver import Solver, Block, StepMeca
>>>
>>> # Define parameters to identify
>>> params = [
>>>     {'name': 'E', 'bounds': (150000, 250000), 'initial': 200000},
>>>     {'name': 'nu', 'bounds': (0.2, 0.4), 'initial': 0.3},
>>> ]
>>>
>>> # Define simulation function
>>> def simulate(param_values):
>>>     E, nu = param_values
>>>     props = np.array([E, nu, 0.0])
>>>     step = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]))
>>>     block = Block(steps=[step], umat_name='ELISO', props=props, nstatev=1)
>>>     solver = Solver(blocks=[block])
>>>     history = solver.solve()
>>>     return {'stress': np.array([h.sigma[0] for h in history])}
>>>
>>> # Create problem
>>> problem = IdentificationProblem(
>>>     parameters=params,
>>>     simulate=simulate,
>>>     exp_data={'stress': exp_stress_data},
>>> )
>>>
>>> # Run identification
>>> result = levenberg_marquardt(problem)
>>> print(f"Identified parameters: {result.x}")
"""

from .problem import (
    IdentificationProblem,
    OptimizationResult,
    ParameterSpec,
)

from .optimizers import (
    levenberg_marquardt,
    differential_evolution,
    hybrid_optimization,
    nelder_mead,
)

from .cost_functions import (
    mse,
    mae,
    r2,
    weighted_mse,
    huber_loss,
    CostFunction,
)

from .sensitivity import (
    compute_sensitivity,
    compute_jacobian,
    correlation_matrix,
)

__all__ = [
    # Problem definition
    'IdentificationProblem',
    'OptimizationResult',
    'ParameterSpec',
    # Optimizers
    'levenberg_marquardt',
    'differential_evolution',
    'hybrid_optimization',
    'nelder_mead',
    # Cost functions
    'mse',
    'mae',
    'r2',
    'weighted_mse',
    'huber_loss',
    'CostFunction',
    # Sensitivity analysis
    'compute_sensitivity',
    'compute_jacobian',
    'correlation_matrix',
]
