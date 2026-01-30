"""
Identification problem definition classes.

This module provides the core data structures for defining material parameter
identification problems.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import (
    Callable, Dict, List, Optional, Tuple, Union, Any
)

import numpy as np
import numpy.typing as npt


@dataclass
class ParameterSpec:
    """
    Specification for a parameter to be identified.

    Attributes
    ----------
    name : str
        Name/identifier of the parameter
    bounds : Tuple[float, float]
        Lower and upper bounds for the parameter
    initial : float, optional
        Initial guess for the parameter. If None, uses midpoint of bounds.
    scale : float, optional
        Scaling factor for normalization. If None, uses (max - min).
    fixed : bool
        If True, parameter is held fixed at initial value during optimization.

    Examples
    --------
    >>> param = ParameterSpec(name='E', bounds=(100000, 300000), initial=200000)
    >>> param.normalized_initial
    0.5
    """

    name: str
    bounds: Tuple[float, float]
    initial: Optional[float] = None
    scale: Optional[float] = None
    fixed: bool = False

    def __post_init__(self):
        if self.initial is None:
            self.initial = 0.5 * (self.bounds[0] + self.bounds[1])
        if self.scale is None:
            self.scale = self.bounds[1] - self.bounds[0]

    @property
    def normalized_initial(self) -> float:
        """Return initial value normalized to [0, 1]."""
        return (self.initial - self.bounds[0]) / (self.bounds[1] - self.bounds[0])

    def denormalize(self, normalized_value: float) -> float:
        """Convert normalized [0, 1] value to actual parameter value."""
        return self.bounds[0] + normalized_value * (self.bounds[1] - self.bounds[0])

    def normalize(self, value: float) -> float:
        """Convert actual parameter value to normalized [0, 1]."""
        return (value - self.bounds[0]) / (self.bounds[1] - self.bounds[0])


@dataclass
class OptimizationResult:
    """
    Container for optimization results.

    Attributes
    ----------
    x : np.ndarray
        Optimal parameter values
    cost : float
        Final cost function value
    success : bool
        Whether optimization converged successfully
    message : str
        Optimization status message
    n_iterations : int
        Number of iterations performed
    n_function_evals : int
        Number of function evaluations
    parameter_names : List[str]
        Names of the optimized parameters
    history : List[Dict], optional
        History of parameter values and costs during optimization
    jacobian : np.ndarray, optional
        Jacobian matrix at the solution
    covariance : np.ndarray, optional
        Parameter covariance matrix (if available)
    """

    x: npt.NDArray[np.float64]
    cost: float
    success: bool
    message: str
    n_iterations: int
    n_function_evals: int
    parameter_names: List[str] = field(default_factory=list)
    history: List[Dict[str, Any]] = field(default_factory=list)
    jacobian: Optional[npt.NDArray[np.float64]] = None
    covariance: Optional[npt.NDArray[np.float64]] = None

    def __repr__(self) -> str:
        params_str = ", ".join(
            f"{name}={val:.6g}"
            for name, val in zip(self.parameter_names, self.x)
        )
        return (
            f"OptimizationResult(\n"
            f"  success={self.success},\n"
            f"  cost={self.cost:.6e},\n"
            f"  parameters={{ {params_str} }},\n"
            f"  iterations={self.n_iterations},\n"
            f"  message='{self.message}'\n"
            f")"
        )

    @classmethod
    def from_scipy(cls, scipy_result, parameter_names: List[str] = None):
        """Create OptimizationResult from scipy.optimize result object."""
        # Handle different scipy result types
        if hasattr(scipy_result, 'nit'):
            n_iterations = scipy_result.nit
        else:
            n_iterations = 0

        if hasattr(scipy_result, 'nfev'):
            n_function_evals = scipy_result.nfev
        elif hasattr(scipy_result, 'njev'):
            n_function_evals = scipy_result.njev
        else:
            n_function_evals = 0

        # Get cost
        if hasattr(scipy_result, 'cost'):
            cost = scipy_result.cost
        elif hasattr(scipy_result, 'fun'):
            cost = float(scipy_result.fun)
        else:
            cost = 0.0

        # Get jacobian
        jacobian = getattr(scipy_result, 'jac', None)

        return cls(
            x=np.asarray(scipy_result.x),
            cost=cost,
            success=scipy_result.success,
            message=str(scipy_result.message),
            n_iterations=n_iterations,
            n_function_evals=n_function_evals,
            parameter_names=parameter_names or [],
            jacobian=jacobian,
        )


class IdentificationProblem:
    """
    Define a material parameter identification/calibration problem.

    This class encapsulates all information needed to perform parameter
    identification: parameters to optimize, experimental data to match,
    and a simulation function that maps parameters to predictions.

    Parameters
    ----------
    parameters : List[Union[ParameterSpec, Dict]]
        List of parameters to identify. Can be ParameterSpec objects or
        dictionaries with keys: 'name', 'bounds', 'initial' (optional).
    simulate : Callable[[np.ndarray], Dict[str, np.ndarray]]
        Function that takes parameter array and returns dict of predictions.
        Keys should match keys in exp_data.
    exp_data : Dict[str, np.ndarray]
        Experimental data to match. Keys are data names, values are arrays.
    weights : Dict[str, float], optional
        Weights for each data type in cost function.
    cost_type : str, optional
        Cost function type: 'mse' (default), 'mae', 'r2', 'huber'

    Attributes
    ----------
    n_params : int
        Number of parameters being optimized

    Examples
    --------
    >>> def my_simulation(params):
    ...     E, nu = params
    ...     # Run simulation and return predictions
    ...     return {'stress': computed_stress, 'strain': computed_strain}
    >>>
    >>> problem = IdentificationProblem(
    ...     parameters=[
    ...         {'name': 'E', 'bounds': (100000, 300000)},
    ...         {'name': 'nu', 'bounds': (0.2, 0.4)},
    ...     ],
    ...     simulate=my_simulation,
    ...     exp_data={'stress': exp_stress},
    ... )
    """

    def __init__(
        self,
        parameters: List[Union[ParameterSpec, Dict]],
        simulate: Callable[[npt.NDArray[np.float64]], Dict[str, npt.NDArray[np.float64]]],
        exp_data: Dict[str, npt.NDArray[np.float64]],
        weights: Optional[Dict[str, float]] = None,
        cost_type: str = 'mse',
    ):
        # Convert dict parameters to ParameterSpec
        self.parameters: List[ParameterSpec] = []
        for p in parameters:
            if isinstance(p, ParameterSpec):
                self.parameters.append(p)
            elif isinstance(p, dict):
                self.parameters.append(ParameterSpec(**p))
            else:
                raise TypeError(f"Parameter must be ParameterSpec or dict, got {type(p)}")

        self.simulate = simulate
        self.exp_data = exp_data
        self.weights = weights or {k: 1.0 for k in exp_data.keys()}
        self.cost_type = cost_type

        # Precompute useful quantities
        self._active_params = [p for p in self.parameters if not p.fixed]
        self._n_params = len(self._active_params)

        # Function evaluation counter
        self._n_evals = 0

    @property
    def n_params(self) -> int:
        """Number of active (non-fixed) parameters."""
        return self._n_params

    @property
    def parameter_names(self) -> List[str]:
        """Names of active parameters."""
        return [p.name for p in self._active_params]

    def get_bounds(self) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """
        Get parameter bounds as arrays.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            Lower and upper bounds arrays
        """
        lower = np.array([p.bounds[0] for p in self._active_params])
        upper = np.array([p.bounds[1] for p in self._active_params])
        return lower, upper

    def get_bounds_list(self) -> List[Tuple[float, float]]:
        """
        Get parameter bounds as list of tuples.

        Returns
        -------
        List[Tuple[float, float]]
            List of (lower, upper) bound tuples
        """
        return [p.bounds for p in self._active_params]

    def get_initial(self) -> npt.NDArray[np.float64]:
        """
        Get initial parameter values.

        Returns
        -------
        np.ndarray
            Initial parameter values
        """
        return np.array([p.initial for p in self._active_params])

    def _expand_params(self, active_params: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Expand active parameters to full parameter vector (including fixed)."""
        full_params = []
        active_idx = 0
        for p in self.parameters:
            if p.fixed:
                full_params.append(p.initial)
            else:
                full_params.append(active_params[active_idx])
                active_idx += 1
        return np.array(full_params)

    def cost_function(self, params: npt.NDArray[np.float64]) -> float:
        """
        Evaluate the cost function for given parameters.

        Parameters
        ----------
        params : np.ndarray
            Parameter values (active parameters only)

        Returns
        -------
        float
            Cost function value
        """
        self._n_evals += 1

        # Expand to full parameter vector
        full_params = self._expand_params(params)

        # Run simulation
        try:
            predictions = self.simulate(full_params)
        except Exception as e:
            # Return large cost on simulation failure
            return 1e10

        # Compute cost
        total_cost = 0.0
        for key, exp in self.exp_data.items():
            if key not in predictions:
                continue

            pred = predictions[key]
            weight = self.weights.get(key, 1.0)

            # Interpolate if lengths don't match
            if len(pred) != len(exp):
                pred = np.interp(
                    np.linspace(0, 1, len(exp)),
                    np.linspace(0, 1, len(pred)),
                    pred
                )

            # Compute cost based on type
            if self.cost_type == 'mse':
                cost = np.mean((pred - exp) ** 2)
            elif self.cost_type == 'mae':
                cost = np.mean(np.abs(pred - exp))
            elif self.cost_type == 'r2':
                ss_res = np.sum((exp - pred) ** 2)
                ss_tot = np.sum((exp - np.mean(exp)) ** 2)
                cost = 1.0 - (ss_res / (ss_tot + 1e-10))
                cost = 1.0 - cost  # Minimize 1 - R^2
            else:
                cost = np.mean((pred - exp) ** 2)

            total_cost += weight * cost

        return total_cost

    def residual_vector(self, params: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Compute residual vector for least-squares optimization.

        Parameters
        ----------
        params : np.ndarray
            Parameter values (active parameters only)

        Returns
        -------
        np.ndarray
            Residual vector (pred - exp) for all data points
        """
        self._n_evals += 1

        # Expand to full parameter vector
        full_params = self._expand_params(params)

        # Run simulation
        try:
            predictions = self.simulate(full_params)
        except Exception:
            # Return large residuals on failure
            total_size = sum(len(v) for v in self.exp_data.values())
            return np.ones(total_size) * 1e5

        # Build residual vector
        residuals = []
        for key, exp in self.exp_data.items():
            if key not in predictions:
                residuals.extend([1e5] * len(exp))
                continue

            pred = predictions[key]
            weight = np.sqrt(self.weights.get(key, 1.0))

            # Interpolate if lengths don't match
            if len(pred) != len(exp):
                pred = np.interp(
                    np.linspace(0, 1, len(exp)),
                    np.linspace(0, 1, len(pred)),
                    pred
                )

            residuals.extend(weight * (pred - exp))

        return np.array(residuals)

    def reset_counter(self):
        """Reset function evaluation counter."""
        self._n_evals = 0

    @property
    def n_function_evals(self) -> int:
        """Number of function evaluations."""
        return self._n_evals
