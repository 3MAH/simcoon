"""
Optimization algorithms for parameter identification.

This module provides wrappers around scipy.optimize functions tailored
for material parameter identification problems.
"""

from __future__ import annotations

from typing import Optional, Dict, Any

import numpy as np
import numpy.typing as npt

from .problem import IdentificationProblem, OptimizationResult


def levenberg_marquardt(
    problem: IdentificationProblem,
    x0: Optional[npt.NDArray[np.float64]] = None,
    ftol: float = 1e-8,
    xtol: float = 1e-8,
    gtol: float = 1e-8,
    max_nfev: Optional[int] = None,
    verbose: int = 0,
    **kwargs,
) -> OptimizationResult:
    """
    Levenberg-Marquardt optimization for least-squares problems.

    This is well-suited for parameter identification where we have
    experimental data and want to minimize the sum of squared residuals.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem to solve
    x0 : np.ndarray, optional
        Initial parameter guess. If None, uses problem's initial values.
    ftol : float
        Tolerance for termination by change of cost function
    xtol : float
        Tolerance for termination by change in parameters
    gtol : float
        Tolerance for termination by norm of gradient
    max_nfev : int, optional
        Maximum number of function evaluations
    verbose : int
        Verbosity level (0, 1, or 2)
    **kwargs
        Additional arguments passed to scipy.optimize.least_squares

    Returns
    -------
    OptimizationResult
        Optimization result containing optimal parameters and diagnostics

    Examples
    --------
    >>> result = levenberg_marquardt(problem)
    >>> print(f"Optimal E: {result.x[0]:.0f}")
    """
    from scipy.optimize import least_squares

    if x0 is None:
        x0 = problem.get_initial()

    bounds = problem.get_bounds()
    problem.reset_counter()

    result = least_squares(
        problem.residual_vector,
        x0,
        bounds=bounds,
        method='trf',  # Trust Region Reflective, handles bounds
        ftol=ftol,
        xtol=xtol,
        gtol=gtol,
        max_nfev=max_nfev,
        verbose=verbose,
        **kwargs,
    )

    opt_result = OptimizationResult.from_scipy(result, problem.parameter_names)

    # Compute covariance if possible
    if result.success and result.jac is not None:
        try:
            # Covariance = inv(J^T * J) * sigma^2
            # sigma^2 estimated from residuals
            jac = result.jac
            residuals = result.fun
            n_data = len(residuals)
            n_params = len(result.x)
            dof = n_data - n_params

            if dof > 0:
                s2 = np.sum(residuals ** 2) / dof
                jtj = jac.T @ jac
                cov = np.linalg.inv(jtj) * s2
                opt_result.covariance = cov
        except Exception:
            pass  # Covariance estimation failed

    return opt_result


def differential_evolution(
    problem: IdentificationProblem,
    strategy: str = 'best1bin',
    maxiter: int = 1000,
    popsize: int = 15,
    tol: float = 1e-7,
    mutation: tuple = (0.5, 1.0),
    recombination: float = 0.7,
    seed: Optional[int] = None,
    workers: int = 1,
    polish: bool = True,
    verbose: bool = False,
    **kwargs,
) -> OptimizationResult:
    """
    Global optimization using differential evolution.

    This is well-suited for problems with many local minima or when
    good initial guesses are not available.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem to solve
    strategy : str
        Differential evolution strategy
    maxiter : int
        Maximum number of generations
    popsize : int
        Population size multiplier
    tol : float
        Convergence tolerance
    mutation : tuple
        Mutation constant range
    recombination : float
        Recombination constant
    seed : int, optional
        Random seed for reproducibility
    workers : int
        Number of parallel workers (-1 for all CPUs)
    polish : bool
        Whether to polish the result with L-BFGS-B
    verbose : bool
        Whether to print progress
    **kwargs
        Additional arguments passed to scipy.optimize.differential_evolution

    Returns
    -------
    OptimizationResult
        Optimization result

    Examples
    --------
    >>> result = differential_evolution(problem, maxiter=500)
    >>> print(f"Global optimum cost: {result.cost:.6f}")
    """
    from scipy.optimize import differential_evolution as scipy_de

    bounds = problem.get_bounds_list()
    problem.reset_counter()

    # Callback for verbose output
    callback = None
    if verbose:
        def callback(xk, convergence):
            cost = problem.cost_function(xk)
            print(f"Generation: cost = {cost:.6e}, convergence = {convergence:.6e}")

    result = scipy_de(
        problem.cost_function,
        bounds,
        strategy=strategy,
        maxiter=maxiter,
        popsize=popsize,
        tol=tol,
        mutation=mutation,
        recombination=recombination,
        seed=seed,
        workers=workers,
        polish=polish,
        callback=callback,
        **kwargs,
    )

    return OptimizationResult.from_scipy(result, problem.parameter_names)


def nelder_mead(
    problem: IdentificationProblem,
    x0: Optional[npt.NDArray[np.float64]] = None,
    maxiter: Optional[int] = None,
    maxfev: Optional[int] = None,
    xatol: float = 1e-4,
    fatol: float = 1e-4,
    adaptive: bool = True,
    verbose: bool = False,
    **kwargs,
) -> OptimizationResult:
    """
    Nelder-Mead simplex optimization.

    Derivative-free method that works well for smooth problems
    with a small number of parameters.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem to solve
    x0 : np.ndarray, optional
        Initial parameter guess
    maxiter : int, optional
        Maximum iterations
    maxfev : int, optional
        Maximum function evaluations
    xatol : float
        Absolute tolerance for parameters
    fatol : float
        Absolute tolerance for function value
    adaptive : bool
        Whether to use adaptive Nelder-Mead
    verbose : bool
        Whether to print progress
    **kwargs
        Additional arguments passed to scipy.optimize.minimize

    Returns
    -------
    OptimizationResult
        Optimization result
    """
    from scipy.optimize import minimize

    if x0 is None:
        x0 = problem.get_initial()

    bounds = problem.get_bounds_list()
    problem.reset_counter()

    options = {
        'maxiter': maxiter,
        'maxfev': maxfev,
        'xatol': xatol,
        'fatol': fatol,
        'adaptive': adaptive,
        'disp': verbose,
    }

    result = minimize(
        problem.cost_function,
        x0,
        method='Nelder-Mead',
        bounds=bounds,
        options=options,
        **kwargs,
    )

    return OptimizationResult.from_scipy(result, problem.parameter_names)


def hybrid_optimization(
    problem: IdentificationProblem,
    global_maxiter: int = 100,
    local_ftol: float = 1e-8,
    n_restarts: int = 3,
    verbose: bool = False,
    seed: Optional[int] = None,
) -> OptimizationResult:
    """
    Hybrid global + local optimization.

    Performs global exploration with differential evolution followed
    by local refinement with Levenberg-Marquardt. Multiple restarts
    help avoid local minima.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem to solve
    global_maxiter : int
        Maximum iterations for global search
    local_ftol : float
        Tolerance for local refinement
    n_restarts : int
        Number of random restarts
    verbose : bool
        Whether to print progress
    seed : int, optional
        Random seed

    Returns
    -------
    OptimizationResult
        Best result across all restarts

    Examples
    --------
    >>> result = hybrid_optimization(problem, n_restarts=5)
    >>> print(f"Best parameters: {result.x}")
    """
    if seed is not None:
        np.random.seed(seed)

    best_result = None
    best_cost = np.inf

    for i in range(n_restarts):
        if verbose:
            print(f"\n--- Restart {i + 1}/{n_restarts} ---")

        # Global search
        if verbose:
            print("Global search (differential evolution)...")

        global_result = differential_evolution(
            problem,
            maxiter=global_maxiter,
            polish=False,
            seed=seed + i if seed else None,
            verbose=verbose,
        )

        if verbose:
            print(f"Global result: cost = {global_result.cost:.6e}")

        # Local refinement
        if verbose:
            print("Local refinement (Levenberg-Marquardt)...")

        local_result = levenberg_marquardt(
            problem,
            x0=global_result.x,
            ftol=local_ftol,
            verbose=2 if verbose else 0,
        )

        if verbose:
            print(f"Local result: cost = {local_result.cost:.6e}")

        # Keep best
        if local_result.cost < best_cost:
            best_cost = local_result.cost
            best_result = local_result

    if verbose:
        print(f"\nBest overall cost: {best_cost:.6e}")

    return best_result
