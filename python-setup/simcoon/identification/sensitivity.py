"""
Sensitivity analysis tools for parameter identification.

This module provides functions for computing parameter sensitivities,
Jacobians, and correlation matrices to assess identifiability and
parameter interactions.
"""

from __future__ import annotations

from typing import Optional, List, Dict, Any

import numpy as np
import numpy.typing as npt

from .problem import IdentificationProblem


def compute_sensitivity(
    problem: IdentificationProblem,
    params: npt.NDArray[np.float64],
    eps: float = 1e-6,
    relative: bool = True,
) -> Dict[str, npt.NDArray[np.float64]]:
    """
    Compute sensitivity of outputs with respect to parameters.

    The sensitivity S_ij = d(output_i) / d(param_j) measures how much
    each output changes when a parameter is perturbed.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem
    params : np.ndarray
        Parameter values at which to compute sensitivity
    eps : float
        Perturbation size for finite differences
    relative : bool
        If True, compute relative sensitivity: (p/y) * dy/dp

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping output names to sensitivity matrices.
        Each matrix has shape (n_outputs, n_params).

    Examples
    --------
    >>> sensitivities = compute_sensitivity(problem, optimal_params)
    >>> print(f"Stress sensitivity to E: {sensitivities['stress'][:, 0]}")
    """
    full_params = problem._expand_params(params)
    n_params = problem.n_params

    # Get baseline predictions
    baseline = problem.simulate(full_params)

    sensitivities = {}

    for output_name, baseline_output in baseline.items():
        n_outputs = len(baseline_output)
        sens_matrix = np.zeros((n_outputs, n_params))

        for j in range(n_params):
            # Perturb parameter j
            params_plus = params.copy()
            params_plus[j] += eps
            full_plus = problem._expand_params(params_plus)

            # Get perturbed prediction
            perturbed = problem.simulate(full_plus)
            perturbed_output = perturbed.get(output_name, baseline_output)

            # Interpolate if needed
            if len(perturbed_output) != len(baseline_output):
                perturbed_output = np.interp(
                    np.linspace(0, 1, len(baseline_output)),
                    np.linspace(0, 1, len(perturbed_output)),
                    perturbed_output
                )

            # Compute sensitivity
            dy = perturbed_output - baseline_output
            dp = eps

            if relative:
                # Relative sensitivity: (p/y) * dy/dp
                scale = params[j] / (baseline_output + 1e-10)
                sens_matrix[:, j] = scale * dy / dp
            else:
                sens_matrix[:, j] = dy / dp

        sensitivities[output_name] = sens_matrix

    return sensitivities


def compute_jacobian(
    problem: IdentificationProblem,
    params: npt.NDArray[np.float64],
    eps: float = 1e-6,
) -> npt.NDArray[np.float64]:
    """
    Compute the Jacobian matrix of residuals with respect to parameters.

    J_ij = d(residual_i) / d(param_j)

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem
    params : np.ndarray
        Parameter values
    eps : float
        Perturbation size

    Returns
    -------
    np.ndarray
        Jacobian matrix of shape (n_residuals, n_params)
    """
    n_params = problem.n_params

    # Get baseline residuals
    residuals_base = problem.residual_vector(params)
    n_residuals = len(residuals_base)

    jacobian = np.zeros((n_residuals, n_params))

    for j in range(n_params):
        params_plus = params.copy()
        params_plus[j] += eps

        residuals_plus = problem.residual_vector(params_plus)

        jacobian[:, j] = (residuals_plus - residuals_base) / eps

    return jacobian


def correlation_matrix(
    problem: IdentificationProblem,
    params: npt.NDArray[np.float64],
    eps: float = 1e-6,
) -> npt.NDArray[np.float64]:
    """
    Compute parameter correlation matrix from the Jacobian.

    High correlations between parameters indicate potential
    identifiability issues (parameter coupling).

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem
    params : np.ndarray
        Parameter values
    eps : float
        Perturbation size for Jacobian computation

    Returns
    -------
    np.ndarray
        Correlation matrix of shape (n_params, n_params)

    Notes
    -----
    The correlation matrix is computed from the covariance matrix:
    corr_ij = cov_ij / sqrt(cov_ii * cov_jj)

    The covariance is estimated as: cov = inv(J^T * J)
    """
    jacobian = compute_jacobian(problem, params, eps)

    # J^T * J
    jtj = jacobian.T @ jacobian

    # Covariance (pseudo-inverse for numerical stability)
    try:
        cov = np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(jtj)

    # Correlation from covariance
    n_params = len(params)
    corr = np.zeros((n_params, n_params))

    for i in range(n_params):
        for j in range(n_params):
            denom = np.sqrt(cov[i, i] * cov[j, j])
            if denom > 1e-10:
                corr[i, j] = cov[i, j] / denom
            else:
                corr[i, j] = 0.0 if i != j else 1.0

    return corr


def parameter_uncertainty(
    problem: IdentificationProblem,
    params: npt.NDArray[np.float64],
    residual_variance: Optional[float] = None,
    eps: float = 1e-6,
) -> npt.NDArray[np.float64]:
    """
    Estimate parameter standard errors from the Jacobian.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem
    params : np.ndarray
        Optimal parameter values
    residual_variance : float, optional
        Variance of residuals. If None, estimated from data.
    eps : float
        Perturbation size for Jacobian

    Returns
    -------
    np.ndarray
        Standard errors for each parameter
    """
    jacobian = compute_jacobian(problem, params, eps)
    residuals = problem.residual_vector(params)

    n_data = len(residuals)
    n_params = len(params)

    # Estimate residual variance if not provided
    if residual_variance is None:
        dof = n_data - n_params
        if dof > 0:
            residual_variance = np.sum(residuals ** 2) / dof
        else:
            residual_variance = np.sum(residuals ** 2) / n_data

    # Covariance = sigma^2 * inv(J^T * J)
    jtj = jacobian.T @ jacobian

    try:
        cov = np.linalg.inv(jtj) * residual_variance
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(jtj) * residual_variance

    # Standard errors are sqrt of diagonal
    std_errors = np.sqrt(np.diag(cov))

    return std_errors


def identifiability_check(
    problem: IdentificationProblem,
    params: npt.NDArray[np.float64],
    eps: float = 1e-6,
    threshold: float = 0.9,
) -> Dict[str, Any]:
    """
    Check parameter identifiability.

    Analyzes the Jacobian to determine if all parameters can be
    uniquely identified from the available data.

    Parameters
    ----------
    problem : IdentificationProblem
        The identification problem
    params : np.ndarray
        Parameter values
    eps : float
        Perturbation size
    threshold : float
        Correlation threshold for flagging issues

    Returns
    -------
    dict
        Dictionary containing:
        - 'identifiable': bool, overall identifiability
        - 'condition_number': float, condition number of J^T*J
        - 'singular_values': np.ndarray, singular values of Jacobian
        - 'high_correlations': list of (i, j, corr) tuples
        - 'recommendations': list of strings
    """
    jacobian = compute_jacobian(problem, params, eps)
    corr = correlation_matrix(problem, params, eps)

    # Singular value decomposition
    U, S, Vt = np.linalg.svd(jacobian, full_matrices=False)

    # Condition number
    if S[-1] > 1e-10:
        condition = S[0] / S[-1]
    else:
        condition = np.inf

    # Find high correlations
    n_params = len(params)
    high_corr = []
    for i in range(n_params):
        for j in range(i + 1, n_params):
            if abs(corr[i, j]) > threshold:
                high_corr.append((
                    problem.parameter_names[i],
                    problem.parameter_names[j],
                    corr[i, j]
                ))

    # Recommendations
    recommendations = []
    identifiable = True

    if condition > 1e6:
        recommendations.append(
            f"High condition number ({condition:.2e}) indicates ill-conditioning. "
            "Consider removing or fixing some parameters."
        )
        identifiable = False

    for p1, p2, c in high_corr:
        recommendations.append(
            f"Parameters '{p1}' and '{p2}' are highly correlated ({c:.3f}). "
            "They may not be independently identifiable."
        )
        identifiable = False

    if S[-1] < 1e-10:
        recommendations.append(
            "Near-zero singular value detected. Some parameter combinations "
            "may not affect the output."
        )
        identifiable = False

    if not recommendations:
        recommendations.append("All parameters appear identifiable.")

    return {
        'identifiable': identifiable,
        'condition_number': condition,
        'singular_values': S,
        'high_correlations': high_corr,
        'correlation_matrix': corr,
        'recommendations': recommendations,
    }
