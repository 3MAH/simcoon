"""
Identification utilities for simcoon.

Provides ``identification()`` (wraps ``scipy.optimize.differential_evolution``)
and ``calc_cost()`` (multi-level weighted cost function for parameter
identification from multiple tests).
"""

from typing import Callable, List, Optional, Union

import numpy as np
from scipy.optimize import differential_evolution


def identification(
    cost_fn: Callable,
    parameters: list,
    method: str = "differential_evolution",
    **kwargs,
):
    """
    Run parameter identification using scipy global optimization.

    This function wraps ``scipy.optimize.differential_evolution`` using the
    bounds defined in the ``Parameter`` objects. After optimization, the
    identified values are written back to each ``Parameter.value``.

    Parameters
    ----------
    cost_fn : callable
        Cost function ``f(x) -> float`` where ``x`` is a 1-D array of
        parameter values (same order as *parameters*).
    parameters : list of Parameter
        Parameter objects whose ``.bounds`` define the search space.
    method : str, optional
        Optimization method. Currently only ``"differential_evolution"``
        is supported. Default: ``"differential_evolution"``.
    **kwargs
        Extra keyword arguments forwarded to the scipy optimizer.
        Common options: ``maxiter``, ``popsize``, ``tol``, ``seed``,
        ``polish`` (default ``True``), ``disp``.

    Returns
    -------
    scipy.optimize.OptimizeResult
        The full scipy result object.  ``result.x`` contains the
        identified parameter values.

    Examples
    --------
    >>> from simcoon.parameter import Parameter
    >>> from simcoon.identify import identification
    >>> params = [
    ...     Parameter(0, bounds=(10000, 200000), key="@Ef",
    ...               sim_input_files=["material.dat"]),
    ... ]
    >>> result = identification(my_cost, params, seed=42)
    >>> print(f"Identified E = {params[0].value:.0f}")
    """
    if method != "differential_evolution":
        raise ValueError(
            f"Unknown method '{method}'. "
            "Currently only 'differential_evolution' is supported."
        )

    bounds = [p.bounds for p in parameters]

    # Sensible defaults
    kwargs.setdefault("polish", True)

    result = differential_evolution(cost_fn, bounds=bounds, **kwargs)

    # Write back identified values
    for p, val in zip(parameters, result.x):
        p.value = val

    return result


def calc_cost(
    y_exp: List[np.ndarray],
    y_num: List[np.ndarray],
    w_test: Optional[np.ndarray] = None,
    w_response: Optional[List[np.ndarray]] = None,
    w_point: Optional[List[np.ndarray]] = None,
    metric: str = "mse",
):
    """
    Compute a weighted cost between experimental and numerical data.

    The data is organized per test (experiment), each test being a 2-D
    array of shape ``(n_points, n_responses)``.  Three levels of weights
    are combined multiplicatively, mirroring the classical identification
    cost function:

    .. math::

        C = \\sum_i \\sum_k \\sum_j
            W^{\\text{test}}_i \\;
            W^{\\text{resp}}_{i,k} \\;
            W^{\\text{pt}}_{i,k,j} \\;
            \\bigl(y^{\\exp}_{i,k,j} - y^{\\num}_{i,k,j}\\bigr)^2

    Parameters
    ----------
    y_exp : list of ndarray, each (n_points, n_responses) or (n_points,)
        Experimental data, one array per test.
    y_num : list of ndarray, same shapes as *y_exp*
        Numerical (model) data, one array per test.
    w_test : ndarray of shape (n_tests,), optional
        Weight per test.  ``None`` means uniform (1.0).
    w_response : list of ndarray, each (n_responses,), optional
        Weight per response type, for each test.  ``None`` means uniform.
    w_point : list of ndarray, each (n_points, n_responses), optional
        Weight per individual data point.  ``None`` means uniform.
    metric : str, optional
        Cost metric.  Built-in (numpy):

        - ``"mse"`` — Mean Squared Error (default)
        - ``"nmse"`` — Normalized MSE (divided by variance of all data)
        - ``"nmse_per_response"`` — NMSE computed per response column,
          then averaged.  Each column is divided by its own
          ``sum(y_exp²)``, balancing responses of different magnitudes
          (e.g., force in N vs displacement in mm).
        - ``"rmse"`` — Root Mean Squared Error
        - ``"mae"`` — Mean Absolute Error

        With scikit-learn installed (``pip install simcoon[identify]``):
        ``"r2"`` and any ``sklearn.metrics`` function that accepts
        ``sample_weight``.

    Returns
    -------
    float
        The scalar cost value.

    Examples
    --------
    Simple MSE over two tensile tests:

    >>> import numpy as np
    >>> from simcoon.identify import calc_cost
    >>> y_exp = [np.array([100, 200, 300]), np.array([150, 250])]
    >>> y_num = [np.array([105, 195, 310]), np.array([148, 260])]
    >>> calc_cost(y_exp, y_num)  # simple MSE
    50.6

    NMSE per response (balances force vs displacement):

    >>> y_exp = [np.column_stack([force_exp, disp_exp])]
    >>> y_num = [np.column_stack([force_num, disp_num])]
    >>> calc_cost(y_exp, y_num, metric='nmse_per_response')
    """
    n_tests = len(y_exp)
    if len(y_num) != n_tests:
        raise ValueError(
            f"y_exp has {n_tests} tests but y_num has {len(y_num)}"
        )

    # Ensure 2-D
    y_exp_2d = [np.atleast_2d(y) if y.ndim == 1 else y for y in y_exp]
    y_num_2d = [np.atleast_2d(y) if y.ndim == 1 else y for y in y_num]
    # atleast_2d on 1-D gives (1, N) — we want (N, 1)
    y_exp_2d = [y.T if y_orig.ndim == 1 else y
                for y, y_orig in zip(y_exp_2d, y_exp)]
    y_num_2d = [y.T if y_orig.ndim == 1 else y
                for y, y_orig in zip(y_num_2d, y_num)]

    # Handle nmse_per_response separately (per-column computation)
    if metric == "nmse_per_response":
        return _nmse_per_response(
            y_exp_2d, y_num_2d, n_tests, w_test, w_response, w_point,
        )

    # Build flat vectors and weight vector
    all_exp = []
    all_num = []
    all_w = []

    for i in range(n_tests):
        exp_i = y_exp_2d[i]
        num_i = y_num_2d[i]
        n_pts, n_resp = exp_i.shape

        if num_i.shape != exp_i.shape:
            raise ValueError(
                f"Test {i}: y_exp shape {exp_i.shape} != "
                f"y_num shape {num_i.shape}"
            )

        W_i = np.ones_like(exp_i)

        # Level 1: per-test weight
        if w_test is not None:
            W_i *= w_test[i]

        # Level 2: per-response weight
        if w_response is not None:
            for k in range(n_resp):
                W_i[:, k] *= w_response[i][k]

        # Level 3: per-point weight
        if w_point is not None:
            W_i *= np.abs(w_point[i])

        all_exp.append(exp_i.ravel())
        all_num.append(num_i.ravel())
        all_w.append(W_i.ravel())

    y_exp_flat = np.concatenate(all_exp)
    y_num_flat = np.concatenate(all_num)
    w_flat = np.concatenate(all_w)

    return _compute_metric(y_exp_flat, y_num_flat, w_flat, metric)


def _nmse_per_response(
    y_exp_2d, y_num_2d, n_tests, w_test, w_response, w_point,
):
    """Compute NMSE independently per response column, then average."""
    # Collect all columns across tests
    # First pass: determine n_resp per test and validate
    col_costs = {}  # col_index -> list of (mse_col, denom_col, weight)

    for i in range(n_tests):
        exp_i = y_exp_2d[i]
        num_i = y_num_2d[i]
        n_pts, n_resp = exp_i.shape

        if num_i.shape != exp_i.shape:
            raise ValueError(
                f"Test {i}: y_exp shape {exp_i.shape} != "
                f"y_num shape {num_i.shape}"
            )

        wt = w_test[i] if w_test is not None else 1.0

        for k in range(n_resp):
            wr = w_response[i][k] if w_response is not None else 1.0

            residuals = exp_i[:, k] - num_i[:, k]
            wp = (np.abs(w_point[i][:, k]) if w_point is not None
                  else np.ones(n_pts))

            w_combined = wt * wr * wp
            weighted_sse = np.sum(w_combined * residuals ** 2)
            denom = np.sum(exp_i[:, k] ** 2)

            if k not in col_costs:
                col_costs[k] = {"sse": 0.0, "denom": 0.0}
            col_costs[k]["sse"] += weighted_sse
            col_costs[k]["denom"] += denom

    # Average NMSE across response columns
    nmse_values = []
    for k in sorted(col_costs):
        denom = col_costs[k]["denom"]
        if denom > 1e-30:
            nmse_values.append(col_costs[k]["sse"] / denom)
        else:
            nmse_values.append(col_costs[k]["sse"])

    return float(np.mean(nmse_values))


def _compute_metric(
    y_exp: np.ndarray,
    y_num: np.ndarray,
    w: np.ndarray,
    metric: str,
) -> float:
    """Compute the requested metric, using sklearn if available."""
    residuals = y_exp - y_num

    # Built-in metrics (numpy only)
    if metric == "mse":
        return float(np.average(residuals ** 2, weights=w))

    if metric == "nmse":
        mse = np.average(residuals ** 2, weights=w)
        var = np.average((y_exp - np.average(y_exp, weights=w)) ** 2, weights=w)
        return float(mse / var) if var > 1e-30 else float(mse)

    if metric == "rmse":
        return float(np.sqrt(np.average(residuals ** 2, weights=w)))

    if metric == "mae":
        return float(np.average(np.abs(residuals), weights=w))

    # sklearn metrics
    try:
        import sklearn.metrics as skm
    except ImportError:
        raise ImportError(
            f"Metric '{metric}' requires scikit-learn. Install it with:\n"
            "  pip install scikit-learn\n"
            "  # or\n"
            "  conda install -c conda-forge scikit-learn"
        )

    if metric == "r2":
        return float(skm.r2_score(y_exp, y_num, sample_weight=w))

    if metric == "mean_squared_error":
        return float(skm.mean_squared_error(y_exp, y_num, sample_weight=w))

    if metric == "mean_absolute_error":
        return float(skm.mean_absolute_error(y_exp, y_num, sample_weight=w))

    if hasattr(skm, metric):
        fn = getattr(skm, metric)
        return float(fn(y_exp, y_num, sample_weight=w))

    raise ValueError(
        f"Unknown metric '{metric}'. Built-in: mse, nmse, rmse, mae. "
        f"With scikit-learn: r2, mean_squared_error, mean_absolute_error, etc."
    )
