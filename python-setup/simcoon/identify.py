"""
Identification utilities for simcoon.

Provides ``identification()`` (bounded global search with
``scipy.optimize.differential_evolution``, or bounded multi-start least
squares with ``scipy.optimize.least_squares``), ``calc_cost()`` (multi-level
weighted cost function for parameter identification from multiple tests) and
``calc_residuals()`` (the residual vector whose sum of squares is that cost,
for the least-squares route).
"""

import os
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List, Optional

import numpy as np
from scipy.optimize import OptimizeResult, differential_evolution, least_squares

#: Optimizers accepted by :func:`identification`
IDENTIFICATION_METHODS = ("differential_evolution", "least_squares")


def identification(
    cost_fn: Callable,
    parameters: list,
    method: str = "differential_evolution",
    *,
    n_starts: int = 1,
    scale="auto",
    n_threads: int = 1,
    **kwargs,
):
    """
    Run parameter identification within the bounds of the ``Parameter`` objects.

    Two optimizers are available:

    - ``method="differential_evolution"`` (default): wraps
      ``scipy.optimize.differential_evolution``. Global and derivative-free,
      but it needs many evaluations (``popsize * n_parameters`` per
      generation) and usually stops on ``maxiter``.
    - ``method="least_squares"``: bounded trust-region reflective least
      squares (``scipy.optimize.least_squares``, ``method="trf"``) from
      ``n_starts`` starting points. ``cost_fn`` returns the vector of
      residuals instead of a scalar (see :func:`calc_residuals`). It converges
      in typically 7 to 20 times fewer evaluations than differential
      evolution on smooth problems, and several starts guard against local
      minima.

    After optimization, the identified values are written back to each
    ``Parameter.value``. For other optimizers (``minimize``, ``shgo``,
    ``dual_annealing``…), call ``scipy.optimize`` directly — the ``Parameter``
    objects expose the bounds and values you need.

    Parameters
    ----------
    cost_fn : callable
        ``f(x, *args)`` where ``x`` is a 1-D array of parameter values (same
        order as *parameters*). Returns a float for
        ``"differential_evolution"``, a 1-D array of residuals of constant
        length for ``"least_squares"`` (the cost minimized is half their sum of
        squares). For ``"least_squares"``, a model failure is best reported by
        residuals set to NaN: the trial step is then rejected and the trust
        region shrunk.
    parameters : list of Parameter
        Parameter objects whose ``.bounds`` define the search space. With
        ``"least_squares"``, the first starting point is ``Parameter.value``
        where it has been set, and the middle of the scaled range otherwise
        (the geometric mean of the bounds for a log-scaled parameter).
    method : {"differential_evolution", "least_squares"}, optional
        Optimizer. Default ``"differential_evolution"``.
    n_starts : int, optional
        ``"least_squares"`` only: number of starting points. The first one is
        described under *parameters*; the others are drawn by Latin
        hypercube sampling inside the (scaled) bounds, reproducibly with
        ``seed``. The best final cost is kept. Default 1.
    scale : {"auto", "linear", "log"} or sequence of them, optional
        ``"least_squares"`` only: how each parameter is mapped onto [0, 1] for
        the optimizer. ``"log"`` (``x = lb (ub/lb)^z``, which requires
        ``lb > 0``) suits parameters spanning decades, such as hardening
        moduli; ``"linear"`` is ``x = lb + z (ub - lb)``. ``"auto"`` (default)
        picks ``"log"`` when the lower bound is strictly positive and
        ``"linear"`` otherwise. One value for all parameters, or one per
        parameter.
    n_threads : int, optional
        Number of threads evaluating ``cost_fn`` concurrently. ``1`` (default)
        is serial; ``0`` uses ``os.cpu_count()``. For
        ``"differential_evolution"``, a generation is evaluated in parallel
        (``workers`` set to a thread pool and ``updating="deferred"``, so the
        trajectory differs from the default serial ``updating="immediate"``,
        but not with the number of threads). For ``"least_squares"``, the
        finite-difference columns of the Jacobian are evaluated in parallel:
        the result does not depend on ``n_threads``. Threads only help when
        ``cost_fn`` spends its time in code releasing the GIL
        (``simcoon.solver.solve`` and ``simcoon.L_eff`` do), and require a
        thread-safe ``cost_fn``.
    **kwargs
        ``"differential_evolution"``: forwarded to
        ``differential_evolution``. Common options: ``maxiter``, ``popsize``,
        ``tol``, ``seed``, ``polish`` (default ``True``), ``disp``, ``args``.

        ``"least_squares"``: ``args`` (extra positional arguments of
        ``cost_fn``), ``seed`` (sampling of the extra starts), ``diff_step``
        (relative finite-difference step in the scaled [0, 1] coordinates,
        default ``sqrt(eps)``; increase it, e.g. to ``1e-6``, when the model
        response is computed to a tolerance), and any option of
        ``least_squares`` but ``x0``, ``bounds``, ``method`` and ``jac``
        (``xtol``, ``ftol``, ``gtol``, ``max_nfev``, ``loss``, ``verbose``…).

    Returns
    -------
    scipy.optimize.OptimizeResult
        The scipy result object; ``result.x`` contains the identified
        parameter values. For ``"least_squares"`` it is the result of the best
        start with ``x`` mapped back to the parameter values, and:

        - ``fun``: the residuals at the solution, ``cost``: half their sum of
          squares (scipy convention);
        - ``z``: the solution in the scaled [0, 1] coordinates;
        - ``nfev``: evaluations of ``cost_fn`` over all starts, Jacobians
          included;
        - ``starts``: one ``OptimizeResult`` per start (``x0``, ``x``,
          ``cost``, ``nfev``, ``status``, ``message``).

    Examples
    --------
    >>> from simcoon.parameter import Parameter
    >>> from simcoon.identify import identification
    >>> params = [Parameter(0, bounds=(10000, 200000))]
    >>> result = identification(my_cost, params, seed=42)
    >>> print(f"Identified E = {params[0].value:.0f}")

    The least-squares route with residuals, three starts:

    >>> def residuals(x):
    ...     return calc_residuals(y_exp, model(x), metric="nmse_per_response")
    >>> result = identification(residuals, params, method="least_squares",
    ...                         n_starts=3, seed=0)

    ``my_cost`` receives the trial values and runs the model in memory
    (:func:`simcoon.solver.solve` takes ``props`` directly). ``key`` and
    ``sim_input_files`` are only needed to drive an external, file-based
    solver by substituting the keys into its input decks.
    """
    if method not in IDENTIFICATION_METHODS:
        raise ValueError(
            f"Unknown method '{method}'. Available: {', '.join(IDENTIFICATION_METHODS)}"
        )
    n_threads = int(n_threads)
    if n_threads < 0:
        raise ValueError(f"n_threads must be >= 0, got {n_threads}")
    if n_threads == 0:
        n_threads = os.cpu_count() or 1

    if method == "least_squares":
        result = _identification_least_squares(
            cost_fn, parameters, n_starts, scale, n_threads, kwargs,
        )
    else:
        bounds = [p.bounds for p in parameters]

        # Sensible defaults
        kwargs.setdefault("polish", True)

        if n_threads > 1:
            if kwargs.get("workers", 1) != 1:
                raise ValueError("Pass either n_threads or workers, not both")
            kwargs["updating"] = "deferred"
            with ThreadPoolExecutor(max_workers=n_threads) as pool:
                kwargs["workers"] = pool.map
                result = differential_evolution(cost_fn, bounds=bounds, **kwargs)
        else:
            result = differential_evolution(cost_fn, bounds=bounds, **kwargs)

    # Write back identified values
    for p, val in zip(parameters, result.x):
        p.value = val

    return result


class _ScaledBounds:
    """Map between parameter values x in [lb, ub] and scaled z in [0, 1]."""

    def __init__(self, parameters, scale):
        self.lb = np.array([p.bounds[0] for p in parameters], dtype=float)
        self.ub = np.array([p.bounds[1] for p in parameters], dtype=float)
        n = len(parameters)
        if n == 0:
            raise ValueError("At least one Parameter is required")
        if not np.all(self.ub > self.lb):
            raise ValueError("Every Parameter needs bounds (lb, ub) with lb < ub")
        if isinstance(scale, str):
            scale = [scale] * n
        if len(scale) != n:
            raise ValueError(f"scale has {len(scale)} entries for {n} parameters")
        self.log = np.zeros(n, dtype=bool)
        for i, s in enumerate(scale):
            if s == "auto":
                self.log[i] = self.lb[i] > 0.0
            elif s == "log":
                if self.lb[i] <= 0.0:
                    raise ValueError(
                        f"scale='log' needs a positive lower bound (parameter {i}: "
                        f"bounds {tuple(parameters[i].bounds)})"
                    )
                self.log[i] = True
            elif s != "linear":
                raise ValueError(f"Unknown scale '{s}': use 'auto', 'linear' or 'log'")
        self._a = np.where(self.log, np.log(np.where(self.log, self.lb, 1.0)), self.lb)
        self._b = np.where(self.log, np.log(np.where(self.log, self.ub, 1.0)), self.ub)

    def to_x(self, z):
        u = self._a + np.asarray(z, dtype=float) * (self._b - self._a)
        x = np.where(self.log, np.exp(u), u)
        # exp/log round trips must not leave the bounds by an ulp
        return np.clip(x, self.lb, self.ub)

    def to_z(self, x):
        x = np.clip(np.asarray(x, dtype=float), self.lb, self.ub)
        u = np.where(self.log, np.log(np.where(self.log, x, 1.0)), x)
        return np.clip((u - self._a) / (self._b - self._a), 0.0, 1.0)


def _identification_least_squares(cost_fn, parameters, n_starts, scale, n_threads, kwargs):
    """Multi-start bounded TRF least squares in scaled coordinates."""
    for key in ("x0", "bounds", "method", "jac"):
        if key in kwargs:
            raise TypeError(f"identification(method='least_squares') does not accept '{key}'")
    args = tuple(kwargs.pop("args", ()))
    seed = kwargs.pop("seed", kwargs.pop("rng", None))
    rel_step = kwargs.pop("diff_step", None)
    rel_step = np.sqrt(np.finfo(float).eps) if rel_step is None else float(rel_step)
    kwargs.setdefault("x_scale", 1.0)
    n_starts = int(n_starts)
    if n_starts < 1:
        raise ValueError(f"n_starts must be >= 1, got {n_starts}")

    sb = _ScaledBounds(parameters, scale)
    n = sb.lb.size
    lock = threading.Lock()
    counter = [0]

    def residuals(z):
        with lock:
            counter[0] += 1
        r = np.atleast_1d(np.asarray(cost_fn(sb.to_x(z), *args), dtype=float))
        if r.ndim != 1:
            raise ValueError("cost_fn must return a 1-D array of residuals")
        return r

    last = {}   # residuals at the last point least_squares evaluated (reused by jac)

    def fun(z):
        r = residuals(z)
        last["z"], last["r"] = z.copy(), r
        return r

    pool = ThreadPoolExecutor(max_workers=n_threads) if n_threads > 1 else None
    mapper = pool.map if pool is not None else map

    def jac(z):
        # Forward differences in the scaled coordinates, one-sided away from
        # the upper bound (as scipy's '2-point' scheme), columns evaluated
        # concurrently when a pool is given. The same arithmetic in both
        # cases, so the result does not depend on n_threads.
        if "z" in last and np.array_equal(last["z"], z):
            r0 = last["r"]
        else:
            r0 = residuals(z)
        h = np.where(z + rel_step <= 1.0, rel_step, -rel_step)
        trial = []
        for i in range(n):
            zi = z.copy()
            zi[i] = z[i] + h[i]
            trial.append(zi)
        cols = list(mapper(residuals, trial))
        J = np.empty((r0.size, n))
        for i, (zi, ri) in enumerate(zip(trial, cols)):
            dz = zi[i] - z[i]   # the step actually representable
            if not np.all(np.isfinite(ri)):
                # model failure on that side: try the other one, else freeze
                zi = z.copy()
                zi[i] = z[i] - h[i]
                ri = residuals(zi) if 0.0 <= zi[i] <= 1.0 else ri
                dz = zi[i] - z[i]
                if not np.all(np.isfinite(ri)):
                    J[:, i] = 0.0
                    continue
            J[:, i] = (ri - r0) / dz
        return J

    # First start: the values set on the Parameters, the middle of the scaled
    # range for those never set (the geometric mean of the bounds in log scale)
    z_set = sb.to_z([p.value for p in parameters])
    starts = [np.array([0.5 if getattr(p, "_value", None) is None else z
                        for p, z in zip(parameters, z_set)])]
    if n_starts > 1:
        # Latin hypercube: one point per stratum of [0.05, 0.95] in each direction
        rng = np.random.default_rng(seed)
        m = n_starts - 1
        strata = np.argsort(rng.random((m, n)), axis=0)
        starts += list(0.05 + 0.9 * (strata + rng.random((m, n))) / m)

    best = None
    runs = []
    try:
        for z0 in starts:
            n0 = counter[0]
            last.clear()
            try:
                res = least_squares(fun, z0, jac=jac, bounds=(0.0, 1.0), method="trf",
                                    **kwargs)
            except ValueError as e:
                if "not finite" not in str(e):
                    raise
                # cost_fn failed at this starting point: try the next one
                runs.append(OptimizeResult(x0=sb.to_x(z0), x=None, cost=np.inf,
                                           nfev=counter[0] - n0, status=-1,
                                           message=str(e), success=False))
                continue
            runs.append(OptimizeResult(x0=sb.to_x(z0), x=sb.to_x(res.x), cost=res.cost,
                                       nfev=counter[0] - n0, status=res.status,
                                       message=res.message, success=res.success))
            if best is None or res.cost < best.cost:
                best = res
    finally:
        if pool is not None:
            pool.shutdown()

    if best is None:
        raise RuntimeError(
            "least_squares failed from every starting point: cost_fn returned "
            "non-finite residuals there"
        )
    best.z = best.x
    best.x = sb.to_x(best.z)
    best.nfev = counter[0]
    best.starts = runs
    return best


def calc_residuals(
    y_exp: List[np.ndarray],
    y_num: List[np.ndarray],
    w_test: Optional[np.ndarray] = None,
    w_response: Optional[List[np.ndarray]] = None,
    w_point: Optional[List[np.ndarray]] = None,
    metric: str = "mse",
) -> np.ndarray:
    """
    Weighted residual vector whose sum of squares is :func:`calc_cost`.

    For the least-squares route of :func:`identification`, which minimizes a
    sum of squared residuals instead of a scalar cost:
    ``np.sum(calc_residuals(...) ** 2) == calc_cost(...)`` with the same
    arguments, for the three metrics that are sums of squares: ``"mse"``,
    ``"nmse"`` and ``"nmse_per_response"``.

    Parameters
    ----------
    y_exp, y_num, w_test, w_response, w_point, metric
        As in :func:`calc_cost`; *metric* must be ``"mse"``, ``"nmse"`` or
        ``"nmse_per_response"``.

    Returns
    -------
    numpy.ndarray
        1-D residuals ``sqrt(weight) * (y_num - y_exp)``, tests concatenated
        in order, each test flattened row by row (point by point).
    """
    if metric not in ("mse", "nmse", "nmse_per_response"):
        raise ValueError(
            f"calc_residuals supports the sum-of-squares metrics 'mse', 'nmse' "
            f"and 'nmse_per_response', not '{metric}'"
        )
    n_tests = len(y_exp)
    if len(y_num) != n_tests:
        raise ValueError(
            f"y_exp has {n_tests} tests but y_num has {len(y_num)}"
        )

    exp_2d, res_2d, w_2d = [], [], []
    for i in range(n_tests):
        exp_i = np.asarray(y_exp[i], dtype=float)
        num_i = np.asarray(y_num[i], dtype=float)
        if exp_i.ndim == 1:
            exp_i = exp_i.reshape(-1, 1)
        if num_i.ndim == 1:
            num_i = num_i.reshape(-1, 1)
        if num_i.shape != exp_i.shape:
            raise ValueError(
                f"Test {i}: y_exp shape {exp_i.shape} != "
                f"y_num shape {num_i.shape}"
            )
        W_i = np.ones_like(exp_i)
        if w_test is not None:
            W_i *= w_test[i]
        if w_response is not None:
            W_i *= np.asarray(w_response[i], dtype=float)[None, :]
        if w_point is not None:
            W_i *= np.abs(w_point[i])
        exp_2d.append(exp_i)
        res_2d.append(num_i - exp_i)
        w_2d.append(W_i)

    if metric == "nmse_per_response":
        # mean over response columns k of sum_k(w r^2) / sum_k(y_exp^2)
        n_cols = max(e.shape[1] for e in exp_2d)
        denom = np.zeros(n_cols)
        for e in exp_2d:
            denom[: e.shape[1]] += np.sum(e ** 2, axis=0)
        denom = np.where(denom > 1e-30, denom, 1.0)
        return np.concatenate([
            (np.sqrt(w / (n_cols * denom[: r.shape[1]])) * r).ravel()
            for r, w in zip(res_2d, w_2d)
        ])

    w_flat = np.concatenate([w.ravel() for w in w_2d])
    r = np.sqrt(w_flat / np.sum(w_flat)) * np.concatenate([r.ravel() for r in res_2d])
    if metric == "nmse":
        e_flat = np.concatenate([e.ravel() for e in exp_2d])
        var = np.average((e_flat - np.average(e_flat, weights=w_flat)) ** 2, weights=w_flat)
        if var > 1e-30:
            r = r / np.sqrt(var)
    return r


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
          then averaged.  Column *k* gives ``sum_tests SSE_k /
          sum_tests sum(y_exp[:, k]²)``: the tests are pooled within a
          column, so the responses of different magnitudes are balanced
          (e.g., force in N vs displacement in mm), but the tests of one
          column are not (a test with larger values weighs more; use
          *w_test* to rebalance them).
        - ``"rmse"`` — Root Mean Squared Error
        - ``"mae"`` — Mean Absolute Error
        - ``"mape"`` — Mean Absolute Percentage Error (same definition as
          ``sklearn.metrics.mean_absolute_percentage_error``: ``|y_exp|`` is
          floored at machine epsilon, so responses crossing zero blow it up —
          prefer ``"wmape"`` for stress/strain histories)
        - ``"wmape"`` — weighted MAPE, ``sum(w |y_exp - y_num|) / sum(w |y_exp|)``

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
    50.8

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
        n_resp = exp_i.shape[1]

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
    """Compute NMSE independently per response column, then average.

    The tests are pooled within a column: NMSE_k = sum_i SSE_ik / sum_i sum(y_exp_ik**2),
    not the mean over the tests of each test's own NMSE.
    """
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


#: Metrics computed with numpy only (also the differentiable set of :func:`simcoon.ml.torch_cost`)
BUILTIN_METRICS = ("mse", "nmse", "nmse_per_response", "rmse", "mae", "mape", "wmape")


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

    if metric == "mape":
        denom = np.maximum(np.abs(y_exp), np.finfo(np.float64).eps)
        return float(np.average(np.abs(residuals) / denom, weights=w))

    if metric == "wmape":
        num = float(np.sum(w * np.abs(residuals)))
        den = float(np.sum(w * np.abs(y_exp)))
        return num / den if den > 1e-30 else num

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
        f"Unknown metric '{metric}'. Built-in: {', '.join(BUILTIN_METRICS)}. "
        f"With scikit-learn: r2, mean_squared_error, mean_absolute_error, etc."
    )
