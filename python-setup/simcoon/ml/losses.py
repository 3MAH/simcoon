"""Differentiable (torch) transcription of the identification cost of :mod:`simcoon.identify`.

:func:`torch_cost` has the same three-level weighting as
:func:`simcoon.identify.calc_cost` (``w_test`` per sequence, ``w_response`` per
component, ``w_point`` per time step) and the same metric definitions, so that
a network trained with ``torch_cost(..., metric="nmse")`` is evaluated with
``calc_cost(..., metric="nmse")`` on identical numbers. The padding mask of a
batch of sequences is simply a ``w_point`` equal to zero on padded steps.
"""

from __future__ import annotations

from typing import Optional

import torch
from torch import Tensor

from ..identify import BUILTIN_METRICS

#: Metrics available for training (differentiable almost everywhere) = the numpy built-ins.
TRAIN_METRICS = BUILTIN_METRICS


def _as3d(y: Tensor) -> Tensor:
    if y.dim() == 2:
        return y.unsqueeze(-1)
    if y.dim() != 3:
        raise ValueError(f"expected (N, T, R) or (N, T) tensors, got shape {tuple(y.shape)}")
    return y


def combined_weights(
    y_true: Tensor,
    w_test: Optional[Tensor] = None,
    w_response: Optional[Tensor] = None,
    w_point: Optional[Tensor] = None,
) -> Tensor:
    """Multiplicative weight tensor ``W (N, T, R)`` (``calc_cost`` semantics).

    ``w_test``: ``(N,)``; ``w_response``: ``(R,)`` or ``(N, R)``; ``w_point``: ``(N, T)``
    or ``(N, T, R)`` (absolute value taken, as in ``calc_cost``). The result is an
    expanded view of the supplied factors (no full-size allocation when possible).
    """
    N, T, R = y_true.shape
    kw = dict(dtype=y_true.dtype, device=y_true.device)
    factors = []
    if w_test is not None:
        factors.append(torch.as_tensor(w_test, **kw).reshape(N, 1, 1))
    if w_response is not None:
        wr = torch.as_tensor(w_response, **kw)
        factors.append(wr.reshape(1, 1, R) if wr.dim() == 1 else wr.reshape(N, 1, R))
    if w_point is not None:
        wp = torch.as_tensor(w_point, **kw).abs()
        factors.append(wp.reshape(N, T, 1) if wp.dim() == 2 else wp.reshape(N, T, R))
    if not factors:
        return torch.ones_like(y_true)
    W = factors[0]
    for f in factors[1:]:
        W = W * f
    return W.expand(N, T, R)


def torch_cost(
    y_true: Tensor,
    y_pred: Tensor,
    w_test: Optional[Tensor] = None,
    w_response: Optional[Tensor] = None,
    w_point: Optional[Tensor] = None,
    metric: str = "mse",
) -> Tensor:
    """Weighted cost between reference and predicted sequences (differentiable).

    Parameters
    ----------
    y_true, y_pred : Tensor (N, T, R)
        Reference and predicted sequences (``N`` sequences, ``T`` steps, ``R`` responses).
    w_test, w_response, w_point
        Weights per sequence, per response and per point (see :func:`combined_weights`).
        A padding mask goes in ``w_point``.
    metric : str
        One of :data:`TRAIN_METRICS`: ``"mse"``, ``"nmse"`` (MSE / weighted variance of
        ``y_true``), ``"nmse_per_response"`` (per-column weighted SSE / unweighted
        ``sum(y_true**2)``, averaged over columns — exactly
        :func:`simcoon.identify.calc_cost`), ``"rmse"``, ``"mae"``, ``"mape"``, ``"wmape"``.

    Returns
    -------
    Tensor
        Scalar cost.
    """
    y_true, y_pred = _as3d(y_true), _as3d(y_pred)
    if y_true.shape != y_pred.shape:
        raise ValueError(f"shape mismatch: y_true {tuple(y_true.shape)} vs y_pred {tuple(y_pred.shape)}")
    W = combined_weights(y_true, w_test, w_response, w_point)
    r = y_true - y_pred
    tiny = torch.as_tensor(1e-30, dtype=y_true.dtype, device=y_true.device)

    if metric == "nmse_per_response":
        sse = (W * r ** 2).sum(dim=(0, 1))
        denom = (y_true ** 2).sum(dim=(0, 1))          # unweighted, as in calc_cost
        nmse_k = torch.where(denom > tiny, sse / torch.clamp(denom, min=1e-300), sse)
        return nmse_k.mean()

    sw = W.sum()
    if metric == "mse":
        return (W * r ** 2).sum() / sw
    if metric == "rmse":
        return torch.sqrt((W * r ** 2).sum() / sw)
    if metric == "mae":
        return (W * r.abs()).sum() / sw
    if metric == "nmse":
        mse = (W * r ** 2).sum() / sw
        mean = (W * y_true).sum() / sw
        var = (W * (y_true - mean) ** 2).sum() / sw
        return torch.where(var > tiny, mse / torch.clamp(var, min=1e-300), mse)
    if metric == "mape":
        eps = torch.finfo(torch.float64).eps
        return (W * r.abs() / torch.clamp(y_true.abs(), min=eps)).sum() / sw
    if metric == "wmape":
        num = (W * r.abs()).sum()
        den = (W * y_true.abs()).sum()
        return torch.where(den > tiny, num / torch.clamp(den, min=1e-300), num)
    raise ValueError(f"Unknown metric '{metric}'. Available: {TRAIN_METRICS}")
