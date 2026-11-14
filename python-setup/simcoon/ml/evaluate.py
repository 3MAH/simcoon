"""Evaluation of a trained :class:`~simcoon.ml.cells.StateModel` cell with the identification metrics.

Predictions are de-standardised (physical units) and scored sequence by sequence
with :func:`simcoon.identify.calc_cost` — the same function, weights and metric
names as the identification workflow (``mse``, ``nmse``, ``nmse_per_response``,
``rmse``, ``mae``, ``mape``, ``wmape``, ``r2`` and any ``sklearn.metrics`` name).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

from ..identify import calc_cost
from .data import SequenceDataset
from .cells import StateModel

DEFAULT_METRICS = ("mse", "nmse", "r2", "mape", "wmape")


@torch.no_grad()
def predict(model: StateModel, ds: SequenceDataset, device=None, batch_size: int = 256
            ) -> np.ndarray:
    """Predicted stress ``(N, T, n_out)`` in physical units (zeros on padded steps)."""
    dev = torch.device(device) if device is not None else next(model.parameters()).device
    model.eval()
    outs = []
    for i in range(0, len(ds), batch_size):
        y, _, _ = model(ds.x[i:i + batch_size].to(dev))
        outs.append(y.cpu().numpy().astype(float))
    y_pred = np.concatenate(outs, axis=0)
    y_pred[~ds.mask.cpu().numpy()] = 0.0
    return y_pred


def _score(y_exp: List[np.ndarray], y_num: List[np.ndarray], metric: str) -> Optional[float]:
    try:
        return calc_cost(y_exp, y_num, metric=metric)
    except ImportError:          # scikit-learn metric requested without scikit-learn
        return None


def evaluate(
    model: StateModel,
    ds: SequenceDataset,
    metrics: Sequence[str] = DEFAULT_METRICS,
    per_component: bool = True,
    per_sequence: bool = False,
    device=None,
) -> Dict:
    """Score a model on a dataset with :func:`simcoon.identify.calc_cost`.

    Returns
    -------
    dict
        ``{metric: value}`` over all sequences/components, plus
        ``"per_component": {name: {metric: value}}`` and, on request,
        ``"per_sequence": {metric: ndarray (N,)}`` (e.g. for NMSE percentiles).
        Metrics needing scikit-learn are ``None`` when it is not installed.
    """
    y_exp = ds.as_lists()
    y_num = ds.as_lists(predict(model, ds, device=device))
    out: Dict = {m: _score(y_exp, y_num, m) for m in metrics}
    names = tuple(ds.meta.get("stress_components") or model.stress_components)
    if per_component:
        out["per_component"] = {}
        for k, name in enumerate(names):
            exp_k = [y[:, [k]] for y in y_exp]
            num_k = [y[:, [k]] for y in y_num]
            out["per_component"][name] = {m: _score(exp_k, num_k, m) for m in metrics}
    if per_sequence:
        out["per_sequence"] = {
            m: np.array([_score([y_exp[i]], [y_num[i]], m) for i in range(len(ds))], dtype=float)
            for m in metrics
        }
    return out
