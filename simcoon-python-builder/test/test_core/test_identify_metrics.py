"""Metrics of simcoon.identify.calc_cost shared with simcoon.ml (mape, wmape, ...)."""

import numpy as np
import pytest

from simcoon.identify import calc_cost


def _data(seed=0):
    rng = np.random.default_rng(seed)
    y_exp = [rng.normal(size=(30, 2)) * 100 + 5, rng.normal(size=(20, 2)) * 50]
    y_num = [y + rng.normal(size=y.shape) * 3 for y in y_exp]
    return y_exp, y_num


def test_wmape_definition():
    y_exp, y_num = _data()
    e = np.concatenate([y.ravel() for y in y_exp])
    n = np.concatenate([y.ravel() for y in y_num])
    expected = np.sum(np.abs(e - n)) / np.sum(np.abs(e))
    assert calc_cost(y_exp, y_num, metric="wmape") == pytest.approx(expected, rel=1e-12)


def test_wmape_weights():
    y_exp, y_num = _data()
    w_test = np.array([2.0, 0.5])
    e0, e1 = y_exp[0].ravel(), y_exp[1].ravel()
    n0, n1 = y_num[0].ravel(), y_num[1].ravel()
    num = 2.0 * np.sum(np.abs(e0 - n0)) + 0.5 * np.sum(np.abs(e1 - n1))
    den = 2.0 * np.sum(np.abs(e0)) + 0.5 * np.sum(np.abs(e1))
    assert calc_cost(y_exp, y_num, w_test=w_test, metric="wmape") == pytest.approx(num / den, rel=1e-12)


def test_mape_definition_matches_sklearn():
    y_exp, y_num = _data()
    e = np.concatenate([y.ravel() for y in y_exp])
    n = np.concatenate([y.ravel() for y in y_num])
    manual = np.mean(np.abs(e - n) / np.maximum(np.abs(e), np.finfo(np.float64).eps))
    assert calc_cost(y_exp, y_num, metric="mape") == pytest.approx(manual, rel=1e-12)
    skm = pytest.importorskip("sklearn.metrics")
    assert calc_cost(y_exp, y_num, metric="mape") == pytest.approx(
        skm.mean_absolute_percentage_error(e, n), rel=1e-12)


def test_unknown_metric_lists_builtins():
    y_exp, y_num = _data()
    with pytest.raises((ValueError, ImportError), match="mape|scikit-learn"):
        calc_cost(y_exp, y_num, metric="no_such_metric")
