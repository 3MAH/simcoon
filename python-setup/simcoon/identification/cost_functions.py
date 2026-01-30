"""
Cost functions for parameter identification.

This module provides standard cost function metrics for comparing
simulation predictions with experimental data.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional, Callable

import numpy as np
import numpy.typing as npt


def mse(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
    weights: Optional[npt.NDArray[np.float64]] = None,
) -> float:
    """
    Mean Squared Error.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth (experimental) values
    y_pred : np.ndarray
        Predicted (simulation) values
    weights : np.ndarray, optional
        Sample weights

    Returns
    -------
    float
        MSE value

    Examples
    --------
    >>> y_true = np.array([1.0, 2.0, 3.0])
    >>> y_pred = np.array([1.1, 2.0, 2.9])
    >>> mse(y_true, y_pred)
    0.006666...
    """
    if weights is None:
        return float(np.mean((y_true - y_pred) ** 2))
    else:
        return float(np.average((y_true - y_pred) ** 2, weights=weights))


def mae(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
    weights: Optional[npt.NDArray[np.float64]] = None,
) -> float:
    """
    Mean Absolute Error.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values
    weights : np.ndarray, optional
        Sample weights

    Returns
    -------
    float
        MAE value
    """
    if weights is None:
        return float(np.mean(np.abs(y_true - y_pred)))
    else:
        return float(np.average(np.abs(y_true - y_pred), weights=weights))


def r2(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
) -> float:
    """
    Coefficient of determination (R-squared).

    Returns 1 - R^2 so that it can be minimized.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values

    Returns
    -------
    float
        1 - R^2 value (for minimization)
    """
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot < 1e-10:
        return 0.0
    r2_val = 1.0 - (ss_res / ss_tot)
    return 1.0 - r2_val  # Return 1 - R^2 for minimization


def weighted_mse(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
    weights: npt.NDArray[np.float64],
) -> float:
    """
    Weighted Mean Squared Error.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values
    weights : np.ndarray
        Sample weights (must be provided)

    Returns
    -------
    float
        Weighted MSE value
    """
    return float(np.average((y_true - y_pred) ** 2, weights=weights))


def huber_loss(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
    delta: float = 1.0,
) -> float:
    """
    Huber loss - robust to outliers.

    For residuals smaller than delta, behaves like MSE.
    For larger residuals, behaves like MAE.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values
    delta : float
        Threshold parameter

    Returns
    -------
    float
        Huber loss value
    """
    residual = np.abs(y_true - y_pred)
    quadratic = np.minimum(residual, delta)
    linear = residual - quadratic
    return float(np.mean(0.5 * quadratic ** 2 + delta * linear))


def rmse(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
) -> float:
    """
    Root Mean Squared Error.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values

    Returns
    -------
    float
        RMSE value
    """
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def nrmse(
    y_true: npt.NDArray[np.float64],
    y_pred: npt.NDArray[np.float64],
    normalization: str = 'range',
) -> float:
    """
    Normalized Root Mean Squared Error.

    Parameters
    ----------
    y_true : np.ndarray
        Ground truth values
    y_pred : np.ndarray
        Predicted values
    normalization : str
        'range' (max-min), 'mean', or 'std'

    Returns
    -------
    float
        NRMSE value
    """
    rmse_val = rmse(y_true, y_pred)

    if normalization == 'range':
        norm_factor = np.max(y_true) - np.min(y_true)
    elif normalization == 'mean':
        norm_factor = np.mean(y_true)
    elif normalization == 'std':
        norm_factor = np.std(y_true)
    else:
        raise ValueError(f"Unknown normalization: {normalization}")

    if norm_factor < 1e-10:
        return 0.0

    return rmse_val / norm_factor


class CostFunction(ABC):
    """
    Abstract base class for custom cost functions.

    Subclass this to implement custom cost functions for specific
    identification problems.

    Examples
    --------
    >>> class MyCost(CostFunction):
    ...     def __call__(self, y_true, y_pred):
    ...         return np.mean((y_true - y_pred) ** 2)
    >>>
    >>> cost_fn = MyCost()
    >>> cost_fn(exp_data, sim_data)
    """

    @abstractmethod
    def __call__(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> float:
        """Compute cost between true and predicted values."""
        pass

    def gradient(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """
        Compute gradient of cost with respect to predictions.

        Default implementation uses numerical differentiation.
        Override for analytical gradients.
        """
        eps = 1e-7
        n = len(y_pred)
        grad = np.zeros(n)
        for i in range(n):
            y_pred_plus = y_pred.copy()
            y_pred_plus[i] += eps
            grad[i] = (self(y_true, y_pred_plus) - self(y_true, y_pred)) / eps
        return grad


class MSECost(CostFunction):
    """MSE cost function as a class."""

    def __call__(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> float:
        return mse(y_true, y_pred)

    def gradient(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        n = len(y_true)
        return 2.0 * (y_pred - y_true) / n


class MAECost(CostFunction):
    """MAE cost function as a class."""

    def __call__(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> float:
        return mae(y_true, y_pred)


class HuberCost(CostFunction):
    """Huber cost function as a class."""

    def __init__(self, delta: float = 1.0):
        self.delta = delta

    def __call__(
        self,
        y_true: npt.NDArray[np.float64],
        y_pred: npt.NDArray[np.float64],
    ) -> float:
        return huber_loss(y_true, y_pred, self.delta)
