"""scikit-learn compatible estimator around :class:`~simcoon.ml.StressLSTM`.

:class:`StressLSTMRegressor` exposes ``fit`` / ``predict`` / ``score`` /
``get_params`` / ``set_params`` so that the stress LSTM plugs into
``sklearn.model_selection`` (``GridSearchCV``, ``cross_val_score``) and
``sklearn.pipeline``. scikit-learn itself is optional: without it the class still
offers ``fit``/``predict``/``score`` with a minimal parameter interface.

Samples are whole sequences: ``X`` has shape ``(n_sequences, T, n_in)`` and
``Y`` ``(n_sequences, T, n_out)`` (zero-padded; pass ``lengths`` for ragged data).
"""

from __future__ import annotations

import inspect
from typing import Optional, Sequence

import numpy as np

try:  # optional
    from sklearn.base import BaseEstimator, RegressorMixin
except ImportError:  # pragma: no cover - minimal stand-ins
    class BaseEstimator:  # type: ignore[no-redef]
        def get_params(self, deep=True):
            return {n: getattr(self, n) for n in inspect.signature(self.__init__).parameters}

        def set_params(self, **params):
            for k, v in params.items():
                setattr(self, k, v)
            return self

    class RegressorMixin:  # type: ignore[no-redef]
        pass

from .data import SequenceDataset
from .evaluate import DEFAULT_METRICS, evaluate, predict
from .law import LSTMLaw
from .cells import VOIGT
from .lstm import StressLSTM
from .train import train


class StressLSTMRegressor(BaseEstimator, RegressorMixin):
    """Stress LSTM as a scikit-learn regressor on sequences.

    Parameters mirror :class:`StressLSTM` (architecture) and :func:`train`
    (optimisation). ``scoring`` is the metric returned by :meth:`score`
    (``"r2"`` as is, any other :func:`simcoon.identify.calc_cost` metric negated so
    that higher is always better, as scikit-learn expects).

    Attributes (after ``fit``)
    --------------------------
    model_ : StressLSTM
    train_losses_, val_losses_ : list of float
    """

    def __init__(
        self,
        components: Sequence[str] = VOIGT,
        stress_components: Optional[Sequence[str]] = None,
        features: Sequence[str] = ("strain",),
        hidden_size: int = 64,
        num_layers: int = 2,
        epochs: int = 200,
        batch_size: int = 64,
        lr: float = 1e-3,
        loss: str = "mse",
        w_response: Optional[Sequence[float]] = None,
        device: str = "auto",
        random_state: Optional[int] = None,
        verbose: bool = False,
        scoring: str = "r2",
    ):
        self.components = components
        self.stress_components = stress_components
        self.features = features
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.epochs = epochs
        self.batch_size = batch_size
        self.lr = lr
        self.loss = loss
        self.w_response = w_response
        self.device = device
        self.random_state = random_state
        self.verbose = verbose
        self.scoring = scoring

    # ------------------------------------------------------------------ data
    @staticmethod
    def _dataset(X, Y=None, lengths=None) -> SequenceDataset:
        X = np.asarray(X, dtype=float)
        if X.ndim != 3:
            raise ValueError(f"X must be (n_sequences, T, n_in), got {X.shape}")
        if Y is None:
            Y = np.zeros(X.shape[:2] + (1,))          # predict(): targets unused
        Y = np.asarray(Y, dtype=float)
        if Y.ndim == 2:
            Y = Y[..., None]
        n = np.full(len(X), X.shape[1]) if lengths is None else np.asarray(lengths, dtype=int)
        return SequenceDataset.from_sequences([X[i, :n[i]] for i in range(len(X))],
                                              [Y[i, :n[i]] for i in range(len(X))])

    # ------------------------------------------------------------------- API
    def fit(self, X, Y, sample_weight=None, lengths=None, X_val=None, Y_val=None):
        if sample_weight is not None:
            raise ValueError("sample_weight is not supported; weight components with w_response")
        ds = self._dataset(X, Y, lengths)
        val = None if X_val is None else self._dataset(X_val, Y_val)
        model = StressLSTM(
            components=self.components, stress_components=self.stress_components,
            features=self.features, hidden_size=self.hidden_size, num_layers=self.num_layers,
        )
        if model.n_in != ds.n_in or model.n_out != ds.n_out:
            raise ValueError(
                f"data have n_in={ds.n_in}, n_out={ds.n_out} but the model expects "
                f"n_in={model.n_in}, n_out={model.n_out} (components/features)")
        self.train_losses_, self.val_losses_ = train(
            model, ds, val, epochs=self.epochs, batch_size=self.batch_size, lr=self.lr,
            loss=self.loss, w_response=self.w_response, device=self.device,
            verbose=self.verbose, seed=self.random_state,
        )
        self.model_ = model
        self.n_features_in_ = ds.n_in
        return self

    def predict(self, X, lengths=None) -> np.ndarray:
        self._check_fitted()
        return predict(self.model_, self._dataset(X, None, lengths))

    def score(self, X, Y, sample_weight=None, lengths=None, metric: Optional[str] = None) -> float:
        """R² (or the negated ``metric``, so that higher is better) over all sequences."""
        self._check_fitted()
        metric = metric or self.scoring
        value = evaluate(self.model_, self._dataset(X, Y, lengths), metrics=(metric,),
                         per_component=False)[metric]
        if value is None:
            raise ImportError(f"metric '{metric}' requires scikit-learn")
        return value if metric == "r2" else -value

    def evaluate(self, X, Y, lengths=None, metrics=DEFAULT_METRICS, **kw) -> dict:
        """Full report (see :func:`simcoon.ml.evaluate`)."""
        self._check_fitted()
        return evaluate(self.model_, self._dataset(X, Y, lengths), metrics=metrics, **kw)

    def law(self, **kwargs) -> LSTMLaw:
        """The fitted model as a simcoon constitutive law (:class:`LSTMLaw`)."""
        self._check_fitted()
        return LSTMLaw(self.model_, **kwargs)

    def _check_fitted(self):
        if not hasattr(self, "model_"):
            raise RuntimeError("StressLSTMRegressor is not fitted yet; call fit(X, Y) first")
