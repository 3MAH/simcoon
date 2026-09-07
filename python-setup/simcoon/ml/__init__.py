"""Machine-learning constitutive models (PyTorch).

Two recurrent constitutive cells share one interface
(:class:`~simcoon.ml.cells.StateModel`) and one UMAT wrapper:

* :class:`StressLSTM`, the stress LSTM of Danoun, Prulière and Chemisky
  (Mech. Mater. 2022, thermodynamically consistent RNN; CMAME 2024, FE-LSTM):
  a gated network whose hidden state plays the role of the internal variables;
* :class:`LMSC`, the linearized minimal state cell of Bonatti and Mohr
  (JMPS 2022), whose update is stationary and self-consistent by construction
  and which carries a closed-form algorithmic tangent.

The subpackage also provides
* :mod:`~simcoon.ml.data` — random non-proportional strain paths, dataset generation
  with the simcoon material-point solver, StressLSTM CSV loader,
* :func:`torch_cost` / :func:`train` / :func:`evaluate` — training and evaluation with the
  identification metrics of :mod:`simcoon.identify` (MSE, NMSE, R², MAPE, wMAPE, ...),
* :class:`StressLSTMRegressor` — a scikit-learn compatible estimator,
* :class:`LSTMLaw` — the trained network as a :class:`simcoon.PythonUMAT`, usable by the
  simcoon solver (``simcoon.solver.solve(blocks, law)``), by ``sim.umat("PYEXT", ...)``
  and, through its batched step, by finite-element codes such as fedoo.

PyTorch is an optional dependency: ``pip install simcoon[ml]``.
"""

try:
    import torch  # noqa: F401
except ImportError as exc:  # pragma: no cover - exercised only without torch
    raise ImportError(
        "simcoon.ml requires PyTorch. In a conda environment install it from conda-forge so "
        "that it shares the environment's OpenMP runtime:\n"
        "  conda install -c conda-forge pytorch scikit-learn\n"
        "otherwise (all-PyPI setup):\n"
        "  pip install simcoon[ml]"
    ) from exc

from .cells import VOIGT, StateModel
from .lstm import StressLSTM
from .lmsc import LMSC
from .data import (
    SequenceDataset,
    generate_dataset,
    load_csv,
    mode_components,
    random_strain_paths,
    split_dataset,
)
from .losses import torch_cost
from .tangent import sequence_tangent
from .train import train
from .evaluate import evaluate
from .law import LSTMLaw, RecurrentLaw

__all__ = [
    "VOIGT", "StateModel", "StressLSTM", "LMSC",
    "SequenceDataset", "generate_dataset", "load_csv", "mode_components",
    "random_strain_paths", "split_dataset",
    "torch_cost", "sequence_tangent",
    "train", "evaluate",
    "LSTMLaw", "RecurrentLaw", "StressLSTMRegressor",
]


def __getattr__(name):
    # scikit-learn is only needed by the estimator wrapper: import it on first use
    if name == "StressLSTMRegressor":
        from .sklearn_api import StressLSTMRegressor
        return StressLSTMRegressor
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
