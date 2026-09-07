"""Common interface of the recurrent constitutive cells.

A *cell* maps a strain history to a stress history through an internal state
vector. :class:`StateModel` is the contract every cell honours so that a single
UMAT wrapper (:class:`simcoon.ml.LSTMLaw`) and a single fedoo law serve all of
them:

* the state is a **flat** vector of ``state_size`` scalars, stored in the
  solver's ``statev`` (a gated cell packs its own ``(h, c)`` in it),
* one increment is :meth:`StateModel.step`, in physical units,
* :meth:`StateModel.analytic_tangent` returns the closed-form
  :math:`\\partial\\sigma/\\partial\\varepsilon` when the cell has one, and ``None``
  when the wrapper should differentiate the step with autograd instead.

Implementations: :class:`simcoon.ml.StressLSTM` (gated, Danoun et al.) and
:class:`simcoon.ml.LMSC` (linearized minimal state cell, Bonatti and Mohr).
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import torch
from torch import Tensor, nn

#: Voigt component names, simcoon order (engineering shear strains).
VOIGT = ("xx", "yy", "zz", "xy", "xz", "yz")

#: Input features that are strain-like vectors (one entry per active component).
_VECTOR_FEATURES = ("strain", "dstrain")
#: Input features that are scalars.
_SCALAR_FEATURES = ("dtime", "temperature")

#: Below this increment norm the loading direction is not defined: a cell whose tangent is
#: direction-dependent has none there, and the UMAT wrapper answers with the elastic predictor.
DIRECTION_MIN = 1e-14

#: Cell classes that :meth:`StateModel.load` can rebuild, filled by ``register_cell``.
_CELL_REGISTRY: Dict[str, type] = {}


def register_cell(klass):
    """Class decorator making a cell loadable by name from a checkpoint."""
    _CELL_REGISTRY[klass.__name__] = klass
    return klass


def voigt_indices(components: Sequence[str]) -> Tuple[int, ...]:
    """Indices in the 6-vector of the given component names."""
    try:
        return tuple(VOIGT.index(c) for c in components)
    except ValueError as exc:
        raise ValueError(f"unknown Voigt component in {tuple(components)}; valid: {VOIGT}") from exc


def n_inputs(features: Sequence[str], n_comp: int) -> int:
    """Width of the input vector for the given features and number of strain components."""
    return sum(n_comp if f in _VECTOR_FEATURES else 1 for f in features)


def concat_features(
    features: Sequence[str],
    strain: Tensor,
    dstrain: Optional[Tensor] = None,
    dtime: Optional[Tensor] = None,
    temperature: Optional[Tensor] = None,
) -> Tensor:
    """Concatenate the input features (physical units) along the last axis.

    The single definition of the input layout, shared by the cells
    (:meth:`StateModel.build_inputs`) and the dataset builders. ``strain``/``dstrain``
    are ``(..., n_comp)``; ``dtime``/``temperature`` are ``(...)`` or ``(..., 1)``.
    A requested optional feature that is ``None`` raises.
    """
    parts = []
    for f in features:
        if f == "strain":
            parts.append(strain)
        elif f == "dstrain":
            if dstrain is None:
                raise ValueError("feature 'dstrain' requested but dstrain is None")
            parts.append(dstrain)
        elif f in _SCALAR_FEATURES:
            v = dtime if f == "dtime" else temperature
            if v is None:
                raise ValueError(f"feature '{f}' requested but {f} is None")
            v = torch.as_tensor(v, dtype=strain.dtype, device=strain.device)
            if v.dim() == strain.dim() - 1:
                v = v.unsqueeze(-1)
            parts.append(v.expand(*strain.shape[:-1], 1))
        else:
            raise ValueError(f"unknown feature '{f}'; valid: {_VECTOR_FEATURES + _SCALAR_FEATURES}")
    return torch.cat(parts, dim=-1)


class StateModel(nn.Module):
    """Base class of the recurrent constitutive cells.

    Parameters
    ----------
    components : sequence of str
        Active strain components (subset of :data:`VOIGT`).
    stress_components : sequence of str, optional
        Predicted stress components (default: same as ``components``).
    features : sequence of str
        Blocks the cell reads from the dataset input vector, in order; always
        starts with ``"strain"``. A cell may declare a feature it does not use in
        its own equations (the LMSC declares ``"dstrain"`` because it needs the
        increment, and ignores the strain block).

    Subclasses implement :attr:`state_size`, :meth:`zero_state`, :meth:`step`
    and :meth:`forward`, and may override :meth:`analytic_tangent`.
    """

    #: Committed-state tolerance the UMAT wrapper uses by default, as a fraction
    #: of the median training increment. 0 = the cell needs no commit rule.
    default_commit_tol: float = 0.0
    #: The tangent depends on the loading direction (rate-independent cell), so it is
    #: undefined at a zero strain increment: the wrapper then returns the elastic
    #: predictor, the way a classical UMAT answers a zero trial increment.
    directional_tangent: bool = False

    def __init__(
        self,
        components: Sequence[str] = VOIGT,
        stress_components: Optional[Sequence[str]] = None,
        features: Sequence[str] = ("strain",),
    ):
        super().__init__()
        components = tuple(components)
        stress_components = tuple(stress_components) if stress_components else components
        features = tuple(features)
        if not features or features[0] != "strain":
            raise ValueError("features must start with 'strain'")
        for f in features:
            if f not in _VECTOR_FEATURES + _SCALAR_FEATURES:
                raise ValueError(f"unknown feature '{f}'; valid: {_VECTOR_FEATURES + _SCALAR_FEATURES}")
        voigt_indices(components)
        voigt_indices(stress_components)

        self.components = components
        self.stress_components = stress_components
        self.features = features
        self.n_comp = len(components)
        self.n_out = len(stress_components)
        self.n_in = n_inputs(features, self.n_comp)

        self.register_buffer("y_mean", torch.zeros(self.n_out))
        self.register_buffer("y_std", torch.ones(self.n_out))
        # median norm of the strain increment between consecutive training steps
        self.register_buffer("median_increment", torch.zeros(()))

    # ------------------------------------------------------------------ meta
    @property
    def hparams(self) -> Dict:
        """Constructor keywords rebuilding this cell (used by :meth:`save`)."""
        raise NotImplementedError

    @property
    def state_size(self) -> int:
        """Number of scalars of the flat state vector stored in ``statev``."""
        raise NotImplementedError

    # ------------------------------------------------------------- features
    def build_inputs(self, strain, dstrain=None, dtime=None, temperature=None) -> Tensor:
        """Input features of this cell in physical units (see :func:`concat_features`)."""
        return concat_features(self.features, strain, dstrain, dtime, temperature)

    def split_features(self, x: Tensor) -> Dict[str, Tensor]:
        """Named blocks of an assembled input vector ``x (..., n_in)``."""
        out, i = {}, 0
        for f in self.features:
            w = self.n_comp if f in _VECTOR_FEATURES else 1
            out[f] = x[..., i:i + w]
            i += w
        return out

    # -------------------------------------------------------------- scalers
    @torch.no_grad()
    def fit_scalers(self, x: Tensor, y: Tensor, mask: Optional[Tensor] = None) -> None:
        """Standardisation buffers and median increment from ``x (N,T,n_in)``, ``y (N,T,n_out)``.

        ``mask (N,T)`` selects the valid (non-padded) steps. Constant responses get
        ``std = 1``.
        """
        m = None if mask is None else mask.bool()
        yf = y.reshape(-1, y.shape[-1]) if m is None else y[m]
        ym, ys = yf.mean(0), yf.std(0, unbiased=False)
        self.y_mean.copy_(ym.to(self.y_mean))
        self.y_std.copy_(torch.where(ys > 0, ys, torch.ones_like(ys)).to(self.y_std))
        strain = x[..., :self.n_comp]
        d = (strain[:, 1:] - strain[:, :-1]).norm(dim=-1)
        valid = (m[:, 1:] & m[:, :-1]) if m is not None else torch.ones_like(d, dtype=torch.bool)
        if bool(valid.any()):
            self.median_increment.copy_(d[valid].median().to(self.median_increment))

    def standardize_y(self, y: Tensor) -> Tensor:
        return (y - self.y_mean) / self.y_std

    def destandardize_y(self, ys: Tensor) -> Tensor:
        return ys * self.y_std + self.y_mean

    # -------------------------------------------------------------- dynamics
    def zero_state(self, batch: int, dtype=None, device=None) -> Tensor:
        """Initial (virgin) state, ``(batch, state_size)``."""
        p = next(self.parameters())
        return torch.zeros(batch, self.state_size, dtype=dtype or p.dtype,
                           device=device or p.device)

    def step(self, x: Tensor, s: Tensor) -> Tuple[Tensor, Tensor]:
        """One increment: ``x (B, n_in)`` assembled features (physical units, see
        :meth:`build_inputs`) and ``s (B, state_size)`` -> ``stress (B, n_out)``,
        ``s' (B, state_size)``."""
        raise NotImplementedError

    def forward(self, x: Tensor, s: Optional[Tensor] = None):
        """Run a whole sequence: ``x (B, T, n_in)`` -> ``y (B, T, n_out)``, ``s (B, state_size)``."""
        raise NotImplementedError

    def analytic_tangent(self, x: Tensor, s: Tensor) -> Optional[Tensor]:
        """Closed-form ``d sigma / d strain`` of :meth:`step`, ``(B, n_out, n_comp)``.

        The derivative is taken with respect to the strain at the end of the increment,
        which is the quantity the solver's Newton loop varies.

        ``None`` (the default) tells the UMAT wrapper to differentiate the step with
        autograd instead.
        """
        return None

    # --------------------------------------------------------------- persist
    def save(self, path) -> None:
        """Save class name, hyper-parameters, weights and scalers (``torch.save``)."""
        torch.save({"class": type(self).__name__, "hparams": self.hparams,
                    "state_dict": self.state_dict()}, path)

    @classmethod
    def load(cls, path, map_location="cpu"):
        """Rebuild a cell saved by :meth:`save` (any registered cell class)."""
        blob = torch.load(path, map_location=map_location, weights_only=True)
        name = blob.get("class", "StressLSTM")
        klass = _CELL_REGISTRY.get(name)
        if klass is None:
            raise ValueError(f"unknown model class '{name}' in {path}; "
                             f"known: {sorted(_CELL_REGISTRY)}")
        if cls is not StateModel and not issubclass(klass, cls):
            raise TypeError(f"{path} holds a {name}, not a {cls.__name__}")
        model = klass(**blob["hparams"])
        # strict=False: models saved before a buffer was added keep that buffer's default
        model.load_state_dict(blob["state_dict"], strict=False)
        model.eval()
        return model
