"""Datasets for the stress LSTM.

* :func:`random_strain_paths` — random non-proportional piecewise-linear strain paths
  (the loading programme of the FE-LSTM / thesis datasets: a few segments towards
  uniformly drawn targets, a fixed number of sub-steps per segment);
* :func:`generate_dataset` — integrate those paths with the simcoon material-point
  solver (any built-in UMAT or :class:`simcoon.PythonUMAT`) and pack the
  strain/stress histories into a :class:`SequenceDataset`;
* :func:`load_csv` — read the StressLSTM CSV format (``simulation_load_id``,
  ``total_strain_*``, ``stress_*``);
* :func:`split_dataset` — train/test split (scikit-learn when available).

Modes follow the ``ndi`` convention of the simcoon kernels (see
:mod:`simcoon.ml.condense`): ``"3D"`` and ``"plane_strain"`` are integrated with
``ndi = 3`` (plane strain = null out-of-plane strains), ``"plane_stress"`` drives the
out-of-plane stress to zero (``ndi = 2``) and ``"uniaxial"`` the lateral stresses
(``ndi = 1``), so the generated data satisfy the constraint the model will be used with.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import torch
from torch.utils.data import Dataset

from .condense import control_from_ndi
from .cells import VOIGT, concat_features, voigt_indices

#: (components, solver control per Voigt component, ndi) for each mode
MODES: Dict[str, Tuple[Tuple[str, ...], Tuple[str, ...], int]] = {
    "3D": (VOIGT, control_from_ndi(3), 3),
    "plane_strain": (("xx", "yy", "xy"), control_from_ndi(3), 3),
    "plane_stress": (("xx", "yy", "xy"), control_from_ndi(2), 2),
    "uniaxial": (("xx",), control_from_ndi(1), 1),
}


def mode_components(mode: str) -> Tuple[Tuple[str, ...], Tuple[str, ...], int]:
    """``(components, control, ndi)`` of a loading mode (see :data:`MODES`)."""
    try:
        return MODES[mode]
    except KeyError as exc:
        raise ValueError(f"unknown mode '{mode}'; valid: {tuple(MODES)}") from exc


# ---------------------------------------------------------------------------
# random loading paths
# ---------------------------------------------------------------------------

def random_strain_paths(
    n_paths: int,
    n_segments: int = 4,
    n_sub: int = 25,
    amplitude: Union[None, float, Sequence[float]] = None,
    components: Sequence[str] = VOIGT,
    seed: Optional[int] = None,
    variable_substeps: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Random piecewise-linear strain paths.

    Each path is a sequence of ``n_segments`` strain targets drawn uniformly in
    ``[-amplitude, amplitude]`` per component, joined by linear ramps of ``n_sub``
    increments (the thesis setting: 4 segments x 25 sub-steps, +/-0.05 on normal and
    +/-0.10 on shear components).

    Parameters
    ----------
    amplitude : float or sequence, optional
        Per-component half-range. Default: 0.05 for normal components, 0.10 for
        engineering shear components.
    variable_substeps : bool
        Draw the number of sub-steps of each segment uniformly in
        ``[max(2, n_sub // 2), 2 * n_sub]`` (data augmentation for variable time steps).

    Returns
    -------
    targets : ndarray (n_paths, n_segments, n_comp)
    ninc : ndarray (n_paths, n_segments) of int
        Number of increments of each segment.
    """
    rng = np.random.default_rng(seed)
    components = tuple(components)
    n_comp = len(components)
    if amplitude is None:
        amp = np.array([0.05 if c in ("xx", "yy", "zz") else 0.10 for c in components])
    else:
        amp = np.broadcast_to(np.asarray(amplitude, dtype=float), (n_comp,)).copy()
    targets = rng.uniform(-1.0, 1.0, size=(n_paths, n_segments, n_comp)) * amp
    if variable_substeps:
        ninc = rng.integers(max(2, n_sub // 2), 2 * n_sub + 1, size=(n_paths, n_segments))
    else:
        ninc = np.full((n_paths, n_segments), int(n_sub), dtype=int)
    return targets, ninc


# ---------------------------------------------------------------------------
# sequence dataset
# ---------------------------------------------------------------------------

class SequenceDataset(Dataset):
    """Padded batch of (input, target) sequences with a validity mask.

    Attributes
    ----------
    x : Tensor (N, T, n_in)
        Input features (physical units), zero-padded.
    y : Tensor (N, T, n_out)
        Target stress, zero-padded.
    mask : BoolTensor (N, T)
        True on valid steps.
    lengths : LongTensor (N,)
    extras : dict of Tensor (N, T, ...)
        Optional per-step quantities recorded with the dataset (``"statev"``,
        ``"tangent"``), padded and sliced like ``x`` and ``y``.
    meta : dict
        ``components``, ``stress_components``, ``features``, ``mode`` when known.
    """

    def __init__(self, x: torch.Tensor, y: torch.Tensor, mask: Optional[torch.Tensor] = None,
                 extras: Optional[Dict[str, torch.Tensor]] = None, **meta):
        if x.dim() != 3 or y.dim() != 3 or x.shape[:2] != y.shape[:2]:
            raise ValueError("x and y must be (N, T, .) tensors with the same N and T")
        self.x = x
        self.y = y
        self.mask = torch.ones(x.shape[:2], dtype=torch.bool) if mask is None else mask.bool()
        self.lengths = self.mask.sum(dim=1)
        #: optional per-step quantities recorded beside (x, y): ``"statev"`` (N, T, nstatev),
        #: ``"tangent"`` (N, T, 6, 6) — see :func:`generate_dataset`.
        self.extras: Dict[str, torch.Tensor] = dict(extras or {})
        self.meta = dict(meta)

    @classmethod
    def from_sequences(cls, xs: Sequence[np.ndarray], ys: Sequence[np.ndarray],
                       dtype=torch.float32, extras: Optional[Dict[str, Sequence[np.ndarray]]] = None,
                       **meta) -> "SequenceDataset":
        """Pad variable-length sequences ``xs[i] (T_i, n_in)``, ``ys[i] (T_i, n_out)``.

        ``extras`` maps a name to one array per sequence, ``(T_i, ...)``, padded the same way.
        """
        n = len(xs)
        if n == 0:
            raise ValueError("no sequences")
        T = max(len(s) for s in xs)
        n_in, n_out = np.shape(xs[0])[1], np.shape(ys[0])[1]
        x = torch.zeros(n, T, n_in, dtype=dtype)
        y = torch.zeros(n, T, n_out, dtype=dtype)
        mask = torch.zeros(n, T, dtype=torch.bool)
        for i, (xi, yi) in enumerate(zip(xs, ys)):
            ti = len(xi)
            if len(yi) != ti:
                raise ValueError(f"sequence {i}: x has {ti} steps, y has {len(yi)}")
            x[i, :ti] = torch.as_tensor(np.asarray(xi), dtype=dtype)
            y[i, :ti] = torch.as_tensor(np.asarray(yi), dtype=dtype)
            mask[i, :ti] = True
        ex = {}
        for name, seqs in (extras or {}).items():
            first = np.asarray(seqs[0])
            buf = torch.zeros((n, T) + first.shape[1:], dtype=dtype)
            for i, v in enumerate(seqs):
                buf[i, :len(v)] = torch.as_tensor(np.asarray(v), dtype=dtype)
            ex[name] = buf
        return cls(x, y, mask, ex, **meta)

    def __len__(self) -> int:
        return self.x.shape[0]

    def __getitem__(self, i):
        return self.x[i], self.y[i], self.mask[i]

    @property
    def n_in(self) -> int:
        return self.x.shape[-1]

    @property
    def n_out(self) -> int:
        return self.y.shape[-1]

    def subset(self, idx) -> "SequenceDataset":
        idx = torch.as_tensor(np.asarray(idx), dtype=torch.long)
        return SequenceDataset(self.x[idx], self.y[idx], self.mask[idx],
                               {k: v[idx] for k, v in self.extras.items()}, **self.meta)

    def to(self, device=None, dtype=None) -> "SequenceDataset":
        x = self.x.to(device=device, dtype=dtype)
        y = self.y.to(device=device, dtype=dtype)
        extras = {k: v.to(device=device, dtype=dtype) for k, v in self.extras.items()}
        return SequenceDataset(x, y, self.mask.to(device=device), extras, **self.meta)

    def as_lists(self, values: Optional[np.ndarray] = None) -> List[np.ndarray]:
        """Unpadded per-sequence arrays, the layout of :func:`simcoon.identify.calc_cost`.

        ``values`` (``(N, T, R)`` array, e.g. predictions) defaults to the targets ``y``.
        """
        src = self.y.detach().cpu().numpy() if values is None else np.asarray(values)
        lengths = self.lengths.cpu().numpy()
        return [src[i, :lengths[i]].astype(float) for i in range(len(self))]


def split_dataset(ds: SequenceDataset, test_size: float = 0.3, seed: int = 0
                  ) -> Tuple[SequenceDataset, SequenceDataset]:
    """Random train/test split of the sequences (``sklearn.model_selection.train_test_split``
    when scikit-learn is installed, a seeded numpy permutation otherwise)."""
    idx = np.arange(len(ds))
    try:
        from sklearn.model_selection import train_test_split
        tr, te = train_test_split(idx, test_size=test_size, random_state=seed)
    except ImportError:
        rng = np.random.default_rng(seed)
        rng.shuffle(idx)
        n_test = int(np.ceil(test_size * len(idx)))      # scikit-learn's convention
        te, tr = idx[:n_test], idx[n_test:]
    return ds.subset(np.sort(tr)), ds.subset(np.sort(te))


# ---------------------------------------------------------------------------
# time resampling (training augmentation against increment-size bias)
# ---------------------------------------------------------------------------

def resample_time(x: torch.Tensor, y: torch.Tensor, mask: torch.Tensor, n_steps: int,
                  features: Sequence[str], n_comp: int,
                  extras: Optional[Dict[str, torch.Tensor]] = None
                  ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, Dict[str, torch.Tensor]]:
    """Rewrite a padded batch of sequences on ``n_steps`` time steps (linear interpolation).

    Every sequence is interpolated over its own valid length onto ``n_steps`` equally
    spaced points (the corners of the loading path are kept as interpolation nodes),
    the derived input features are rebuilt from the resampled strain (``dstrain`` as the
    backward difference, ``dtime`` scaled by the change of step count) and the batch is
    returned unpadded (all sequences have ``n_steps`` valid steps). Used by
    :func:`simcoon.ml.train` to expose the network to many increment sizes, so that it does
    not associate one discretisation with the constitutive response.
    """
    N, T, _ = x.shape
    extras = extras or {}
    lengths = mask.sum(dim=1)
    grid = torch.linspace(0.0, 1.0, n_steps, device=x.device, dtype=x.dtype)
    xs, ys = [], []
    exs: Dict[str, list] = {k: [] for k in extras}
    for i in range(N):
        L = int(lengths[i])
        src = torch.linspace(0.0, 1.0, L, device=x.device, dtype=x.dtype)
        # index of the left node for each target point
        j = torch.clamp(torch.searchsorted(src, grid, right=True) - 1, 0, max(L - 2, 0))
        wt = ((grid - src[j]) / (src[j + 1] - src[j])) if L > 1 else torch.zeros(n_steps, dtype=x.dtype, device=x.device)

        def interp(a, _w=wt, _j=j, _L=L):
            """Linear interpolation along axis 0, whatever the trailing shape."""
            if _L <= 1:
                return a[_j]
            return a[_j] + _w.reshape(-1, *([1] * (a.dim() - 1))) * (a[_j + 1] - a[_j])
        strain = interp(x[i, :L, :n_comp])
        parts = [strain]
        col = n_comp
        for f in features[1:]:
            if f == "dstrain":
                parts.append(torch.diff(strain, dim=0, prepend=torch.zeros_like(strain[:1])))
                col += n_comp
            elif f == "dtime":
                # per-row increments: rescale by the row count so the total duration of the
                # sequence is preserved ((L-1)/(n-1) would stretch or shrink it)
                parts.append(interp(x[i, :L, col:col + 1]) * (L / n_steps))
                col += 1
            else:                                   # temperature: plain interpolation
                parts.append(interp(x[i, :L, col:col + 1]))
                col += 1
        xs.append(torch.cat(parts, dim=-1))
        ys.append(interp(y[i, :L]))
        for k, v in extras.items():
            exs[k].append(interp(v[i, :L]))
    return (torch.stack(xs), torch.stack(ys),
            torch.ones(N, n_steps, dtype=torch.bool, device=x.device),
            {k: torch.stack(v) for k, v in exs.items()})


# ---------------------------------------------------------------------------
# feature assembly (numpy front-end of the model's input layout)
# ---------------------------------------------------------------------------

def assemble_features(features: Sequence[str], strain: np.ndarray,
                      dtime: Optional[np.ndarray] = None,
                      temperature: Optional[np.ndarray] = None) -> np.ndarray:
    """Input array ``(T, n_in)`` of a strain history ``(T, n_comp)`` for the given features
    (same layout as :func:`simcoon.ml.cells.concat_features`); the strain increment
    feature is the backward difference of the history, starting from zero."""
    strain_t = torch.as_tensor(np.asarray(strain, dtype=float))
    dstrain_t = None
    if "dstrain" in features:
        dstrain_t = torch.diff(strain_t, dim=0, prepend=torch.zeros_like(strain_t[:1]))
    x = concat_features(features, strain_t, dstrain_t,
                        None if dtime is None else np.asarray(dtime, dtype=float),
                        None if temperature is None else np.asarray(temperature, dtype=float))
    return x.numpy()


# ---------------------------------------------------------------------------
# generation with the simcoon solver
# ---------------------------------------------------------------------------

def generate_dataset(
    umat,
    props,
    nstatev,
    targets: np.ndarray,
    ninc: Union[int, np.ndarray] = 25,
    mode: str = "3D",
    components: Optional[Sequence[str]] = None,
    control: Optional[Sequence[str]] = None,
    stress_components: Optional[Sequence[str]] = None,
    features: Sequence[str] = ("strain",),
    T_init: float = 293.15,
    step_time: float = 1.0,
    record: Sequence[str] = (),
    dtype=torch.float32,
    progress: bool = False,
    **solve_kw,
) -> SequenceDataset:
    """Integrate strain paths with the simcoon material-point solver into a dataset.

    Parameters
    ----------
    umat, props, nstatev
        Constitutive model as accepted by :func:`simcoon.solver.solve` (built-in name
        with ``props``/``nstatev``, or a :class:`simcoon.PythonUMAT` with ``None``).
    targets : ndarray (n_paths, n_segments, n_comp)
        Strain targets of the active components (see :func:`random_strain_paths`).
    ninc : int or ndarray (n_paths, n_segments)
        Increments per segment.
    mode : str
        ``"3D"``, ``"plane_strain"``, ``"plane_stress"`` or ``"uniaxial"`` — sets the
        active components and the solver control (strain-driven active components,
        zero-stress-driven out-of-plane directions for plane stress / uniaxial).
    components, control, stress_components
        Override the mode defaults (``control`` = 6 entries ``"strain"``/``"stress"``,
        e.g. free lateral faces for a 6-component model).
    features : sequence of str
        Input features stored in ``x`` (see :class:`~simcoon.ml.StressLSTM`).
    step_time : float
        Duration of each segment (sets the ``dtime`` feature).
    record : sequence of str
        Extra per-step quantities to keep in ``dataset.extras``: ``"statev"`` (internal
        variables, e.g. to weight the loss by the plastic increment — see
        :func:`simcoon.ml.weights_from_increment`) and ``"tangent"`` (the reference
        tangent operator, for a tangent-aware loss — see :func:`simcoon.ml.train`).
    **solve_kw
        Forwarded to :func:`simcoon.solver.solve` (e.g. ``tangent_mode``).

    Returns
    -------
    SequenceDataset
        ``x = features``, ``y = stress`` of ``stress_components``; ``meta`` records the
        components, features and mode.
    """
    from ..solver import Block, StepMeca, solve

    m_comp, m_ctrl, ndi = mode_components(mode)
    components = tuple(components) if components else m_comp
    control = tuple(control) if control else m_ctrl
    stress_components = tuple(stress_components) if stress_components else components
    if len(control) != 6:
        raise ValueError("control must have 6 entries")
    idx_in = list(voigt_indices(components))
    idx_out = list(voigt_indices(stress_components))
    targets = np.asarray(targets, dtype=float)
    n_paths, n_seg, n_comp = targets.shape
    if n_comp != len(components):
        raise ValueError(f"targets have {n_comp} components, expected {len(components)}")
    ninc_arr = np.broadcast_to(np.asarray(ninc, dtype=int), (n_paths, n_seg))

    record = tuple(record)
    for name in record:
        if name not in ("statev", "tangent"):
            raise ValueError(f"unknown record '{name}'; valid: 'statev', 'tangent'")
    xs, ys = [], []
    ex: Dict[str, list] = {name: [] for name in record}
    it = range(n_paths)
    if progress:
        try:
            from tqdm import tqdm
            it = tqdm(it, desc="generate_dataset")
        except ImportError:
            pass
    for i in it:
        steps = []
        for k in range(n_seg):
            value = np.zeros(6)
            value[idx_in] = targets[i, k]
            steps.append(StepMeca(control=list(control), value=value, ninc=int(ninc_arr[i, k]),
                                  time=step_time))
        res = solve(Block(steps=steps), umat, props, nstatev, T_init=T_init, **solve_kw)
        time = np.asarray(res["Time"], dtype=float)
        xs.append(assemble_features(features, res["Strain"][idx_in].T,
                                    dtime=np.diff(time, prepend=0.0),
                                    temperature=np.asarray(res["Temp"], dtype=float)))
        ys.append(res["Stress"][idx_out].T)
        if "statev" in ex:
            ex["statev"].append(np.asarray(res["Statev"], dtype=float).T)
        if "tangent" in ex:
            # (6, 6, N) -> (N, n_out, n_comp): the block the Newton loop consumes
            ex["tangent"].append(res["TangentMatrix"][np.ix_(idx_out, idx_in)].transpose(2, 0, 1))
    return SequenceDataset.from_sequences(
        xs, ys, dtype=dtype, extras=ex or None, components=components,
        stress_components=stress_components, features=tuple(features), mode=mode, ndi=ndi,
    )


# ---------------------------------------------------------------------------
# StressLSTM CSV format
# ---------------------------------------------------------------------------

def load_csv(
    path,
    inputs: Sequence[str] = ("total_strain_xx", "total_strain_yy", "total_strain_xy"),
    targets: Sequence[str] = ("stress_xx", "stress_yy", "stress_xy"),
    id_col: str = "simulation_load_id",
    time_col: Optional[str] = "timestep",
    drop_first: bool = True,
    components: Sequence[str] = ("xx", "yy", "xy"),
    features: Sequence[str] = ("strain",),
    dtype=torch.float32,
) -> SequenceDataset:
    """Read a StressLSTM-style CSV (one row per time step, sequences grouped by
    ``id_col``; the initial ``timestep == 0`` row is dropped as in the reference
    implementation). Requires pandas."""
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError("load_csv requires pandas: pip install pandas") from exc
    df = pd.read_csv(path)
    if drop_first and time_col is not None and time_col in df.columns:
        df = df[df[time_col] != 0]
    xs, ys = [], []
    for _, g in df.groupby(id_col, sort=True):
        if time_col is not None and time_col in g.columns:
            g = g.sort_values(time_col)
        xs.append(assemble_features(features, g[list(inputs)].to_numpy(dtype=float)))
        ys.append(g[list(targets)].to_numpy(dtype=float))
    return SequenceDataset.from_sequences(
        xs, ys, dtype=dtype, components=tuple(components), stress_components=tuple(components),
        features=tuple(features), mode="csv",
    )
