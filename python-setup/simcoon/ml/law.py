"""A trained recurrent cell as a constitutive law (:class:`simcoon.PythonUMAT`).

:class:`RecurrentLaw` wraps any :class:`~simcoon.ml.cells.StateModel`
(:class:`~simcoon.ml.StressLSTM`, :class:`~simcoon.ml.LMSC`) so that the C++ solver
(``simcoon.solver.solve(blocks, law)``), the batch entry point
``sim.umat("PYEXT", ...)`` and finite-element codes (through
:meth:`RecurrentLaw.step_batch`) drive it like any other kernel:

* the cell's flat state is stored in ``statev`` — the solver's start-of-increment
  rollback therefore rewinds the network correctly on Newton retrials and step cuts;
* the tangent is the cell's closed form when it has one
  (:meth:`~simcoon.ml.cells.StateModel.analytic_tangent`, the LMSC case) and the
  autograd Jacobian of the step at frozen state otherwise (the gated case);
  ``tangent_mode = 0`` returns the elastic operator identified at the origin;
* ``ndi`` follows the classical convention (3 = 3D / plane strain, 2 = plane stress,
  1 = uniaxial): a model whose components do not include the stress-free directions
  is used as is (it was trained on data satisfying the constraint); a 3D model is
  condensed by the local Newton of :mod:`simcoon.ml.condense`.

Limits: the model is a small-strain law. Under finite strain the solver feeds the
corotational logarithmic strain (log_R route) like any small-strain kernel, but the
recurrent state cannot be rotated by ``DR``; the law is therefore meant for moderate
rotations.
"""

from __future__ import annotations

import copy
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from torch import Tensor

from ..pyumat import PythonUMAT, StepCut
from .cells import DIRECTION_MIN, StateModel, voigt_indices
from .condense import ZERO_STRESS, CondensationError, condensed_newton


class RecurrentLaw(PythonUMAT):
    """A recurrent constitutive cell served as a simcoon UMAT.

    Parameters
    ----------
    model : StateModel
        Trained cell. The law works on a private copy (moved to ``device``/``dtype``,
        evaluation mode, gradients disabled); the caller's model stays trainable.
    elastic_L : ndarray (6, 6), optional
        Elastic operator returned for ``tangent_mode = 0``; its diagonal also serves as
        the fallback stiffness of the components the model does not know (kept so that
        the solver's Jacobian stays regular when such a component is stress-driven).
        Default: tangent of the cell at zero strain and zero state. Zero diagonal
        entries are completed with the largest active stiffness in both cases.
    device, dtype
        Torch device and floating type used for inference (float64 by default so that
        the tangent matches the solver's precision).
    newton_tol, newton_maxiter
        Tolerance and iteration cap of the plane-stress / uniaxial condensation.
    commit_tol : float, optional
        Committed-state tolerance (strain norm over the active components). The recurrent
        state advances only when the strain moved by at least ``commit_tol`` since the last
        committed strain; a smaller move is answered with a trial stress computed from the
        committed state, which is left untouched (the trial-state contract of an implicit
        material routine). It bounds the increment-size sensitivity of a *gated* cell.
        Default: ``model.default_commit_tol`` times the median training increment recorded
        by ``fit_scalers`` — zero for a self-consistent cell such as the LMSC, which needs
        no such rule.
    max_increment : float
        Norm of the trial strain increment beyond which :meth:`integrate` raises a
        :class:`simcoon.StepCut` instead of evaluating the cell. A diverging Newton iterate
        of the solver (``|dE| ~ 1e154`` has been observed under stress control) is outside
        the domain of any small-strain law; answering it with a step cut lets the solver
        restart the increment smaller, where evaluating it overflows the energy or the
        tangent into ``inf``/``nan``. The default (a 100 % strain increment) never triggers
        on a sane trial.

    Notes
    -----
    ``statev`` layout: ``[s (state_size), strain_c (n_comp)]`` where ``strain_c`` is the
    committed strain: the active strain at which the recurrent state was last advanced
    (input of the ``"dstrain"`` feature, start of the condensation).
    """

    props = np.zeros(0)

    def __init__(
        self,
        model: StateModel,
        elastic_L: Optional[np.ndarray] = None,
        device: str = "cpu",
        dtype=torch.float64,
        newton_tol: float = 1e-8,
        newton_maxiter: int = 25,
        commit_tol: Optional[float] = None,
        max_increment: float = 1.0,
    ):
        self.device = torch.device(device)
        self.dtype = dtype
        self.model = copy.deepcopy(model).to(device=self.device, dtype=dtype).eval()
        for p in self.model.parameters():
            p.requires_grad_(False)
        self.components = tuple(model.components)
        self.stress_components = tuple(model.stress_components)
        self.features = tuple(model.features)
        self.idx_in = list(voigt_indices(self.components))
        self.idx_out = list(voigt_indices(self.stress_components))
        self.n_comp = len(self.components)
        self.state_size = model.state_size
        self._sl_s = slice(0, self.state_size)
        self._sl_prev = slice(self.state_size, self.state_size + self.n_comp)
        self.nstatev = self.state_size + self.n_comp
        self.newton_tol = float(newton_tol)
        self.newton_maxiter = int(newton_maxiter)
        if commit_tol is None:
            commit_tol = (model.default_commit_tol
                          * float(getattr(model, "median_increment", torch.zeros(()))))
        self.commit_tol = float(commit_tol)
        self.max_increment = float(max_increment)
        # components the model does not know (row or column): fallback diagonal stiffness
        self._fill_idx = [i for i in range(6) if i not in self.idx_out or i not in self.idx_in]
        # stress-free directions to condense, by ndi (positions in the model's components)
        self._zero_pos: Dict[int, List[int]] = {
            ndi: [self.components.index(c) for c in names if c in self.components]
            for ndi, names in ZERO_STRESS.items()
        }
        if any(self._zero_pos.values()) and self.stress_components != self.components:
            raise ValueError("plane-stress/uniaxial condensation requires stress_components == components")
        self._batched_grad: Optional[bool] = None   # decided on the first Jacobian
        L = self._identify_elastic_L() if elastic_L is None else np.array(elastic_L, dtype=float)
        self.elastic_L = self._complete_diagonal(L)
        # elastic operator per ndi (tangent_mode = 0): statically condensed on the stress-free
        # directions, the same algebra as el_pred(ndi) for a linear law
        self._elastic_by_ndi: Dict[int, np.ndarray] = {}
        for ndi, names in ZERO_STRESS.items():
            S = list(voigt_indices(names))
            Lc = self.elastic_L.copy()
            if S:
                F = [i for i in range(6) if i not in S]
                Q = Lc[np.ix_(F, F)] - Lc[np.ix_(F, S)] @ np.linalg.solve(Lc[np.ix_(S, S)], Lc[np.ix_(S, F)])
                Lc = np.zeros((6, 6))
                Lc[np.ix_(F, F)] = Q
            self._elastic_by_ndi[ndi] = Lc

    # ------------------------------------------------------------------ utils
    def _t(self, a) -> Tensor:
        a = np.asarray(a, dtype=float)
        if not a.flags.writeable:            # broadcast/read-only inputs (torch would warn)
            a = a.copy()
        return torch.as_tensor(a, dtype=self.dtype, device=self.device)

    def _scalar_field(self, v, n: int) -> Optional[Tensor]:
        if v is None:
            return None
        if np.isscalar(v):
            return torch.full((n,), float(v), dtype=self.dtype, device=self.device)
        return self._t(np.broadcast_to(v, (n,)))

    def zero_state(self, n: int) -> Tensor:
        return self.model.zero_state(n, dtype=self.dtype, device=self.device)

    @staticmethod
    def _complete_diagonal(L: np.ndarray) -> np.ndarray:
        diag = np.diag(L)
        fill = float(np.max(np.abs(diag))) if np.any(diag) else 1.0
        for i in range(6):
            if L[i, i] == 0.0:
                L[i, i] = fill
        return L

    def _jacobian(self, y: Tensor, strain: Tensor, rows: Optional[List[int]]) -> Tensor:
        """``d y[:, rows] / d strain`` as ``(N, len(rows), n_comp)``; one vectorised backward
        pass when torch supports it for this network, one pass per row otherwise."""
        N, n_out = y.shape
        rows = list(range(n_out)) if rows is None else list(rows)
        if self._batched_grad is not False:
            try:
                eye = torch.zeros(len(rows), N, n_out, dtype=y.dtype, device=y.device)
                eye[torch.arange(len(rows)), :, torch.as_tensor(rows, device=y.device)] = 1.0
                # keep the graph while undecided: a failed batched attempt frees the saved
                # tensors otherwise, and the per-row fallback below could not run
                (g,) = torch.autograd.grad(y, strain, grad_outputs=eye, is_grads_batched=True,
                                           retain_graph=(self._batched_grad is None))
                self._batched_grad = True
                return g.permute(1, 0, 2)
            except RuntimeError:
                if self._batched_grad:
                    raise
                self._batched_grad = False      # vmap over this backward unsupported (e.g. MPS): loop
        cols = []
        for k, r in enumerate(rows):
            (g,) = torch.autograd.grad(y[:, r].sum(), strain, retain_graph=(k < len(rows) - 1))
            cols.append(g)
        return torch.stack(cols, dim=1)

    def _step(self, strain: Tensor, s: Tensor, strain_prev: Tensor,
              dtime: Optional[Tensor], temperature: Optional[Tensor],
              rows: Optional[List[int]] = None, need_jac: bool = True):
        """Stress, Jacobian rows (or ``None``) and updated state, at frozen incoming state.

        The cell's closed-form tangent is used when it offers one; otherwise the step is
        replayed with a graph and differentiated by autograd.
        """
        def inputs(e):
            d = (e - strain_prev) if "dstrain" in self.features else None
            return self.model.build_inputs(e, d, dtime, temperature)

        e0 = strain.detach()
        if need_jac:
            x0 = inputs(e0)
            with torch.enable_grad():       # jacrev must not run under an outer no_grad
                J = self.model.analytic_tangent(x0, s)
            if J is not None:
                with torch.no_grad():
                    y, s1 = self.model.step(x0, s)
                J = J.detach()
                return y, (J if rows is None else J[:, list(rows), :]), s1

        e = e0.requires_grad_(need_jac)
        x = inputs(e)
        with torch.set_grad_enabled(need_jac):
            y, s1 = self.model.step(x, s)
        J = self._jacobian(y, e, rows).detach() if need_jac else None
        return y.detach(), J, s1.detach()

    def _identify_elastic_L(self) -> np.ndarray:
        """Tangent at zero strain and zero state, in a direction with no preferred axis."""
        s = self.zero_state(1)
        z = torch.zeros(1, self.n_comp, dtype=self.dtype, device=self.device)
        # a directional cell has no tangent at a strictly zero increment: probe it with a
        # tiny isotropic-in-index increment, which is the elastic branch at the virgin state
        prev = z if not getattr(self.model, "directional_tangent", False) else z - 1e-9
        _, J, _ = self._step(z, s, prev, self._scalar_field(1.0, 1), self._scalar_field(293.15, 1))
        L = np.zeros((6, 6))
        L[np.ix_(self.idx_out, self.idx_in)] = J[0].cpu().numpy()
        return L

    # -------------------------------------------------------------- batched
    def step_batch(
        self,
        strain: np.ndarray,
        state: np.ndarray,
        strain_prev: Optional[np.ndarray] = None,
        dtime=None,
        temperature=None,
        ndi: int = 3,
        tangent: bool = True,
    ) -> Tuple[np.ndarray, Optional[np.ndarray], np.ndarray, np.ndarray]:
        """One increment for ``N`` material points at once (fedoo-style layout).

        Parameters
        ----------
        strain : ndarray (6, N)
            Total strain at the end of the increment (engineering shear), all six
            components (inactive ones ignored).
        state : ndarray (state_size, N)
            Cell state at the beginning of the increment (zeros at start).
        strain_prev : ndarray (n_comp, N), optional
            Committed active strain (state last advanced there; ``"dstrain"`` feature /
            condensation start). Default: zeros.
        dtime, temperature : float or ndarray (N,), optional
        ndi : int
            3 (3D / plane strain), 2 (plane stress), 1 (uniaxial).
        tangent : bool
            Compute the tangent (``False`` for a cheaper residual-only evaluation; the
            condensation always needs it).

        Returns
        -------
        stress (6, N), Lt (6, 6, N) or None, state_new (state_size, N), strain_c (n_comp, N)
            New committed state: for points that moved less than ``commit_tol`` from
            ``strain_prev`` the returned state and committed strain are the inputs
            (trial evaluation, nothing committed).
        """
        strain = np.asarray(strain, dtype=float)
        N = strain.shape[1]
        eps = self._t(strain[self.idx_in].T)
        prev = torch.zeros_like(eps) if strain_prev is None else self._t(np.asarray(strain_prev).T)
        s0 = self._t(state).T.contiguous()
        dt = self._scalar_field(dtime, N)
        Tt = self._scalar_field(temperature, N)
        zero_pos = self._zero_pos.get(int(ndi), [])

        if zero_pos:
            last = {}

            def fn(e, rows):
                y, J, s1 = self._step(e, s0, prev, dt, Tt, rows=rows)
                last["state"] = s1
                return y, J

            try:
                eps, sig, J = condensed_newton(fn, eps, zero_pos, self.newton_tol, self.newton_maxiter)
            except CondensationError as exc:
                raise StepCut(msg=str(exc))
            s1 = last["state"]                # state of the last (converged) evaluation
        else:
            sig, J, s1 = self._step(eps, s0, prev, dt, Tt, need_jac=tangent)

        # committed-state rule: below the tolerance the evaluation is a trial, the state
        # and the committed strain are not advanced
        moved = (eps - prev).norm(dim=1)
        if self.commit_tol > 0.0:
            keep = moved < self.commit_tol                               # (N,)
            if bool(keep.any()):
                s1 = torch.where(keep.view(-1, 1), s0, s1)
                eps = torch.where(keep.view(-1, 1), prev, eps)

        stress = np.zeros((6, N))
        stress[self.idx_out] = sig.cpu().numpy().T
        Lt = None
        if J is not None:
            Lt = np.zeros((6, 6, N), order="F")
            Lt[np.ix_(self.idx_out, self.idx_in)] = J.cpu().numpy().transpose(1, 2, 0)
            if self._fill_idx:
                Lt[self._fill_idx, self._fill_idx, :] = np.diag(self.elastic_L)[self._fill_idx][:, None]
            # a directional tangent does not exist at a zero increment: elastic predictor,
            # the way a classical UMAT answers a zero trial increment
            if getattr(self.model, "directional_tangent", False):
                flat = (moved < DIRECTION_MIN).cpu().numpy()
                if flat.any():
                    Lt[:, :, flat] = self._elastic_by_ndi.get(int(ndi), self.elastic_L)[:, :, None]
        return stress, Lt, s1.cpu().numpy().T, eps.cpu().numpy().T

    # ------------------------------------------------------------- UMAT API
    def integrate(self, *, Etot, DEtot, sigma, statev, Wm, DTime, T, DT, ndi, start,
                  tangent_mode, **_):
        # no reset on `start`: the zero state IS the zero statev the solver provides, and a
        # coupler that keeps Time == 0 across increments must not lose the history
        statev = np.asarray(statev, dtype=float)
        DEtot = np.asarray(DEtot, dtype=float)
        # a diverging Newton iterate is outside any small-strain law's domain: ask for a
        # smaller increment before the energy or the tangent overflow
        norm_d = float(np.linalg.norm(DEtot))
        if not norm_d <= self.max_increment:
            raise StepCut(msg=f"trial strain increment |dE| = {norm_d:.2e} beyond "
                              f"max_increment = {self.max_increment:g}")
        eps6 = (np.asarray(Etot, dtype=float) + DEtot)[:, None]
        stress, Lt, s1, used = self.step_batch(
            eps6, statev[self._sl_s][:, None],
            strain_prev=statev[self._sl_prev][:, None], dtime=float(DTime),
            temperature=float(T + DT), ndi=int(ndi), tangent=(tangent_mode != 0))
        sig = stress[:, 0]
        L = self._elastic_by_ndi.get(int(ndi), self.elastic_L)
        if tangent_mode == 0:
            Lt = None                              # elastic operator (condensed for ndi < 3)
        Wm = np.array(Wm, dtype=float)
        Wm[0] += 0.5 * (np.asarray(sigma, dtype=float) + sig) @ DEtot
        new = np.empty(self.nstatev)
        new[self._sl_s] = s1[:, 0]
        new[self._sl_prev] = used[:, 0]
        return sig, (L if Lt is None else Lt[:, :, 0]), new, Wm, L

    # ------------------------------------------------------------- persist
    def save(self, path) -> None:
        self.model.save(path)

    @classmethod
    def load(cls, path, **kwargs) -> "RecurrentLaw":
        return cls(StateModel.load(path), **kwargs)


#: Backwards-compatible name: the law is not tied to the LSTM cell any more.
LSTMLaw = RecurrentLaw
