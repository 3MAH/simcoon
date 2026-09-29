"""Linearized Minimal State Cell (LMSC) of Bonatti and Mohr.

Reference
---------
C. Bonatti, D. Mohr, *On the importance of self-consistency in recurrent neural
network models representing elasto-plastic solids*, J. Mech. Phys. Solids 158
(2022) 104697. The equations below are Eqs. (20)-(26) of Section 2.2. The MSC
proper is the earlier architecture (Bonatti and Mohr, Sci. Adv. 7, eabf3658,
2021); the LMSC is its self-consistent linearization.

Transition function, for a state :math:`\\chi` and a strain increment
:math:`\\Delta\\varepsilon`::

    nu    = ||d_eps||                                            (20)
    l_0   = [ chi ; d_eps / nu ]                                 (21)
    l_i   = tanh(Wa_i l_{i-1} + ba_i) * tanh(Wb_i l_{i-1} + bb_i) (22)
    alpha = exp(W_A l_d + b_A)                                   (23)
    beta  = tanh(W_B l_d + b_B)                                  (24)
    chi'  = exp(-nu alpha) * (chi - beta) + beta                  (25)
    sigma = W_s chi'                                             (26)

Three consequences that separate this cell from a gated one:

* **stationarity is exact**: ``nu = 0`` gives ``exp(0) = 1``, hence ``chi' = chi``
  and an unchanged stress, to machine precision. This is the condition
  (Eq. 11 of the paper) that LSTM and GRU cells cannot satisfy;
* **self-consistency is exact at frozen coefficients**, since
  ``exp(-nu1 a) exp(-nu2 a) = exp(-(nu1+nu2) a)``: splitting an increment changes
  nothing as long as ``alpha`` and ``beta`` do not vary along it;
* **rate independence is structural**: no time increment enters the equations.

The total strain is *not* an input: only the state and the increment direction
are. The output map is linear and unbiased so that ``chi = 0`` is exactly the
stress-free state, which is why :meth:`LMSC.fit_scalers` forces ``y_mean = 0``.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import torch
from torch import Tensor, nn

from .cells import DIRECTION_MIN, VOIGT, StateModel, register_cell

__all__ = ["LMSC"]


def _phi(t: Tensor) -> Tensor:
    """``(1 - exp(-t)) / t``, regular at ``t = 0`` (``phi(0) = 1``)."""
    small = t.abs() < 1e-6
    safe = torch.where(small, torch.ones_like(t), t)
    return torch.where(small, 1.0 - 0.5 * t + t * t / 6.0, -torch.expm1(-safe) / safe)


@register_cell
class LMSC(StateModel):
    """Linearized minimal state cell (Bonatti and Mohr, JMPS 2022).

    Parameters
    ----------
    components : sequence of str
        Active strain components (subset of :data:`~simcoon.ml.VOIGT`).
    stress_components : sequence of str, optional
        Predicted stress components (default: same as ``components``).
    n_state : int
        Size of the state vector. Six is the theoretical minimum for 3D von Mises
        plasticity and the size of the model the authors deploy; they report that
        the minimum is only reachable when training on long sequences of small
        increments, and that short sequences of large increments need excess state
        variables.
    depth, width : int
        Number ``d`` and width ``w`` of the quadratic layers of Eq. (22)
        (``d = 3``, ``w = 25`` in the deployed model). ``depth = 0`` maps the input
        straight to ``alpha`` and ``beta``.
    bias_alpha_init : float
        Constant initialisation of ``b_A`` (Eq. 23). The authors use 3, which
        starts the rates at ``exp(3)`` and keeps the initial dynamics slow.

    Notes
    -----
    ``features`` is fixed to ``("strain", "dstrain")``: the cell needs the
    increment and ignores the strain block, which the dataset carries anyway.
    The increment norm is taken on the engineering-shear Voigt components, the
    same convention as the rest of simcoon; training and inference share it.
    """

    #: Exact by construction: no committed-state rule needed.
    default_commit_tol: float = 0.0
    #: The tangent depends on the loading direction, undefined at a zero increment.
    directional_tangent: bool = True

    def __init__(
        self,
        components: Sequence[str] = VOIGT,
        stress_components: Optional[Sequence[str]] = None,
        n_state: int = 20,
        depth: int = 3,
        width: int = 40,
        bias_alpha_init: float = 3.0,
    ):
        super().__init__(components, stress_components, features=("strain", "dstrain"))
        self.n_state = int(n_state)
        self.depth = int(depth)
        self.width = int(width)
        self.bias_alpha_init = float(bias_alpha_init)

        # Eq. (22): d quadratic layers, each the product of two tanh branches
        self.qa, self.qb = nn.ModuleList(), nn.ModuleList()
        n = self.n_state + self.n_comp
        for _ in range(self.depth):
            self.qa.append(nn.Linear(n, self.width))
            self.qb.append(nn.Linear(n, self.width))
            n = self.width
        self.to_alpha = nn.Linear(n, self.n_state)          # Eq. (23), before exp
        self.to_beta = nn.Linear(n, self.n_state)           # Eq. (24), before tanh
        self.to_stress = nn.Linear(self.n_state, self.n_out, bias=False)   # Eq. (26)
        nn.init.constant_(self.to_alpha.bias, self.bias_alpha_init)

    # ------------------------------------------------------------------ meta
    @property
    def hparams(self) -> Dict:
        return dict(
            components=self.components, stress_components=self.stress_components,
            n_state=self.n_state, depth=self.depth, width=self.width,
            bias_alpha_init=self.bias_alpha_init,
        )

    @property
    def state_size(self) -> int:
        return self.n_state

    # -------------------------------------------------------------- scalers
    @torch.no_grad()
    def fit_scalers(self, x: Tensor, y: Tensor, mask: Optional[Tensor] = None) -> None:
        """Response scaling only, with ``y_mean`` forced to zero.

        The output map (Eq. 26) has no bias so that ``chi = 0`` is the stress-free
        state; a non-zero mean would destroy that property.
        """
        super().fit_scalers(x, y, mask)
        self.y_mean.zero_()

    # -------------------------------------------------------------- dynamics
    def _alpha_beta(self, chi: Tensor, n: Tensor) -> Tuple[Tensor, Tensor]:
        """Eqs. (21)-(24): rates ``alpha > 0`` and targets ``beta``, both ``(B, n_state)``."""
        z = torch.cat([chi, n], dim=-1)
        for a, b in zip(self.qa, self.qb):
            z = torch.tanh(a(z)) * torch.tanh(b(z))
        return torch.exp(self.to_alpha(z)), torch.tanh(self.to_beta(z))

    @staticmethod
    def _split(dstrain: Tensor) -> Tuple[Tensor, Tensor]:
        """Amplitude ``nu (B, 1)`` and direction ``n (B, n_comp)`` of the increment."""
        nu = dstrain.norm(dim=-1, keepdim=True)
        return nu, dstrain / nu.clamp_min(DIRECTION_MIN)

    def update(self, dstrain: Tensor, chi: Tensor) -> Tensor:
        """Eq. (25): new state from ``dstrain (B, n_comp)`` and ``chi (B, n_state)``.

        Written in increment form, ``chi' = chi - (1 - exp(-nu alpha)) (chi - beta)``,
        which is algebraically Eq. (25) but returns ``chi`` bit for bit at ``nu = 0``
        (stationarity exact, not merely to rounding) and avoids the cancellation of
        ``e (chi - beta) + beta`` for small increments.
        """
        nu, n = self._split(dstrain)
        alpha, beta = self._alpha_beta(chi, n)
        return chi + torch.expm1(-nu * alpha) * (chi - beta)

    def stress(self, chi: Tensor) -> Tensor:
        """Eq. (26) in physical units; ``chi = 0`` gives exactly zero stress."""
        return self.destandardize_y(self.to_stress(chi))

    def step(self, x: Tensor, s: Tensor) -> Tuple[Tensor, Tensor]:
        """One increment: ``x (B, n_in)``, ``s (B, n_state)`` -> ``sigma``, ``s'``."""
        chi = self.update(self.split_features(x)["dstrain"], s)
        return self.stress(chi), chi

    def forward(self, x: Tensor, s: Optional[Tensor] = None):
        """Run a whole sequence ``x (B, T, n_in)`` -> ``y (B, T, n_out)``, ``s (B, n_state)``.

        The recurrence is sequential by construction (the state update is not a
        cuDNN kernel); T steps of batched tensor algebra.
        """
        d = self.split_features(x)["dstrain"]
        chi = self.zero_state(x.shape[0], dtype=x.dtype, device=x.device) if s is None else s
        out = []
        for t in range(x.shape[1]):
            chi = self.update(d[:, t, :], chi)
            out.append(chi)
        return self.stress(torch.stack(out, dim=1)), chi, None

    # -------------------------------------------------------------- tangent
    def _dab_dn(self, chi: Tensor, n: Tensor) -> Tuple[Tensor, Tensor]:
        """Jacobians ``d alpha / d n`` and ``d beta / d n``, both ``(B, n_state, n_comp)``."""
        def f(chi_i, n_i):
            a, b = self._alpha_beta(chi_i.unsqueeze(0), n_i.unsqueeze(0))
            return a.squeeze(0), b.squeeze(0)

        try:
            return torch.func.vmap(torch.func.jacrev(f, argnums=1))(chi, n)
        except (RuntimeError, NotImplementedError):   # backend without a batching rule
            da, db = zip(*(torch.func.jacrev(f, argnums=1)(chi[i], n[i])
                           for i in range(chi.shape[0])))
            return torch.stack(da), torch.stack(db)

    def analytic_tangent(self, x: Tensor, s: Tensor) -> Optional[Tensor]:
        """Closed-form ``d sigma / d eps`` of :meth:`step`, ``(B, n_out, n_comp)``.

        Differentiating Eq. (25) with ``e = exp(-nu alpha)``, ``u = chi - beta`` and
        ``phi(t) = (1 - exp(-t)) / t``::

            d chi'/d eps = -diag(e*u) [ alpha n^T + (d alpha/d n)(I - n n^T) ]
                           + diag(alpha * phi(nu alpha)) (d beta/d n)(I - n n^T)
            d sigma/d eps = diag(y_std) W_s  d chi'/d eps

        The ``phi`` form keeps the expression regular as ``nu -> 0``; the limit stays
        direction-dependent, which is the correct behaviour of a rate-independent law
        (elastic stiffness on unloading, plastic tangent on loading). At a strictly zero
        increment the direction does not exist and the value returned here is
        meaningless; :attr:`directional_tangent` tells the wrapper to replace those
        points by the elastic predictor.
        """
        dstrain = self.split_features(x)["dstrain"]
        nu, n = self._split(dstrain)
        alpha, beta = self._alpha_beta(s, n)
        e = torch.exp(-nu * alpha)
        u = s - beta
        da_dn, db_dn = self._dab_dn(s, n)
        # projector onto the sphere, (B, n_comp, n_comp)
        proj = torch.eye(self.n_comp, dtype=x.dtype, device=x.device) - n.unsqueeze(-1) * n.unsqueeze(-2)
        dchi = (-(e * u).unsqueeze(-1) * (alpha.unsqueeze(-1) * n.unsqueeze(-2) + da_dn @ proj)
                + (alpha * _phi(nu * alpha)).unsqueeze(-1) * (db_dn @ proj))
        return self.y_std.view(1, -1, 1) * (self.to_stress.weight @ dchi)
