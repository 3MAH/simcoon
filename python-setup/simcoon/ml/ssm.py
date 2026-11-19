"""Arc-length state-space cell (ArcSSM): the LMSC update with input-selected coefficients.

In the vocabulary of state-space models the LMSC of Bonatti and Mohr
(:mod:`simcoon.ml.lmsc`) is a diagonal *selective* recurrence discretised in arc
length: ``chi' = chi + expm1(-nu alpha) (chi - beta)`` is the zero-order-hold
solution of ``d chi / d s = -alpha (chi - beta)`` over an arc-length increment
``nu``, the rate ``alpha`` and the target ``beta`` being selected by the input.
What makes the LMSC sequential at training time is that ``alpha`` and ``beta``
also depend on ``chi`` itself.

This cell keeps the update, hence its guarantees, and moves the state feedback
from *within* a layer to *across* layers, the way structured state-space models
(S4/S5, Mamba) and the constitutive state-space model of Barreira et al. [1]_ do::

    nu = ||d_eps||,  n = d_eps / nu
    layer 1:      z_1 = n                     -> alpha_1, beta_1
    layer l > 1:  z_l = [h'_{l-1} ; n]        -> alpha_l, beta_l
    h'_l  = h_l + expm1(-nu alpha_l) * (h_l - beta_l)
    sigma = W_sigma [h'_1 ; ... ; h'_L]        (linear, no bias)

with ``alpha_l = exp(.)``, ``beta_l = tanh(.)`` read off the quadratic layers of
the LMSC (Eq. 22 of [2]_). Within a layer the coefficients do not depend on that
layer's own state, so the recurrence is *linear* with time-varying coefficients
and a whole sequence is integrated by a parallel scan (:func:`linear_scan`)
instead of a Python loop with the coefficient network inside it. Kept by
construction, exactly as for the LMSC:

* **stationarity**: ``nu = 0`` leaves every layer's state bit for bit;
* **rate independence**: no time increment enters the equations;
* **self-consistency at frozen coefficients**, from the additivity of the
  exponential — and, for the first layer, *exact* self-consistency along any
  straight strain segment, since its coefficients depend on the direction only;
* ``h = 0`` is exactly the stress-free state (unbiased output map,
  :meth:`ArcSSM.fit_scalers` forces ``y_mean = 0``).

What is given up is the intra-layer feedback of the LMSC (the coefficients of a
layer reacting to its own state); the layers above see the layers below, which
is where the non-linearity in the state now lives. ``top_lmsc=True`` restores the
feedback on the last layer only: that layer becomes a genuine LMSC fed by the
layers below and is integrated step by step, the others stay parallel.

References
----------
.. [1] L. Barreira, A. Soydan, F. Scipione, M. A. Bessa, D. Mohr, *Constitutive
   state-space modeling of path-dependent plasticity: a resolution-consistent and
   parallelizable computational framework*, arXiv:2609.07294 (2026).
.. [2] C. Bonatti, D. Mohr, J. Mech. Phys. Solids 158 (2022) 104697.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import torch
import torch.nn.functional as F
from torch import Tensor, nn

from .cells import DIRECTION_MIN, VOIGT, StateModel, register_cell

__all__ = ["ArcSSM", "linear_scan"]


def linear_scan(a: Tensor, b: Tensor, h0: Tensor) -> Tensor:
    """All states of the linear recurrence ``h_t = a_t h_{t-1} + b_t`` at once.

    ``a, b (B, T, n)``, ``h0 (B, n)`` -> ``h (B, T, n)``. Hillis-Steele prefix
    composition of the affine maps ``h -> a_t h + b_t``: ``log2(T)`` rounds of
    element-wise products on the whole sequence, no Python loop over the steps.
    Stable whenever ``|a_t| <= 1`` (products can only shrink), which the
    arc-length update guarantees.
    """
    T = a.shape[1]
    # fold the initial state into the first step: h_0 = a_0 h0 + b_0
    b = torch.cat([b[:, :1] + a[:, :1] * h0.unsqueeze(1), b[:, 1:]], dim=1)
    k = 1
    while k < T:
        # compose each map with the composite k steps earlier; the identity (1, 0)
        # stands in for the steps before the start of the sequence
        a_prev = F.pad(a[:, :-k], (0, 0, k, 0), value=1.0)
        b_prev = F.pad(b[:, :-k], (0, 0, k, 0), value=0.0)
        b = a * b_prev + b
        a = a * a_prev
        k *= 2
    return b


class _Coefficients(nn.Module):
    """Quadratic layers of the LMSC (Eq. 22) mapping ``z`` to ``alpha > 0`` and ``beta``."""

    def __init__(self, n_in: int, n_state: int, depth: int, width: int, bias_alpha_init: float):
        super().__init__()
        self.qa, self.qb = nn.ModuleList(), nn.ModuleList()
        n = n_in
        for _ in range(depth):
            self.qa.append(nn.Linear(n, width))
            self.qb.append(nn.Linear(n, width))
            n = width
        self.to_alpha = nn.Linear(n, n_state)
        self.to_beta = nn.Linear(n, n_state)
        nn.init.constant_(self.to_alpha.bias, bias_alpha_init)

    def forward(self, z: Tensor) -> Tuple[Tensor, Tensor]:
        for a, b in zip(self.qa, self.qb):
            z = torch.tanh(a(z)) * torch.tanh(b(z))
        return torch.exp(self.to_alpha(z)), torch.tanh(self.to_beta(z))


@register_cell
class ArcSSM(StateModel):
    """Arc-length state-space cell: stacked LMSC layers with input-selected coefficients.

    Parameters
    ----------
    components : sequence of str
        Active strain components (subset of :data:`~simcoon.ml.VOIGT`).
    stress_components : sequence of str, optional
        Predicted stress components (default: same as ``components``).
    n_state : int
        State variables per layer; the flat state holds ``n_layers * n_state``.
    n_layers : int
        Number of layers. The first reads the increment direction only, each
        following one reads the updated state of the layer below and the direction.
    depth, width : int
        Number and width of the quadratic layers of every coefficient network
        (Eq. 22 of Bonatti and Mohr).
    bias_alpha_init : float
        Constant initialisation of the rate bias (``alpha`` starts at ``exp`` of it).
    top_lmsc : bool
        Feed the last layer its own state as well (an LMSC on top of the stack).
        That layer is then integrated step by step; the layers below stay parallel.

    Notes
    -----
    ``features`` is fixed to ``("strain", "dstrain")``: the cell needs the increment
    and ignores the strain block. The increment norm is taken on the engineering-shear
    Voigt components, like the LMSC.
    """

    #: Exact by construction: no committed-state rule needed.
    default_commit_tol: float = 0.0
    #: The tangent depends on the loading direction, undefined at a zero increment.
    directional_tangent: bool = True

    def __init__(
        self,
        components: Sequence[str] = VOIGT,
        stress_components: Optional[Sequence[str]] = None,
        n_state: int = 16,
        n_layers: int = 2,
        depth: int = 2,
        width: int = 32,
        bias_alpha_init: float = 3.0,
        top_lmsc: bool = False,
    ):
        super().__init__(components, stress_components, features=("strain", "dstrain"))
        if int(n_layers) < 1:
            raise ValueError("n_layers must be >= 1")
        self.n_state = int(n_state)
        self.n_layers = int(n_layers)
        self.depth = int(depth)
        self.width = int(width)
        self.bias_alpha_init = float(bias_alpha_init)
        self.top_lmsc = bool(top_lmsc)

        self.coef = nn.ModuleList()
        for layer in range(self.n_layers):
            n_in = self.n_comp + (self.n_state if layer > 0 else 0)
            if self.top_lmsc and layer == self.n_layers - 1:
                n_in += self.n_state                       # its own state, LMSC-style
            self.coef.append(_Coefficients(n_in, self.n_state, self.depth, self.width,
                                           self.bias_alpha_init))
        self.to_stress = nn.Linear(self.n_layers * self.n_state, self.n_out, bias=False)

    # ------------------------------------------------------------------ meta
    @property
    def hparams(self) -> Dict:
        return dict(
            components=self.components, stress_components=self.stress_components,
            n_state=self.n_state, n_layers=self.n_layers, depth=self.depth,
            width=self.width, bias_alpha_init=self.bias_alpha_init, top_lmsc=self.top_lmsc,
        )

    @property
    def state_size(self) -> int:
        return self.n_layers * self.n_state

    # -------------------------------------------------------------- scalers
    @torch.no_grad()
    def fit_scalers(self, x: Tensor, y: Tensor, mask: Optional[Tensor] = None) -> None:
        """Response scaling only, with ``y_mean`` forced to zero (unbiased output map)."""
        super().fit_scalers(x, y, mask)
        self.y_mean.zero_()

    # -------------------------------------------------------------- dynamics
    @staticmethod
    def _split(dstrain: Tensor) -> Tuple[Tensor, Tensor]:
        """Amplitude ``nu (..., 1)`` and direction ``n (..., n_comp)`` of the increment.

        The norm is computed on the increment scaled by its largest component so that it
        stays finite for any finite increment: a diverging Newton iterate of the solver
        (``|d_eps| ~ 1e154``) would otherwise overflow the sum of squares to ``inf`` and turn
        the autograd tangent into ``inf * 0 = nan``, where the closed form of the LMSC returns
        zero and lets the solver cut the step.
        """
        scale = dstrain.detach().abs().amax(dim=-1, keepdim=True).clamp_min(DIRECTION_MIN)
        nu = scale * (dstrain / scale).norm(dim=-1, keepdim=True)
        return nu, dstrain / nu.clamp_min(DIRECTION_MIN)

    def _is_top_lmsc(self, layer: int) -> bool:
        return self.top_lmsc and layer == self.n_layers - 1

    def _coefficients(self, layer: int, n: Tensor, below: Optional[Tensor],
                      own: Optional[Tensor] = None) -> Tuple[Tensor, Tensor]:
        """``alpha, beta`` of one layer from the direction, the layer below and (top LMSC) itself."""
        parts: List[Tensor] = [] if own is None else [own]
        if below is not None:
            parts.append(below)
        parts.append(n)
        return self.coef[layer](torch.cat(parts, dim=-1) if len(parts) > 1 else parts[0])

    def update(self, dstrain: Tensor, s: Tensor) -> Tensor:
        """One increment: ``dstrain (B, n_comp)``, ``s (B, state_size)`` -> ``s'``.

        Written in increment form, ``h' = h + expm1(-nu alpha) (h - beta)``, so that a
        zero increment returns ``h`` bit for bit (see :mod:`simcoon.ml.lmsc`).
        """
        nu, n = self._split(dstrain)
        h = s.reshape(s.shape[0], self.n_layers, self.n_state)
        out, below = [], None
        for layer in range(self.n_layers):
            own = h[:, layer]
            alpha, beta = self._coefficients(layer, n, below, own if self._is_top_lmsc(layer) else None)
            below = own + torch.expm1(-nu * alpha) * (own - beta)
            out.append(below)
        return torch.cat(out, dim=-1)

    def stress(self, s: Tensor) -> Tensor:
        """Physical stress of a state; ``s = 0`` gives exactly zero stress."""
        return self.destandardize_y(self.to_stress(s))

    def step(self, x: Tensor, s: Tensor) -> Tuple[Tensor, Tensor]:
        """One increment: ``x (B, n_in)``, ``s (B, state_size)`` -> ``sigma``, ``s'``."""
        s1 = self.update(self.split_features(x)["dstrain"], s)
        return self.stress(s1), s1

    def forward(self, x: Tensor, s: Optional[Tensor] = None):
        """Whole sequence ``x (B, T, n_in)`` -> ``y (B, T, n_out)``, ``s_T``, ``None``.

        Each layer's coefficients are evaluated for every step at once, then the
        linear recurrence is integrated by :func:`linear_scan`; only a top LMSC
        layer runs step by step.
        """
        nu, n = self._split(self.split_features(x)["dstrain"])            # (B,T,1), (B,T,n_comp)
        B, T = x.shape[0], x.shape[1]
        h0 = (self.zero_state(B, dtype=x.dtype, device=x.device) if s is None else s)
        h0 = h0.reshape(B, self.n_layers, self.n_state)
        out, below = [], None
        for layer in range(self.n_layers):
            if self._is_top_lmsc(layer):
                h, seq = h0[:, layer], []
                for t in range(T):
                    alpha, beta = self._coefficients(layer, n[:, t], None if below is None else below[:, t], h)
                    h = h + torch.expm1(-nu[:, t] * alpha) * (h - beta)
                    seq.append(h)
                below = torch.stack(seq, dim=1)
            else:
                alpha, beta = self._coefficients(layer, n, below)          # (B,T,n_state)
                g = torch.expm1(-nu * alpha)                               # in (-1, 0]
                # h_t = h_{t-1} + g_t (h_{t-1} - beta_t) = (1 + g_t) h_{t-1} - g_t beta_t
                below = linear_scan(1.0 + g, -g * beta, h0[:, layer])
            out.append(below)
        states = torch.cat(out, dim=-1)
        return self.stress(states), states[:, -1], None
