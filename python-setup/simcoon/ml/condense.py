"""Static condensation for reduced stress states (``ndi`` = 2 plane stress, ``ndi`` = 1 uniaxial).

A network trained in 3D knows all six components. Used under plane stress, the
out-of-plane strain is unknown and must be found so that the out-of-plane stress
vanishes — the same problem the elastic predictor of the simcoon kernels solves in
closed form (``el_pred(ndi=2)``: ``Q = L_FF - L_FS L_SS^-1 L_SF``). For a general
(non-linear, history-dependent) model the strain of the stress-free directions
``S`` is obtained by a local Newton iteration on ``sigma_S(eps_F, eps_S) = 0`` at
frozen internal state, using the autograd Jacobian; the tangent returned to the
caller is the condensed operator over the free directions ``F``.

This module also owns the ``ndi`` convention (which Voigt components are stress-free):
:func:`control_from_ndi` builds the solver control of :func:`simcoon.ml.generate_dataset`
from the same table :func:`condensed_newton` condenses on, so data generation and
condensation cannot drift apart.

The routine is vectorised over material points (``N``) so that the same code
serves the single-point UMAT and a batched finite-element evaluation.
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence, Tuple

import torch
from torch import Tensor

from .cells import VOIGT

#: Stress-free Voigt components by ``ndi`` (classical convention: 3 = 3D / plane strain,
#: 2 = plane stress, 1 = uniaxial).
ZERO_STRESS = {3: (), 2: ("zz",), 1: ("yy", "zz")}

def control_from_ndi(ndi: int) -> Tuple[str, ...]:
    """Solver control per Voigt component (``"strain"``/``"stress"``) generating data
    consistent with ``ndi``: the stress-free directions are stress-driven (to zero), all
    other components strain-driven."""
    zero = ZERO_STRESS[int(ndi)]
    return tuple("stress" if c in zero else "strain" for c in VOIGT)


class CondensationError(RuntimeError):
    """Raised when the local Newton iteration on the stress-free strains fails."""


def condensed_newton(
    fn: Callable[[Tensor, Optional[Sequence[int]]], Tuple[Tensor, Tensor]],
    eps: Tensor,
    zero_idx: Sequence[int],
    tol: float = 1e-8,
    maxiter: int = 25,
) -> Tuple[Tensor, Tensor, Tensor]:
    """Find ``eps[:, zero_idx]`` such that ``sigma[:, zero_idx] = 0`` (batched Newton).

    Parameters
    ----------
    fn : callable
        ``fn(eps (N, n), rows) -> (sigma (N, n), J (N, len(rows), n))`` with
        ``J = d sigma[rows] / d eps`` at frozen internal state (``rows=None`` = all rows).
        Only the stress-free rows are differentiated during the iterations; the full
        Jacobian is requested once, at the converged strain. Stress and strain
        components must be in the same order (``stress_components == components``).
    eps : Tensor (N, n)
        Strain; the entries at ``zero_idx`` are the starting guess.
    zero_idx : sequence of int
        Indices (in ``eps``/``sigma``) of the stress-free directions ``S``.
    tol : float
        Convergence: ``max|sigma_S| <= tol * max(1, max|sigma|)``.

    Returns
    -------
    eps : Tensor (N, n)
        Strain with the converged stress-free entries.
    sigma : Tensor (N, n)
        Stress (exactly zero at ``zero_idx``).
    Lt : Tensor (N, n, n)
        Condensed tangent: ``L_FF - L_FS L_SS^-1 L_SF`` on the free directions, zero rows
        and columns at ``zero_idx``.
    """
    if maxiter < 1:
        raise ValueError(f"maxiter must be >= 1, got {maxiter}")
    eps = eps.clone()
    N, n = eps.shape
    S = list(zero_idx)
    F = [i for i in range(n) if i not in S]
    for _ in range(maxiter):
        sigma, J_S = fn(eps, S)                      # J_S: (N, |S|, n)
        r = sigma[:, S]
        scale = torch.clamp(sigma.abs().amax(dim=1, keepdim=True), min=1.0)
        if bool((r.abs() <= tol * scale).all()):
            break
        try:
            deps = torch.linalg.solve(J_S[:, :, S], -r.unsqueeze(-1)).squeeze(-1)
        except RuntimeError as exc:
            raise CondensationError(f"singular condensation Jacobian: {exc}") from exc
        eps[:, S] = eps[:, S] + deps
    else:
        raise CondensationError(
            f"plane-stress/uniaxial condensation did not converge in {maxiter} iterations "
            f"(max residual {float(r.abs().max()):.3e})")
    sigma, J = fn(eps, None)                         # full Jacobian at the converged strain
    JF, JS = J[:, F], J[:, S]
    Q = JF[:, :, F] - JF[:, :, S] @ torch.linalg.solve(JS[:, :, S], JS[:, :, F])
    Lt = torch.zeros_like(J)
    Fi = torch.as_tensor(F, device=J.device)
    Lt[:, Fi.unsqueeze(1), Fi.unsqueeze(0)] = Q
    sigma = sigma.clone()
    sigma[:, S] = 0.0
    return eps, sigma, Lt
