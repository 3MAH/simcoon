"""Block-diagonal tangent of a sequence model, as a diagnostic.

The quantity the material-point (or finite-element) Newton loop consumes is the
derivative of the stress of one increment with respect to the strain of that increment,
at frozen start state — the *total* derivative through the internal variables, which for
a classical model is assembled analytically as
:math:`\\mathbf{D}^\\varepsilon = \\mathbf{L} - \\sum_j \\boldsymbol{\\kappa}^j
\\mathbf{P}^j_\\varepsilon` (simcoon theory manual, numerical methods) and which the
solver records in ``SolverResults["TangentMatrix"]``.

For a recurrent surrogate the analogue is
:math:`\\partial \\boldsymbol{\\sigma}_t / \\partial \\boldsymbol{\\varepsilon}_t` with
:math:`\\mathbf{s}_{t-1}` held fixed, which is exactly what
:class:`simcoon.ml.RecurrentLaw` hands to the solver. :func:`sequence_tangent` evaluates it
for every step of a batch of sequences at once, so it can be compared with the reference
operator the solver recorded (``generate_dataset(..., record=("tangent",))``) or with a
cell's closed form (:meth:`~simcoon.ml.cells.StateModel.analytic_tangent`).
"""

from __future__ import annotations

from typing import Dict, Optional

import torch
from torch import Tensor

from .cells import StateModel

#: whether ``is_grads_batched`` works for this cell's backward, per device type (MPS has
#: no batching rule for the LSTM one); decided on first use, one entry per backend
_BATCHED_GRAD: Dict[str, bool] = {}


@torch.no_grad()
def _states_before(model: StateModel, x: Tensor) -> Tensor:
    """State ``s_{t-1}`` seen at every step, flattened over ``(N, T)``."""
    N, T, _ = x.shape
    s = model.zero_state(N, dtype=x.dtype, device=x.device)
    seen = []
    for t in range(T):
        seen.append(s)
        _, s = model.step(x[:, t], s)
    # (N*T, state_size), matching x.reshape(N*T, .)
    return torch.stack(seen, dim=1).reshape(N * T, model.state_size).contiguous()


def sequence_tangent(model: StateModel, x: Tensor, strain_slice: Optional[slice] = None) -> Tensor:
    """``d sigma_t / d eps_t`` at frozen state ``s_{t-1}``, for every step of a batch.

    Parameters
    ----------
    model : StateModel
    x : Tensor (N, T, n_in)
        Input features in physical units.
    strain_slice : slice, optional
        Restrict the derivative to these columns of ``x``. By default **every strain-like
        block moves together**: perturbing the strain at the end of an increment changes
        the ``"strain"`` feature and the ``"dstrain"`` feature by the same amount, and the
        derivative is the sum of both contributions. Restricting to one block answers a
        different question and is rarely what a solver consumes.

    Returns
    -------
    Tensor (N, T, n_out, n_comp)

    Notes
    -----
    The states are computed once without gradients and then treated as constants, which
    makes the ``N x T`` steps independent: one batched forward and ``n_out`` batched
    backward passes give every block, instead of one backward per step.
    """
    N, T, n_in = x.shape
    s = _states_before(model, x)
    # selector writing a strain perturbation into every strain-like feature block
    sel = torch.zeros(model.n_comp, n_in, dtype=x.dtype, device=x.device)
    if strain_slice is None:
        i = 0
        for f in model.features:
            w = model.n_comp if f in ("strain", "dstrain") else 1
            if f in ("strain", "dstrain"):
                sel[:, i:i + w] += torch.eye(model.n_comp, dtype=x.dtype, device=x.device)
            i += w
    else:
        w = len(range(*strain_slice.indices(n_in)))
        sel[:w, strain_slice] = torch.eye(w, dtype=x.dtype, device=x.device)[:w, :w]
    # a derivative is wanted even when the caller is inside torch.no_grad() (scoring a
    # validation batch, for instance): the graph built here is local to this function
    with torch.enable_grad():
        return _blocks(model, x, s, sel, N, T, n_in)


def _blocks(model, x, s, sel, N, T, n_in):
    flat = x.reshape(N * T, n_in)
    delta = torch.zeros(N * T, model.n_comp, dtype=x.dtype, device=x.device).requires_grad_(True)
    y, _ = model.step(flat.detach() + delta @ sel, s)         # (N*T, n_out)

    n_out = y.shape[1]
    dev = y.device.type
    if _BATCHED_GRAD.get(dev) is not False:
        try:
            eye = torch.zeros(n_out, N * T, n_out, dtype=y.dtype, device=y.device)
            eye[torch.arange(n_out), :, torch.arange(n_out)] = 1.0
            # keep the graph while undecided: a failed batched attempt frees the saved
            # tensors otherwise, and the per-row fallback below could not run
            (g,) = torch.autograd.grad(y, delta, grad_outputs=eye, is_grads_batched=True,
                                       create_graph=True, retain_graph=True)
            _BATCHED_GRAD[dev] = True
            return g.permute(1, 0, 2).reshape(N, T, n_out, -1)
        except RuntimeError:
            if _BATCHED_GRAD.get(dev):
                raise
            _BATCHED_GRAD[dev] = False   # no vmap rule here (MPS): one pass per row
    cols = [torch.autograd.grad(y[:, k].sum(), delta, create_graph=True, retain_graph=True)[0]
            for k in range(n_out)]
    return torch.stack(cols, dim=1).reshape(N, T, n_out, -1)
