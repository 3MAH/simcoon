"""Stress LSTM: a gated recurrent network mapping a strain history to the stress.

Architecture of StressLSTM (Guevara Garban / Danoun et al.): a stacked
``torch.nn.LSTM`` followed by a linear projection head. Inputs and outputs are
standardised per component; the scalers are registered buffers of the module so
that a saved model carries them and :meth:`StressLSTM.forward` works directly on
physical quantities (strain in, stress out).

Known limitation, established by Bonatti and Mohr (*J. Mech. Phys. Solids* 158
(2022) 104697, Eqs. 9-11): a gated transition function
:math:`\\chi' = a(\\chi, \\Delta\\varepsilon) \\odot \\chi + b(\\chi, \\Delta\\varepsilon)`
has no mechanism forcing :math:`\\chi \\odot (1 - a(\\chi, 0)) = b(\\chi, 0)`, so it
almost surely loses information on zero-norm increments and its response depends
on the increment size. :class:`simcoon.ml.LMSC` fixes this by construction; here
the effect is only *bounded*, by the committed-state rule of the UMAT wrapper
(:class:`simcoon.ml.LSTMLaw`, ``commit_tol``).
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple, Sequence

import torch
from torch import Tensor, nn

from .cells import (
    VOIGT, StateModel, concat_features, n_inputs, register_cell, voigt_indices,
)

__all__ = ["VOIGT", "StressLSTM", "concat_features", "n_inputs", "voigt_indices"]


@register_cell
class StressLSTM(StateModel):
    """Stacked LSTM + linear head predicting stress components from a strain history.

    Parameters
    ----------
    components : sequence of str
        Active strain components fed to the network (subset of :data:`VOIGT`), e.g.
        ``("xx", "yy", "xy")`` for a 2D model, all six for 3D.
    stress_components : sequence of str, optional
        Predicted stress components (default: same as ``components``; a plane-strain
        model may add ``"zz"``).
    features : sequence of str
        Input features, in order: ``"strain"`` (always), optionally ``"dstrain"``
        (strain increment), ``"dtime"`` (time increment) and ``"temperature"``.
        The faithful StressLSTM setting is ``("strain",)``; add ``"dstrain"``
        when the model must be robust to a variable increment size.
    hidden_size, num_layers : int
        LSTM width and depth (64 and 2 in the reference implementation).
    psi_head : bool
        Add a scalar free-energy head (thermodynamic-consistency option of the
        ThC-RNN of Danoun et al.; off by default). When set, :meth:`forward`
        returns ``psi`` as a third value.

    Notes
    -----
    The flat state of :class:`~simcoon.ml.cells.StateModel` is the concatenation
    ``[h, c]``, each ``num_layers * hidden_size`` long.
    """

    #: The gated update is not self-consistent: the wrapper advances the state only
    #: past this fraction of the median training increment.
    default_commit_tol: float = 0.1

    def __init__(
        self,
        components: Sequence[str] = VOIGT,
        stress_components: Optional[Sequence[str]] = None,
        features: Sequence[str] = ("strain",),
        hidden_size: int = 64,
        num_layers: int = 2,
        psi_head: bool = False,
    ):
        super().__init__(components, stress_components, features)
        self.hidden_size = int(hidden_size)
        self.num_layers = int(num_layers)

        self.lstm = nn.LSTM(self.n_in, self.hidden_size, self.num_layers, batch_first=True)
        self.head = nn.Linear(self.hidden_size, self.n_out)
        self.psi_head = nn.Linear(self.hidden_size, 1) if psi_head else None

        self.register_buffer("x_mean", torch.zeros(self.n_in))
        self.register_buffer("x_std", torch.ones(self.n_in))

    # ------------------------------------------------------------------ meta
    @property
    def hparams(self) -> Dict:
        return dict(
            components=self.components, stress_components=self.stress_components,
            features=self.features, hidden_size=self.hidden_size,
            num_layers=self.num_layers, psi_head=self.psi_head is not None,
        )

    @property
    def state_size(self) -> int:
        """Scalars of the flat state ``[h, c]``: ``2 * num_layers * hidden_size``."""
        return 2 * self.num_layers * self.hidden_size

    # ------------------------------------------------------------ state pack
    def pack_state(self, h: Tensor, c: Tensor) -> Tensor:
        """``h``/``c (num_layers, B, hidden)`` -> flat state ``(B, state_size)``."""
        n = self.num_layers * self.hidden_size
        return torch.cat([h.permute(1, 0, 2).reshape(-1, n),
                          c.permute(1, 0, 2).reshape(-1, n)], dim=-1)

    def unpack_state(self, s: Tensor) -> Tuple[Tensor, Tensor]:
        """Flat state ``(B, state_size)`` -> ``h``, ``c (num_layers, B, hidden)``."""
        n = self.num_layers * self.hidden_size
        h = s[..., :n].reshape(-1, self.num_layers, self.hidden_size).permute(1, 0, 2)
        c = s[..., n:].reshape(-1, self.num_layers, self.hidden_size).permute(1, 0, 2)
        return h.contiguous(), c.contiguous()

    # -------------------------------------------------------------- scalers
    @torch.no_grad()
    def fit_scalers(self, x: Tensor, y: Tensor, mask: Optional[Tensor] = None) -> None:
        """Per-feature standardisation, on top of the response scalers of the base class."""
        super().fit_scalers(x, y, mask)
        m = None if mask is None else mask.bool()
        xf = x.reshape(-1, x.shape[-1]) if m is None else x[m]
        xm, xs = xf.mean(0), xf.std(0, unbiased=False)
        self.x_mean.copy_(xm.to(self.x_mean))
        self.x_std.copy_(torch.where(xs > 0, xs, torch.ones_like(xs)).to(self.x_std))

    def standardize_x(self, x: Tensor) -> Tensor:
        return (x - self.x_mean) / self.x_std

    # -------------------------------------------------------------- forward
    def forward(self, x: Tensor, s: Optional[Tensor] = None):
        """Run the network on a physical input sequence.

        Parameters
        ----------
        x : Tensor (B, T, n_in)
            Input features in physical units (see :meth:`build_inputs`).
        s : Tensor (B, state_size), optional
            Flat recurrent state; zeros when omitted.

        Returns
        -------
        y : Tensor (B, T, n_out)
            Stress in physical units.
        s : Tensor (B, state_size)
            Updated flat state.
        psi : Tensor (B, T, 1) or None
            Free energy when ``psi_head`` is enabled, ``None`` otherwise.
        """
        state = None if s is None else self.unpack_state(s)
        out, (h, c) = self.lstm(self.standardize_x(x), state)
        y = self.destandardize_y(self.head(out))
        psi = self.psi_head(out) if self.psi_head is not None else None
        return y, self.pack_state(h, c), psi

    def step(self, x: Tensor, s: Tensor) -> Tuple[Tensor, Tensor]:
        """One increment: ``x (B, n_in)``, ``s (B, state_size)`` -> ``y (B, n_out)``, ``s'``."""
        y, s1, _ = self.forward(x.unsqueeze(1), s)
        return y[:, 0, :], s1

    def step_psi(self, x: Tensor, s: Tensor):
        """:meth:`step` also returning ``psi (B, 1)`` (``None`` without ``psi_head``)."""
        y, s1, psi = self.forward(x.unsqueeze(1), s)
        return y[:, 0, :], s1, (psi[:, 0, :] if psi is not None else None)
