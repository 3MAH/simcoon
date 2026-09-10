"""Training loop for the recurrent constitutive cells (PyTorch).

The loss is :func:`~simcoon.ml.torch_cost`, the differentiable twin of the
identification cost :func:`simcoon.identify.calc_cost`, evaluated on standardised
targets (StressLSTM setting: MSE on standardised data, Adam). Padding is handled
through the point weights. The datasets are moved to the device once and batched by
index permutation (they are dense padded tensors, no collation needed).
"""

from __future__ import annotations

from typing import Callable, List, Optional, Sequence, Tuple

import torch

from .data import SequenceDataset
from .losses import torch_cost
from .cells import StateModel


def resolve_device(device: str = "auto") -> torch.device:
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device)


def _dissipation_penalty(model: StateModel, x: torch.Tensor, y_pred: torch.Tensor,
                         psi: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Discrete Clausius-Duhem penalty ``mean(relu(-(sigma . deps - dpsi))**2)`` on valid
    steps (isothermal; engineering shear strains make ``sigma . deps`` the work
    double contraction). Requires ``stress_components == components``."""
    if model.stress_components != model.components:
        raise ValueError("dissipation penalty requires stress_components == components")
    strain = x[..., :model.n_comp]
    prev = torch.cat([torch.zeros_like(strain[:, :1]), strain[:, :-1]], dim=1)
    deps = strain - prev
    prev_psi = torch.cat([torch.zeros_like(psi[:, :1]), psi[:, :-1]], dim=1)
    D = (y_pred * deps).sum(-1, keepdim=True) - (psi - prev_psi)
    scale = (y_pred.abs() * deps.abs()).sum(-1, keepdim=True).mean().clamp(min=1e-30)
    pen = torch.relu(-D / scale) ** 2
    m = mask.unsqueeze(-1).to(pen.dtype)
    return (pen * m).sum() / m.sum().clamp(min=1.0)


def train(
    model: StateModel,
    train_ds: SequenceDataset,
    val_ds: Optional[SequenceDataset] = None,
    epochs: int = 2000,
    batch_size: int = 64,
    lr: float = 1e-3,
    loss: str = "mse",
    w_response: Optional[Sequence[float]] = None,
    device: str = "auto",
    dissipation_weight: float = 0.0,
    fit_scalers: bool = True,
    optimizer: Optional[torch.optim.Optimizer] = None,
    scheduler=None,
    clip_grad_norm: Optional[float] = None,
    log_every: int = 50,
    verbose: bool = True,
    callback: Optional[Callable[[int, float, Optional[float]], None]] = None,
    seed: Optional[int] = None,
) -> Tuple[List[float], List[float]]:
    """Train a :class:`StressLSTM` on a :class:`SequenceDataset`.

    Parameters
    ----------
    model, train_ds, val_ds
        Network and datasets (``val_ds`` optional, evaluated every epoch).
    epochs, batch_size, lr
        Adam settings (StressLSTM: 2000 epochs, batch 64/128, lr 1e-3).
    loss : str
        Metric of :func:`torch_cost` used as loss (``"mse"``, ``"nmse"``,
        ``"nmse_per_response"``, ``"rmse"``, ``"mae"``, ``"wmape"``, ...), computed on
        standardised targets.
    w_response : sequence, optional
        Per-component weights (``w_response`` of ``calc_cost``).
    dissipation_weight : float
        Weight of the discrete Clausius-Duhem penalty (thermodynamic consistency,
        needs ``psi_head=True``; 0 = off).
    fit_scalers : bool
        Fit the model's standardisation buffers on ``train_ds`` before training.
    optimizer, scheduler
        Custom optimiser (default Adam) and optional LR scheduler stepped every epoch.
    clip_grad_norm : float, optional
        Cap on the gradient norm (``torch.nn.utils.clip_grad_norm_``). The reference
        recipe of the LMSC uses 1e-3.
    callback : callable(epoch, train_loss, val_loss), optional
    seed : int, optional
        ``torch.manual_seed`` for reproducible initialisation/shuffling.

    Returns
    -------
    train_losses, val_losses : list of float
        Per-epoch losses (``val_losses`` empty without ``val_ds``).
    """
    if seed is not None:
        torch.manual_seed(seed)
    dev = resolve_device(device)
    model.to(dev)
    train_ds = train_ds.to(dev)
    val_ds = None if val_ds is None else val_ds.to(dev)
    if fit_scalers:
        model.fit_scalers(train_ds.x, train_ds.y, train_ds.mask)
    wr = None if w_response is None else torch.as_tensor(w_response, dtype=train_ds.x.dtype)
    if dissipation_weight > 0.0 and getattr(model, "psi_head", None) is None:
        raise ValueError("dissipation_weight > 0 requires a model built with psi_head=True")
    optimizer = optimizer or torch.optim.Adam(model.parameters(), lr=lr)

    def batch_cost(ds: SequenceDataset, idx: torch.Tensor) -> torch.Tensor:
        x, y, mask = ds.x[idx], ds.y[idx], ds.mask[idx]
        weights = mask.to(x.dtype)
        y_pred, _, psi = model(x)
        cost = torch_cost(model.standardize_y(y), model.standardize_y(y_pred),
                          w_response=wr, w_point=weights, metric=loss)
        if dissipation_weight > 0.0:
            cost = cost + dissipation_weight * _dissipation_penalty(model, x, y_pred, psi, mask)
        return cost

    def epoch_loss(ds: SequenceDataset, training: bool) -> float:
        n = len(ds)
        order = torch.randperm(n, device=ds.x.device) if training else torch.arange(n, device=ds.x.device)
        total = 0.0
        for start in range(0, n, batch_size):
            idx = order[start:start + batch_size]
            if training:
                optimizer.zero_grad()
                cost = batch_cost(ds, idx)
                cost.backward()
                if clip_grad_norm is not None:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), clip_grad_norm)
                optimizer.step()
            else:
                with torch.no_grad():
                    cost = batch_cost(ds, idx)
            total += float(cost.detach()) * len(idx)
        return total / max(n, 1)

    train_losses: List[float] = []
    val_losses: List[float] = []
    for epoch in range(1, epochs + 1):
        model.train()
        tl = epoch_loss(train_ds, training=True)
        train_losses.append(tl)
        vl = None
        if val_ds is not None:
            model.eval()
            vl = epoch_loss(val_ds, training=False)
            val_losses.append(vl)
        if scheduler is not None:
            scheduler.step()
        if callback is not None:
            callback(epoch, tl, vl)
        if verbose and (epoch % log_every == 0 or epoch == 1 or epoch == epochs):
            msg = f"Epoch [{epoch}/{epochs}] train {loss}: {tl:.6e}"
            if vl is not None:
                msg += f"  val {loss}: {vl:.6e}"
            print(msg)
    model.eval()
    return train_losses, val_losses
