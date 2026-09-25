"""
A stress LSTM trained on EPICP paths
====================================

Train the recurrent network of :mod:`simcoon.ml` (Danoun et al.) on random
non-proportional strain paths integrated with the elastoplastic ``EPICP`` model,
evaluate it with the identification metrics, then run it **as a constitutive
law** in the simcoon solver: under strain control (as in training) and under
stress control, where the Newton loop uses the autograd tangent of the network.

Requires PyTorch (conda-forge ``pytorch`` in a conda environment, ``pip install simcoon[ml]``
otherwise). The training budget below is kept
small for the documentation build; increase ``N_PATHS`` and ``EPOCHS`` for an
accurate surrogate (the reference setting is 7000 paths, 2000 epochs).
"""

import matplotlib.pyplot as plt
import torch

from simcoon import ml
from simcoon.solver import StepMeca, solve

torch.manual_seed(0)
N_PATHS, EPOCHS = 200, 150
EPICP_PROPS = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]   # E, nu, alpha, sigma_Y, k, m
EPICP_NSTATEV = 8

###############################################################################
# 1. Data: random strain paths integrated with EPICP
# ---------------------------------------------------
# Four linear segments of 25 increments towards uniformly drawn targets
# (+/- 2 % normal, +/- 4 % shear here), integrated by the material-point solver.

targets, ninc = ml.random_strain_paths(N_PATHS, n_segments=4, n_sub=25, seed=0,
                                       amplitude=[0.02, 0.02, 0.02, 0.04, 0.04, 0.04])
ds = ml.generate_dataset("EPICP", EPICP_PROPS, EPICP_NSTATEV, targets=targets, ninc=ninc, mode="3D")
train_ds, test_ds = ml.split_dataset(ds, test_size=0.3, seed=0)
print(f"sequences: {len(train_ds)} train / {len(test_ds)} test, {ds.x.shape[1]} steps")

###############################################################################
# 2. Training
# -----------
# StressLSTM architecture (2 LSTM layers, linear head, per-component
# standardisation), MSE on standardised data with Adam. The loss is the
# differentiable twin of :func:`simcoon.identify.calc_cost`.

model = ml.StressLSTM(hidden_size=64, num_layers=2)
train_losses, val_losses = ml.train(model, train_ds, test_ds, epochs=EPOCHS, batch_size=32,
                                    lr=2e-3, loss="mse", log_every=50)

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.semilogy(train_losses, label="train")
ax.semilogy(val_losses, label="test")
ax.set_xlabel("epoch")
ax.set_ylabel("MSE (standardised)")
ax.legend()
ax.grid(True, which="both", ls=":")
fig.tight_layout()

###############################################################################
# 3. Evaluation with the identification metrics
# ---------------------------------------------

report = ml.evaluate(model, test_ds, metrics=("nmse", "wmape"), per_component=True)
print("test NMSE = %.3e, wMAPE = %.2f %%" % (report["nmse"], 100 * report["wmape"]))
for name, m in report["per_component"].items():
    print(f"  {name}: NMSE = {m['nmse']:.3e}, wMAPE = {100 * m['wmape']:.2f} %")

###############################################################################
# 4. The LSTM as a constitutive law in the solver
# -----------------------------------------------
# :class:`simcoon.ml.LSTMLaw` serves the network under the ``PYEXT`` name: the
# hidden state lives in ``statev`` and the tangent is the autograd Jacobian of the
# network. First a strain-controlled random path of the test set ...

law = ml.LSTMLaw(model)
k = 0
x_test = test_ds.x[k].numpy()
steps = []
for seg in range(4):
    steps.append(StepMeca(control=["strain"] * 6, value=x_test[25 * seg + 24], ninc=25))
res_lstm = solve(steps, law)
res_ref = solve(steps, "EPICP", EPICP_PROPS, EPICP_NSTATEV)

fig, axes = plt.subplots(2, 3, figsize=(10, 5.5), sharex=True)
for i, (ax, name) in enumerate(zip(axes.ravel(), ml.VOIGT)):
    ax.plot(res_ref["Stress"][i], "k-", lw=2, label="EPICP")
    ax.plot(res_lstm["Stress"][i], "r--", lw=1.5, label="LSTM (solver)")
    ax.set_title(r"$\sigma_{%s}$" % name)
    ax.grid(True)
axes[0, 0].legend()
axes[1, 0].set_xlabel("increment")
fig.suptitle("Strain-controlled test path")
fig.tight_layout()

###############################################################################
# ... then a uniaxial tension **under stress control** (axial stress prescribed,
# lateral strains free): the Newton loop of the solver iterates on the strain
# with the tangent of the network.

tension = StepMeca(control=["stress"] * 6, value=[450.0, 0, 0, 0, 0, 0], ninc=50)
res_lstm = solve(tension, law)
res_ref = solve(tension, "EPICP", EPICP_PROPS, EPICP_NSTATEV)

fig, ax = plt.subplots(figsize=(5.5, 4))
ax.plot(res_ref["Strain"][0], res_ref["Stress"][0], "k-", lw=2, label="EPICP")
ax.plot(res_lstm["Strain"][0], res_lstm["Stress"][0], "r--", lw=1.5, label="LSTM (solver, stress control)")
ax.set_xlabel(r"$\varepsilon_{11}$")
ax.set_ylabel(r"$\sigma_{11}$ (MPa)")
ax.grid(True)
ax.legend()
fig.tight_layout()
plt.show()
