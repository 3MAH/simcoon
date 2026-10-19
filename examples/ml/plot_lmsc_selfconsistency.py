"""
Self-consistency: LMSC versus LSTM
==================================

A recurrent surrogate of a rate-independent law should answer the same thing whether a
strain path is described in 25 steps or 200. A gated cell has no mechanism enforcing that
(Bonatti and Mohr, *J. Mech. Phys. Solids* 158 (2022) 104697, Eqs. 9-11), while the
Linearized Minimal State Cell builds it into its update rule.

This example trains both cells on the same elasto-plastic paths, with the same budget,
then measures the two properties that separate them: **stationarity** (a held strain must
leave the stress unchanged) and **self-consistency** (refining the discretisation of a
path must not change the answer).
"""

import matplotlib.pyplot as plt
import numpy as np
import torch

from simcoon import ml
from simcoon.solver import StepMeca, solve

# %%
# Reference material and dataset: J2 plasticity with power-law isotropic hardening,
# random non-proportional strain paths integrated by the simcoon material-point solver.

E, NU = 70000.0, 0.3
PROPS = [E, NU, 1.0e-5, 300.0, 1000.0, 0.3]      # E nu alpha sigma_Y k m
NSTATEV = 8
FEATURES = ("strain", "dstrain")

torch.manual_seed(0)
targets, ninc = ml.random_strain_paths(240, n_segments=4, n_sub=25, amplitude=0.012, seed=1)
ds = ml.generate_dataset("EPICP", PROPS, NSTATEV, targets=targets, ninc=ninc, features=FEATURES)
train_ds, test_ds = ml.split_dataset(ds, test_size=0.25, seed=0)
print(f"{len(ds)} paths of {ds.x.shape[1]} steps, |sigma| max {ds.y.abs().max():.0f} MPa")

# %%
# Two cells, one budget. The LMSC carries 20 state variables against the LSTM's 256,
# and about fifteen times fewer parameters.

cells = {
    "LSTM": ml.StressLSTM(features=FEATURES, hidden_size=64, num_layers=2),
    "LMSC": ml.LMSC(n_state=20, depth=3, width=40),
}
laws = {}
for name, model in cells.items():
    ml.train(model, train_ds, epochs=250, batch_size=64,
             lr=2e-3 if name == "LSTM" else 5e-3, verbose=False, seed=0)
    rep = ml.evaluate(model, test_ds, metrics=("nmse", "wmape"))
    laws[name] = ml.LSTMLaw(model)
    print(f"{name}: {sum(p.numel() for p in model.parameters()):6d} parameters, "
          f"state {model.state_size:3d}, test NMSE {rep['nmse']:.2e}, "
          f"wMAPE {100 * rep['wmape']:.1f} %, commit_tol {laws[name].commit_tol:.2e}")

# %%
# Stationarity
# ------------
# Load a path, then hold the strain for a hundred further increments, the protocol of the
# reference. A stationary cell answers the same stress throughout the hold; a gated cell
# does not, because nothing in its training forces its transition function to be the
# identity at a zero increment.

# the gated cell is also shown without its committed-state rule, to separate what the
# architecture provides from what the wrapper patches
probes = dict(laws)
probes["LSTM (no commit rule)"] = ml.LSTMLaw(cells["LSTM"], commit_tol=0.0)

TARGET = [6e-3, -1e-3, 0.0, 2e-3, 0.0, 0.0]
hold = [StepMeca(control=["strain"] * 6, value=TARGET, ninc=50),
        StepMeca(control=["strain"] * 6, value=TARGET, ninc=100)]     # second block: held
print("\nstationarity: stress drift during 100 held-strain increments (MPa)")
ref = solve(hold, "EPICP", PROPS, NSTATEV)["Stress"][:, 50:]
print(f"  EPICP: {np.abs(ref - ref[:, [0]]).max():.3e}")
for name, law in probes.items():
    held = solve(hold, law, raise_on_abort=False)["Stress"][:, 50:]
    print(f"  {name}: {np.abs(held - held[:, [0]]).max():.3e}")

# %%
# Self-consistency
# ----------------
# The same monotonic tension, described with more and more increments. The reference law
# converges immediately; a surrogate should too.

ninc_list = [25, 50, 100, 200, 400]
curves = {name: [] for name in probes}
curves["reference"] = []
for n in ninc_list:
    step = StepMeca(control=["strain"] * 6,
                    value=[0.01, -0.003, -0.003, 0.004, 0.0, 0.0], ninc=n)
    curves["reference"].append(solve(step, "EPICP", PROPS, NSTATEV)["Stress"][0, -1])
    for name, law in probes.items():
        curves[name].append(solve(step, law, raise_on_abort=False)["Stress"][0, -1])

print("\nself-consistency: spread of the final stress over the refinement range (MPa)")
for name, v in curves.items():
    print(f"  {name}: {max(v) - min(v):8.2f}  ({', '.join(f'{x:.0f}' for x in v)})")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
for name, style in (("reference", "k-o"), ("LSTM", "C0--s"),
                    ("LSTM (no commit rule)", "C1:v"), ("LMSC", "C3-.^")):
    ax1.plot(ninc_list, curves[name], style, label=name)
ax1.set_xscale("log")
ax1.set_xlabel("increments describing the same path")
ax1.set_ylabel(r"final $\sigma_{11}$ (MPa)")
ax1.set_title("self-consistency")
ax1.grid(True, which="both", ls=":")
ax1.legend(fontsize=8)

# %%
# And the response itself, under uniaxial tension driven in mixed control.

st = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.01, 0, 0, 0, 0, 0], ninc=100)
ref = solve(st, "EPICP", PROPS, NSTATEV)
ax2.plot(100 * ref["Strain"][0], ref["Stress"][0], "k-", lw=2, label="EPICP")
for name, style in (("LSTM", "C0--"), ("LMSC", "C3-.")):
    r = solve(st, laws[name], raise_on_abort=False)
    ax2.plot(100 * r["Strain"][0], r["Stress"][0], style, lw=1.5, label=name)
ax2.set_xlabel(r"$\varepsilon_{11}$ (%)")
ax2.set_ylabel(r"$\sigma_{11}$ (MPa)")
ax2.set_title("uniaxial tension, mixed control")
ax2.grid(True)
ax2.legend()
fig.tight_layout()
plt.show()

# %%
# What the numbers say, at this deliberately small training budget:
#
# * **stationarity** is delivered either by the architecture (LMSC, exactly zero) or by
#   the wrapper's committed-state rule (LSTM, also zero once the rule is on); the gated
#   cell left to itself drifts by tens of MPa;
# * **self-consistency** separates them. Over a sixteenfold refinement the LMSC spreads
#   about as little as the reference law itself, while the gated cell moves by tens of MPa
#   even with the committed-state rule, which bounds the departure without removing it.
#
# The price is accuracy per epoch: the LMSC is a much smaller model with a constrained
# update, and it needs a larger training budget to reach the stress accuracy a gated cell
# gets quickly. It also smooths the elastic-plastic knee, because the linearization behind
# its update neglects second derivatives that are unbounded at the onset of plastic flow
# (Section 7.2 of the reference).
