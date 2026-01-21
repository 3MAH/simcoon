"""
Stress transfer library examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example demonstrates the stress transfer helpers exposed in
:mod:`simcoon.simmit`.

The goal is to show (in a single place) how to use:

* :func:`simcoon.simmit.stress_convert` for stress-measure conversions
* the batch conventions expected by the Python bindings

**Batch conventions**

Throughout this script we use the shapes expected by :func:`stress_convert`:

* stress batches: ``(6, N)``
* deformation gradient batches: ``(3, 3, N)``

We deliberately use a *general* deformation gradient (stretch + shear) to avoid
special cases where some measures can accidentally coincide.

"""

###################################################################################
# Notes
#
# * The Voigt order in simcoon is assumed to be ``[11, 22, 33, 12, 23, 13]``.
# * The conversion keys are validated by the bindings (unknown keys raise).


###################################################################################
# Imports

from __future__ import annotations
import numpy as np
from simcoon import simmit as sim

###################################################################################
# Helper: pretty printing


def _fmt6(v: np.ndarray) -> str:
    v = np.asarray(v, dtype=float).reshape(6)
    return np.array2string(v, precision=6, suppress_small=True)


###################################################################################
# Prescribed deformation gradient(s)
#
# We use a small family of deformation gradients with both stretch and shear.
# This keeps the example generic.

N = 4


def _make_F(lam1: float, lam2: float, lam3: float, gamma12: float) -> np.ndarray:
    return np.array(
        [
            [lam1, gamma12, 0.0],
            [0.0, lam2, 0.0],
            [0.0, 0.0, lam3],
        ],
        dtype=float,
    )


lam1_vals = np.linspace(1.05, 1.25, N)
lam2_vals = np.linspace(0.95, 1.05, N)
lam3_vals = np.linspace(1.00, 1.10, N)
gamma12_vals = np.linspace(0.00, 0.15, N)

F_batch = np.zeros((3, 3, N), dtype=float)
for i in range(N):
    F_batch[:, :, i] = _make_F(
        float(lam1_vals[i]),
        float(lam2_vals[i]),
        float(lam3_vals[i]),
        float(gamma12_vals[i]),
    )

J_vals = np.linalg.det(np.transpose(F_batch, (2, 0, 1)))
print("J range:", float(J_vals.min()), "..", float(J_vals.max()))

###################################################################################
# Define a Cauchy stress batch
#
# Here we provide a synthetic Cauchy stress batch (6,N). This keeps the example
# focused on the *transfer* API rather than on a particular constitutive law.
#
# The only requirement is that it has the right shape.

sigma_cauchy_batch = np.zeros((6, N), dtype=float)

# Give it a simple, varying pattern across the batch.
sigma_cauchy_batch[0, :] = np.linspace(100.0, 150.0, N)  # s11
sigma_cauchy_batch[1, :] = np.linspace(10.0, 15.0, N)  # s22
sigma_cauchy_batch[2, :] = np.linspace(5.0, 7.0, N)  # s33
sigma_cauchy_batch[3, :] = np.linspace(20.0, 30.0, N)  # s12
sigma_cauchy_batch[4, :] = np.linspace(0.0, 2.0, N)  # s23
sigma_cauchy_batch[5, :] = np.linspace(1.0, 3.0, N)  # s13

print("\nInput Cauchy batch (6,N):")
print(sigma_cauchy_batch)

###################################################################################
# stress_convert: conversion keys
#
# The bindings validate keys. The following list reflects the keys exposed by the
# bindings (including Kirchhoff spelling aliases).

keys = [
    "Cauchy2PKI",
    "Cauchy2PKII",
    "Cauchy2Kirchhoff",
    "PKI2Cauchy",
    "PKII2Cauchy",
    "Kirchhoff2Cauchy",
    "Kirchhoff2PKI",
    "Kirchhoff2PKII",
    "PKI2Kirchhoff",
    "PKII2Kirchhoff",
]

###################################################################################
# Demonstrate batch conversions

converted: dict[str, np.ndarray] = {}
for key in keys:
    converted[key] = np.asarray(
        sim.stress_convert(sigma_cauchy_batch, F_batch, key), dtype=float
    )

print("\nBatch conversion results (showing column 0 only):")
for key in keys:
    print(f"- {key:15s} -> {_fmt6(converted[key][:, 0])}")

###################################################################################
# Round-trip sanity checks (informal)
#
# These are not strict asserts (we keep the script lightweight), but printing the
# max error helps confirm that key pairs are consistent.

sigma_back_from_pki = np.asarray(
    sim.stress_convert(converted["Cauchy2PKI"], F_batch, "PKI2Cauchy"), dtype=float
)
err_pki = np.max(np.abs(sigma_back_from_pki - sigma_cauchy_batch))
print("\nRound-trip error: Cauchy -> PKI -> Cauchy =", float(err_pki))

sigma_back_from_pk2 = np.asarray(
    sim.stress_convert(converted["Cauchy2PKII"], F_batch, "PKII2Cauchy"), dtype=float
)
err_pk2 = np.max(np.abs(sigma_back_from_pk2 - sigma_cauchy_batch))
print("Round-trip error: Cauchy -> PKII -> Cauchy =", float(err_pk2))

sigma_back_from_tau = np.asarray(
    sim.stress_convert(converted["Cauchy2Kirchhoff"], F_batch, "Kirchhoff2Cauchy"),
    dtype=float,
)
err_tau = np.max(np.abs(sigma_back_from_tau - sigma_cauchy_batch))
print("Round-trip error: Cauchy -> Kirchhoff -> Cauchy =", float(err_tau))

###################################################################################
# Tensor <-> Voigt helper demo
#
# Show how to convert one Voigt vector to a 3x3 tensor and back.

sigma0_v = np.ascontiguousarray(sigma_cauchy_batch[:, 0])
sigma0_t = sim.v2t_stress(sigma0_v)
sigma0_v2 = np.asarray(sim.t2v_stress(sigma0_t), dtype=float).reshape(6)

print("\nVoigt/tensor helper check (sample 0):")
print("sigma0_v :", _fmt6(sigma0_v))
print("sigma0_t :\n", np.asarray(sigma0_t, dtype=float))
print("back to v:", _fmt6(sigma0_v2))
