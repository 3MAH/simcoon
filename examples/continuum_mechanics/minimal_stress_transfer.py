"""
Minimal stress transfer example
===============================

Minimal hyperelastic example highlighting stress transfer helpers.

This script computes the Cauchy stress for a simple Neo-Hookean model under
uniaxial stretch and demonstrates:

- Converting the 3x3 Cauchy stress to a 6-component Voigt vector.
- Converting back from Voigt to a 3x3 matrix.
- Performing a batch conversion of multiple Voigt stresses to Cauchy using
  `stress_convert` (example shows PKII -> Cauchy conversion as an example key).

"""

import numpy as np
from simcoon import simmit as sim

# Simple Neo-Hookean parameters
mu = 0.5
kappa = 1000.0

# Single uniaxial stretch
lam1 = 1.5
F = np.diag([lam1, 1.0, 1.0])
J = np.linalg.det(F)

# Compute right Cauchy-Green and left Cauchy-Green as helpers
b = sim.L_Cauchy_Green(F)

# compute derivatives for Neo-Hookean: dW/dI1_bar = mu/2, dW/dI2_bar = 0, dU/dJ = kappa*ln(J)
dWdI_1_bar = 0.5 * mu
dWdI_2_bar = 0.0
dUdJ = kappa * np.log(J)

# Compute isochoric and volumetric parts of the Cauchy stress
sigma_iso = sim.sigma_iso_hyper_invariants(
    float(dWdI_1_bar), float(dWdI_2_bar), b, J, False
)
sigma_vol = sim.sigma_vol_hyper(dUdJ, b, J, False)

# Total Cauchy stress (3x3)
sigma = sigma_iso + sigma_vol
print("Cauchy stress (3x3):")
print(sigma)

# Convert to Voigt (6,) vector
sigma_voigt = sim.t2v_stress(sigma)
print("\nCauchy stress in Voigt (6,):")
print(sigma_voigt)

# Ensure voigt is 1D (6,) not a column (6,1) to avoid extra singleton axes
sigma_voigt = np.squeeze(sigma_voigt)

# Convert back to 3x3 (sanity check)
sigma_mat = sim.v2t_stress(sigma_voigt)
print("\nConverted back to 3x3:")
print(sigma_mat)

# Batch example: create N identical Voigt columns and a matching F stack
# Make N configurable and robust; default to 4 for more realistic batches
N = 4
# sigma_batch must be shaped (6, N)
sigma_batch = np.repeat(sigma_voigt[:, np.newaxis], N, axis=1)  # shape (6, N)
# F_batch must be shaped (3, 3, N) as required by stress_convert
F_batch = np.repeat(F[:, :, np.newaxis], N, axis=2)  # shape (3, 3, N)

# Sanity checks to match the C++ wrapper expectations
assert (
    sigma_batch.ndim == 2 and sigma_batch.shape[0] == 6 and sigma_batch.shape[1] == N
), f"sigma_batch must have shape (6,N), got {sigma_batch.shape}"
assert (
    F_batch.ndim == 3
    and F_batch.shape[0] == 3
    and F_batch.shape[1] == 3
    and F_batch.shape[2] == N
), f"F_batch must have shape (3,3,N), got {F_batch.shape}"

# Use stress_convert for a batch conversion (here we assume input is PKII and want Cauchy)
# In this minimal example the input is actually Cauchy; we just demonstrate the call shape.
converted = sim.stress_convert(sigma_batch, F_batch, "PKII2Cauchy")
print("\nBatch converted shape:", converted.shape)
print(converted)
