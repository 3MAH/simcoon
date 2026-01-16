"""
Yield Criteria Examples
~~~~~~~~~~~~~~~~~~~~~~~~~
This example demonstrates different yield criteria available in Simcoon.
"""

import numpy as np
from simcoon import simmit as sim
import matplotlib.pyplot as plt

###################################################################################
# Yield criteria are used to determine when a material begins to yield under
# a general stress state. Simcoon provides several yield criteria including
# von Mises, Tresca, Drucker, and anisotropic Hill criteria.
#
# The general form of a yield criterion is:
#
# .. math::
#
#    \Phi(\boldsymbol{\sigma}) = \bar{\sigma} - \sigma_Y \leq 0
#
# where :math:`\bar{\sigma}` is the equivalent stress and :math:`\sigma_Y` is the
# yield stress.

###################################################################################
# Von Mises equivalent stress
# ---------------------------
# The von Mises equivalent stress is defined as:
#
# .. math::
#
#    \bar{\sigma}_{vM} = \sqrt{\frac{3}{2} \mathbf{s} : \mathbf{s}}
#
# where :math:`\mathbf{s}` is the deviatoric stress tensor.

# Define a uniaxial stress state in Voigt notation [s11, s22, s33, s12, s13, s23]
sigma_uniaxial = np.array([400.0, 0.0, 0.0, 0.0, 0.0, 0.0])

sigma_Mises = sim.Mises_stress(sigma_uniaxial)
print(f"Uniaxial stress: σ11 = {sigma_uniaxial[0]} MPa")
print(f"von Mises equivalent stress: {sigma_Mises:.2f} MPa")

# The von Mises stress equals the uniaxial stress for this case
print(
    f"Verification: σ_vM = σ11 for uniaxial tension ✓"
    if abs(sigma_Mises - 400) < 1e-6
    else ""
)

###################################################################################
# Tresca equivalent stress
# ------------------------
# The Tresca criterion is based on the maximum shear stress:
#
# .. math::
#
#    \bar{\sigma}_{T} = \sigma_1 - \sigma_3
#
# where :math:`\sigma_1` and :math:`\sigma_3` are the maximum and minimum
# principal stresses.

sigma_Tresca = sim.Tresca_stress(sigma_uniaxial)
print(f"\nTresca equivalent stress: {sigma_Tresca:.2f} MPa")
print("For uniaxial tension, Tresca = von Mises = σ11")

# Compute the gradient (flow direction)
dsigma_Tresca = sim.dTresca_stress(sigma_uniaxial)
print(f"Tresca flow direction: {np.array_str(dsigma_Tresca, precision=4)}")

###################################################################################
# Drucker equivalent stress
# -------------------------
# The Drucker criterion generalizes the von Mises criterion:
#
# .. math::
#
#    \bar{\sigma}_{D} = \left[ (J_2)^3 - b (J_3)^2 \right]^{1/n}
#
# where :math:`J_2` and :math:`J_3` are the second and third invariants of the
# deviatoric stress, :math:`b` is a material parameter, and :math:`n` is an exponent.

b = 0.0  # When b=0 and n=2, Drucker reduces to von Mises
n = 2.0
props_Drucker = np.array([b, n])

sigma_Drucker = sim.Drucker_stress(sigma_uniaxial, props_Drucker)
print(f"\nDrucker equivalent stress (b={b}, n={n}): {sigma_Drucker:.2f} MPa")
print("With b=0 and n=2, Drucker ≈ von Mises")

# Compute the gradient
dsigma_Drucker = sim.dDrucker_stress(sigma_uniaxial, props_Drucker)
print(f"Drucker flow direction: {np.array_str(dsigma_Drucker, precision=4)}")

###################################################################################
# Hill anisotropic yield criterion
# --------------------------------
# The Hill criterion extends von Mises to orthotropic materials:
#
# .. math::
#
#    \bar{\sigma}_{H}^2 = F(\sigma_{22}-\sigma_{33})^2 + G(\sigma_{33}-\sigma_{11})^2 + H(\sigma_{11}-\sigma_{22})^2 + 2L\sigma_{23}^2 + 2M\sigma_{13}^2 + 2N\sigma_{12}^2
#
# For isotropic materials: F = G = H = 0.5 and L = M = N = 1.5

# Isotropic Hill parameters (equivalent to von Mises)
F, G, H = 0.5, 0.5, 0.5
L, M, N = 1.5, 1.5, 1.5
props_Hill = np.array([F, G, H, L, M, N])

sigma_Hill = sim.Hill_stress(sigma_uniaxial, props_Hill)
print(f"\nHill equivalent stress (isotropic params): {sigma_Hill:.2f} MPa")
print("With isotropic parameters, Hill = von Mises")

# Anisotropic case: material is stronger in direction 1
F_aniso, G_aniso, H_aniso = 0.3, 0.6, 0.6  # Lower F means σ22-σ33 matters less
L_aniso, M_aniso, N_aniso = 1.5, 1.5, 1.5
props_Hill_aniso = np.array([F_aniso, G_aniso, H_aniso, L_aniso, M_aniso, N_aniso])

sigma_Hill_aniso = sim.Hill_stress(sigma_uniaxial, props_Hill_aniso)
print(f"Hill equivalent stress (anisotropic): {sigma_Hill_aniso:.2f} MPa")

###################################################################################
# Hill configurational tensor
# ---------------------------
# The Hill criterion can be written in matrix form:
#
# .. math::
#
#    \bar{\sigma}_{H}^2 = \boldsymbol{\sigma}^T \mathbf{P} \boldsymbol{\sigma}
#
# where :math:`\mathbf{P}` is the Hill configurational tensor.

P_Hill = sim.P_Hill(props_Hill)
print("\nHill configurational tensor P (isotropic):")
print(np.array_str(P_Hill, precision=4, suppress_small=True))

###################################################################################
# Yield surface visualization
# ---------------------------
# Let's visualize the yield surfaces in the σ11-σ22 plane (plane stress).

theta = np.linspace(0, 2 * np.pi, 360)
sigma_Y = 400.0  # Yield stress

# von Mises yield surface (ellipse in principal stress space)
s11_vM = sigma_Y * np.cos(theta)
s22_vM = sigma_Y * (np.cos(theta) * 0.5 + np.sin(theta) * np.sqrt(3) / 2)

# Tresca yield surface (hexagon)
# Computing Tresca surface points
s11_T = []
s22_T = []
for t in theta:
    # Sample stress states on a circle and find where Tresca = sigma_Y
    s11 = sigma_Y * np.cos(t)
    s22 = sigma_Y * np.sin(t)
    sigma_test = np.array([s11, s22, 0.0, 0.0, 0.0, 0.0])
    sig_T = sim.Tresca_stress(sigma_test)
    if sig_T > 1e-6:
        scale = sigma_Y / sig_T
        s11_T.append(s11 * scale)
        s22_T.append(s22 * scale)

fig, ax = plt.subplots(figsize=(10, 10))

# Plot von Mises (approximated as ellipse)
n_points = 100
s11_plot = []
s22_plot = []
for i in range(n_points):
    t = 2 * np.pi * i / n_points
    s11 = 1.2 * sigma_Y * np.cos(t)
    s22 = 1.2 * sigma_Y * np.sin(t)
    sigma_test = np.array([s11, s22, 0.0, 0.0, 0.0, 0.0])
    sig_vM = sim.Mises_stress(sigma_test)
    if sig_vM > 1e-6:
        scale = sigma_Y / sig_vM
        s11_plot.append(s11 * scale)
        s22_plot.append(s22 * scale)

ax.plot(s11_plot, s22_plot, "b-", linewidth=2, label="von Mises")
ax.plot(s11_T, s22_T, "r--", linewidth=2, label="Tresca")

ax.set_xlabel(r"$\sigma_{11}$ (MPa)", fontsize=12)
ax.set_ylabel(r"$\sigma_{22}$ (MPa)", fontsize=12)
ax.set_title("Yield Surfaces in Plane Stress", fontsize=14)
ax.legend(fontsize=11)
ax.set_aspect("equal")
ax.grid(True, alpha=0.3)
ax.axhline(y=0, color="k", linewidth=0.5)
ax.axvline(x=0, color="k", linewidth=0.5)
plt.tight_layout()
plt.show()
