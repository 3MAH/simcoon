"""
2 phase composite homogenization example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial studies the evolution of the mechanical properties of a
composite material considering spherical and ellipsoidal reinforcements.
"""

import numpy as np
import simcoon as sim
import matplotlib.pyplot as plt


###################################################################################
# In this tutorial we will study mechanical properties of a 2-phase composite
# material considering spherical reinforcement. We will also explore the effect
# of inclusion aspect ratio on both the Eshelby tensor and the effective properties.
#
# First we define the number of state variables (if any) at the macroscopic level
# and the material properties.

nstatev = 0  # None here

int1 = 50  # Integration points of the Eshelby integrals, first direction
int2 = 50  # Integration points of the Eshelby integrals, second direction
n_matrix = 0  # Index of the matrix phase in the list of phases

# Mori-Tanaka props: [int1, int2, n_matrix]; the phases themselves are passed as objects
props = np.array([int1, int2, n_matrix], dtype="float")

###############################################################################
# The two phases are described in memory, as objects: a matrix and a spherical
# reinforcement, each with its own constitutive model and properties.

from simcoon.solver.micromechanics import Ellipsoid

matrix = Ellipsoid(
    number=0, umat_name="ELISO", save=1, concentration=0.8, nstatev=1,
    props=np.array([2250.0, 0.19, 8.8e-5]),
)
reinforcement = Ellipsoid(
    number=1, umat_name="ELISO", save=1, concentration=0.2, nstatev=1,
    props=np.array([73000.0, 0.19, 0.5e-6]),
)
phases = [matrix, reinforcement]

###############################################################################
# There is a possibility to consider a misorientation between the test frame
# of reference and the composite material orientation, using Euler angles
# with the z–x–z convention. The stiffness tensor will be returned in the
# test frame of reference.
#
# We also specify which micromechanical scheme will be used:
#
# - `MIMTN` → Mori-Tanaka scheme
# - `MISCN` → Self-consistent scheme

psi_rve = 0.0
theta_rve = 0.0
phi_rve = 0.0

umat_name = "MIMTN"  # Micromechanical scheme (Mori-Tanaka here)

###############################################################################
# Run the simulation to obtain the effective isotropic elastic properties.
# Since both phases are isotropic and the reinforcements are spherical, the
# result is also isotropic.

L_eff = sim.L_eff(
    umat_name, props, nstatev, orientation=(psi_rve, theta_rve, phi_rve), phases=phases
)
p = sim.L_iso_props(L_eff).flatten()
np.set_printoptions(precision=3, suppress=True)

print("L_eff = ", L_eff)
print("Effective isotropic properties (E, nu):", p)

###############################################################################
# Note:
# Even though both materials share the same Poisson ratio, the stiffness
# mismatch leads to a different *effective* Poisson ratio.

###############################################################################
# Effect of aspect ratio on Eshelby tensor components
# ---------------------------------------------------
# The Eshelby tensor :math:`\mathbf{S}` relates the constrained strain inside an
# ellipsoidal inclusion to the eigenstrain. It depends on the inclusion shape
# and the matrix Poisson ratio. Let's visualize how its components vary with
# aspect ratio.

nu = 0.3  # Poisson ratio of the matrix

aspect_ratios = np.logspace(-1, 1, 50)  # From 0.1 to 10
S_11 = []
S_22 = []
S_33 = []

for ar in aspect_ratios:
    if ar > 1:
        S = sim.Eshelby_prolate(nu, ar)
    else:
        S = sim.Eshelby_oblate(nu, ar)
    S_11.append(S[0, 0])
    S_22.append(S[1, 1])
    S_33.append(S[2, 2])

fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogx(aspect_ratios, S_11, "b-", label=r"$S_{1111}$", linewidth=2)
ax.semilogx(aspect_ratios, S_22, "r--", label=r"$S_{2222}$", linewidth=2)
ax.semilogx(aspect_ratios, S_33, "g-.", label=r"$S_{3333}$", linewidth=2)
ax.axvline(x=1, color="k", linestyle=":", alpha=0.5, label="Sphere (ar=1)")
ax.set_xlabel("Aspect ratio $a_1/a_3$", fontsize=12)
ax.set_ylabel("Eshelby tensor components", fontsize=12)
ax.set_title(
    f"Eshelby tensor diagonal components vs aspect ratio ($\\nu$={nu})", fontsize=14
)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
