"""
2 phase composite homogenization example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This tutorial studies the evolution of the mechanical properties of a
composite material considering spherical reinforcement.
"""

import numpy as np
from simcoon import simmit as sim
import os


###################################################################################
# In this tutorial we will study mechanical properties of a 2-phase composite material considering spherical reinforcement (see figure below)
#
# First we define the number of state variables (if any) at the macroscopic level
# and the material properties.

nstatev = 0  # None here

nphases = 2  # Number of phases
num_file = 0  # Index of the file that contains the subphases
int1 = 50  # Number of integration points along the long axis
int2 = 50  # Number of integration points along the lat axis
n_matrix = 0  # Phase number for the matrix

props = np.array([nphases, num_file, int1, int2, n_matrix], dtype="float")

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

L_eff = sim.L_eff(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve)
p = sim.L_iso_props(L_eff).flatten()
np.set_printoptions(precision=3, suppress=True)

print("L_eff = ", L_eff)
print("Effective isotropic properties (E, nu:", p)

###############################################################################
# Note:
# Even though both materials share the same Poisson ratio, the stiffness
# mismatch leads to a different *effective* Poisson ratio.
