"""
Constitutive Library Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# sphinx_gallery_thumbnail_path = 'images/thumb/sphx_glr_constitutive_relations_thumb.png'

###################################################################################
# This example demonstrates how to compute isotropic, cubic, transverse isotropic,
# and orthotropic stiffness and compliance tensors using Simcoon.
# It also shows how to check symmetries and extract properties from each tensor.

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim

###################################################################################
# L_iso
# -----
# Provides the elastic stiffness tensor for an isotropic material.
# The two first arguments are a couple of elastic properties.
# The third argument specifies which couple has been provided.
# Returns a NumPy ndarray.

E = 70000.0
nu = 0.3
L = sim.L_iso([E, nu], "Enu")
print("L_iso:\n", np.array_str(L, precision=2, suppress_small=True))

d = sim.check_symetries(L, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Properties:", d["props"])

x = sim.L_iso_props(L)
print("L_iso_props:", x)

###################################################################################
# M_iso
# -----
# Provides the elastic compliance tensor for an isotropic material.
# The two first arguments are a couple of elastic properties.
# The third argument specifies which couple has been provided.
# Returns a NumPy ndarray.

M = sim.M_iso([E, nu], "Enu")
print("M_iso:\n", np.array_str(M, suppress_small=True))

L_inv = np.linalg.inv(M)
d = sim.check_symetries(L_inv, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Properties:", d["props"])

x = sim.M_iso_props(M)
print("M_iso_props:", x)

###################################################################################
# L_cubic
# -------
# Provides the elastic stiffness tensor for a cubic material.
# Arguments are the stiffness coefficients C11, C12, C44, or elastic constants E, nu, G.

G = 23000.0
L = sim.L_cubic([E, nu, G], "EnuG")
print("L_cubic:\n", np.array_str(L, precision=2, suppress_small=True))

d = sim.check_symetries(L, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Properties:", d["props"])

x = sim.L_cubic_props(L)
print("L_cubic_props:", x)

###################################################################################
# M_cubic
# -------
# Provides the elastic compliance tensor for a cubic material.
# Arguments are the stiffness coefficients C11, C12, C44, or elastic constants E, nu, G.

M = sim.M_cubic([E, nu, G], "EnuG")
print("M_cubic:\n", np.array_str(M, suppress_small=True))

L = np.linalg.inv(M)
d = sim.check_symetries(L, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Properties:", d["props"])

x = sim.M_cubic_props(M)
print("M_cubic_props:", x)

###################################################################################
# L_isotrans
# ----------
# Provides the elastic stiffness tensor for a transverse isotropic material.
# Arguments: longitudinal modulus EL, transverse modulus ET, Poisson ratios nuTL, nuTT,
# shear modulus GLT, and axis of symmetry.

EL = 70000.0
ET = 20000.0
nuTL = 0.08
nuTT = 0.3
GLT = 12000.0
axis = 3

L = sim.L_isotrans([EL, ET, nuTL, nuTT, GLT], axis)
print("L_isotrans:\n", np.array_str(L, precision=2, suppress_small=True))

d = sim.check_symetries(L, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Axis:", d["axis"])
print("Properties:", np.array_str(d["props"], precision=2))

x = sim.L_isotrans_props(L, axis)
print("L_isotrans_props:", np.array_str(x, precision=2))

###################################################################################
# M_isotrans
# ----------
# Provides the elastic compliance tensor for a transverse isotropic material.

M = sim.M_isotrans([EL, ET, nuTL, nuTT, GLT], axis)
print("M_isotrans:\n", np.array_str(M, suppress_small=True))

x = sim.M_isotrans_props(M, axis)
print("M_isotrans_props:", np.array_str(x, precision=2))

###################################################################################
# L_ortho
# -------
# Provides the elastic stiffness tensor for an orthotropic material.
# Arguments: Young moduli E1, E2, E3; Poisson ratios nu12, nu13, nu23;
# shear moduli G12, G13, G23; or alternatively the list of Cii coefficients.

E_1, E_2, E_3 = 4500.0, 2300.0, 2700.0
nu_12, nu_13, nu_23 = 0.06, 0.08, 0.3
G_12, G_13, G_23 = 2200.0, 2100.0, 2400.0

L = sim.L_ortho([E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23], "EnuG")
print("L_ortho:\n", np.array_str(L, precision=2, suppress_small=True))

d = sim.check_symetries(L, 1.0e-2)
print("Symmetry type:", d["umat_type"])
print("Axis:", d["axis"])
print("Properties:", np.array_str(d["props"], precision=2))

x = sim.L_ortho_props(L)
print("L_ortho_props:", np.array_str(x, precision=2, suppress_small=True))

###################################################################################
# M_ortho
# -------
# Provides the elastic compliance tensor for an orthotropic material.

M = sim.M_ortho([E_1, E_2, E_3, nu_12, nu_13, nu_23, G_12, G_13, G_23], "EnuG")
print("M_ortho:\n", np.array_str(M, suppress_small=True))

x = sim.M_ortho_props(M)
print("M_ortho_props:", np.array_str(x, precision=4, suppress_small=True))
