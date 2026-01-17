"""
Eshelby Tensor Examples
~~~~~~~~~~~~~~~~~~~~~~~~~
This example demonstrates Eshelby tensor computations for various inclusion shapes.
"""

import numpy as np
from simcoon import simmit as sim
import matplotlib.pyplot as plt

###################################################################################
# The Eshelby tensor :math:`\mathbf{S}` relates the constrained strain inside an
# ellipsoidal inclusion to the eigenstrain (stress-free strain) when the inclusion
# is embedded in an infinite elastic matrix.
#
# For an eigenstrain :math:`\varepsilon^*`, the constrained strain inside the inclusion is:
#
# .. math::
#
#    \varepsilon^c = \mathbf{S} : \varepsilon^*
#
# The Eshelby tensor depends on the inclusion shape (aspect ratios) and the matrix
# Poisson ratio. Simcoon provides analytical solutions for special cases (sphere,
# cylinder, prolate, oblate) and numerical integration for general ellipsoids.

###################################################################################
# Eshelby tensor for a sphere
# ---------------------------
# For a spherical inclusion, the Eshelby tensor has a simple analytical form
# that depends only on the Poisson ratio :math:`\nu` of the matrix.

nu = 0.3  # Poisson ratio of the matrix

S_sphere = sim.Eshelby_sphere(nu)
print("Eshelby tensor for sphere (nu=0.3):")
print(np.array_str(S_sphere, precision=4, suppress_small=True))

###################################################################################
# Eshelby tensor for an infinite cylinder
# ---------------------------------------
# For a cylindrical inclusion (infinite along axis 1), the tensor also has
# an analytical form.

S_cylinder = sim.Eshelby_cylinder(nu)
print("\nEshelby tensor for infinite cylinder (nu=0.3):")
print(np.array_str(S_cylinder, precision=4, suppress_small=True))

###################################################################################
# Eshelby tensor for prolate ellipsoid
# ------------------------------------
# A prolate ellipsoid is elongated along one axis (a1 > a2 = a3).
# The aspect ratio is defined as ar = a1/a3.

ar = 5.0  # Aspect ratio a1/a3
S_prolate = sim.Eshelby_prolate(nu, ar)
print(f"\nEshelby tensor for prolate ellipsoid (ar={ar}, nu={nu}):")
print(np.array_str(S_prolate, precision=4, suppress_small=True))

###################################################################################
# Eshelby tensor for oblate ellipsoid
# -----------------------------------
# An oblate ellipsoid is flattened along one axis (a1 < a2 = a3).
# The aspect ratio is defined as ar = a1/a3 (< 1).

ar = 0.2  # Aspect ratio a1/a3
S_oblate = sim.Eshelby_oblate(nu, ar)
print(f"\nEshelby tensor for oblate ellipsoid (ar={ar}, nu={nu}):")
print(np.array_str(S_oblate, precision=4, suppress_small=True))

###################################################################################
# Eshelby tensor for penny-shaped crack
# -------------------------------------
# A penny-shaped (crack) inclusion represents the limiting case of an oblate
# ellipsoid when the aspect ratio :math:`ar \to 0`. This corresponds to a flat
# disc or crack perpendicular to axis 1.
#
# The analytical limit has:
#
# - :math:`S_{1111} = 1`
# - :math:`S_{1122} = S_{1133} = \frac{\nu}{1-\nu}`
# - :math:`S_{1212} = S_{1313} = \frac{1}{2}` (factor 2 in Voigt notation)

S_penny = sim.Eshelby_penny(nu)
print(f"\nEshelby tensor for penny-shaped crack (nu={nu}):")
print(np.array_str(S_penny, precision=4, suppress_small=True))

# Verify oblate converges to penny limit
print("\nConvergence of oblate to penny limit:")
for ar_test in [0.1, 0.01, 0.001]:
    S_obl = sim.Eshelby_oblate(nu, ar_test)
    print(f"  ar={ar_test}: S_1111={S_obl[0, 0]:.4f}, S_1122={S_obl[0, 1]:.4f}")
print(f"  Penny limit: S_1111={S_penny[0, 0]:.4f}, S_1122={S_penny[0, 1]:.4f}")

###################################################################################
# Numerical Eshelby tensor for general ellipsoid
# ----------------------------------------------
# For anisotropic matrices or general ellipsoidal shapes, Simcoon uses numerical
# integration over the unit sphere.

# Define an isotropic stiffness tensor
E = 70000.0  # Young's modulus (MPa)
L = sim.L_iso([E, nu], "Enu")

# Ellipsoid semi-axes
a1, a2, a3 = 1.0, 1.0, 1.0  # Sphere for verification
mp, np_int = 50, 50  # Integration points

S_numerical = sim.Eshelby(L, a1, a2, a3, mp, np_int)
print("\nNumerical Eshelby tensor for sphere:")
print(np.array_str(S_numerical, precision=4, suppress_small=True))

# Verify against analytical solution
error = np.linalg.norm(S_numerical - S_sphere)
print(f"\nError vs analytical sphere solution: {error:.2e}")

###################################################################################
# Hill's interaction tensor
# -------------------------
# The Hill interaction tensor :math:`\mathbf{T}` is related to the Eshelby tensor
# through :math:`\mathbf{S} = \mathbf{T} : \mathbf{L}`, where :math:`\mathbf{L}`
# is the stiffness tensor of the matrix.

T_II = sim.T_II(L, a1, a2, a3, mp, np_int)
print("\nHill interaction tensor T_II:")
print(np.array_str(T_II, precision=6, suppress_small=True))

# Verify the relation S = T:L
S_from_T = T_II @ L
error_TL = np.linalg.norm(S_from_T - S_numerical)
print(f"\nVerification S = T:L, error: {error_TL:.2e}")

###################################################################################
# .. seealso::
#
#    - :ref:`sphx_glr_examples_heterogeneous_effective_props.py` - Micromechanics
#      example showing effective properties vs volume fraction and aspect ratio
#    - :ref:`sphx_glr_examples_analysis_eshelby_numerical_vs_analytical.py` -
#      Detailed comparison of numerical vs analytical Eshelby tensor computations
