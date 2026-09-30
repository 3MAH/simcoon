"""
Shape Memory Alloy - Transformation and Reorientation
=====================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
from simcoon.solver import Block, StepMeca

plt.rcParams["figure.figsize"] = (18, 10)

###################################################################################
# The ``SMRAI`` constitutive law adds a third rate-independent mechanism to the
# superelastic transformation model of :ref:`SMA_T <sphx_glr_gallery_mechanical_SMA_T.py>`:
# the reorientation of the martensite variants, which is what a non-proportional
# stress path activates once the material has transformed. Its four variants are
#
# - ``SMRDI`` / ``SMRDC``: isotropic / cubic elasticity, Drucker criteria,
# - ``SMRAI`` / ``SMRAC``: isotropic / cubic elasticity, anisotropic Drucker
#   (Deshpande-Fleck-Ashby) criteria.
#
# The properties are those of the transformation-only law, followed by the seven
# DFA anisotropy parameters :math:`F, G, H, L, M, N, K` for the ``SMRA*`` variants,
# then seven reorientation parameters:
#
# 1. :math:`Y^{Reo}` : stress limit for the onset of reorientation
# 2. :math:`H^{Reo}` : reorientation kinematic hardening modulus (MPa)
# 3. :math:`E_T^{Reo,max}` : maximum reorientation back-strain magnitude
# 4-7. :math:`c_{\lambda}, p_{0,\lambda}, n_{\lambda}, \alpha_{\lambda}` : penalty
#    function of the reorientation saturation
#
# With :math:`F = G = H = 1/2`, :math:`L = M = N = 3/2` and :math:`K = 0` the DFA
# operator is the von Mises one and ``SMRAI`` coincides with ``SMRDI``. Choosing
# :math:`Y^{Reo}` far above any reachable stress switches reorientation off, which
# recovers the transformation-only response of ``SMADI``.
#
# The thirty state variables hold the transformation block of ``SMADI`` (17 values),
# the cumulative reorientation multiplier :math:`p^{TR}`, the reorientation back-strain
# :math:`\mathbf{v}^{re}` (6) and the macroscopic reorientation strain
# :math:`\mathbf{E}^{Reo}` (6).

umat_name = "SMRAI"
nstatev = 30
T_init = 353.15  # K, above A_f: superelastic regime

# Transformation block: NiTi-like parameters (same order as in SMA_T.py)
flagT = 0
E_A, E_M = 70000.0, 70000.0
nu_A, nu_M = 0.3, 0.3
alphaA, alphaM = 1.0e-6, 1.0e-6
Hmin, Hmax, k1, sigmacrit = 0.0, 0.05, 0.021, 0.0
C_A, C_M = 6.0, 5.0
Ms0, Mf0, As0, Af0 = 293.15, 273.15, 313.15, 333.15
n1 = n2 = n3 = n4 = 0.2
sigmacaliber = 300.0
b_prager, n_prager = 1.4, 2.0
c_lambda, p0_lambda, n_lambda, alpha_lambda = 1.0e-6, 1.0e-3, 1.0, 1.0e8

props_T = [
    flagT, E_A, E_M, nu_A, nu_M, alphaA, alphaM,
    Hmin, Hmax, k1, sigmacrit,
    C_A, C_M, Ms0, Mf0, As0, Af0,
    n1, n2, n3, n4,
    sigmacaliber, b_prager, n_prager,
    c_lambda, p0_lambda, n_lambda, alpha_lambda,
]

# Anisotropic Drucker (DFA) parameters: F, G, H, L, M, N, K. These values reduce the
# criterion to the isotropic Drucker one; change them to introduce anisotropy.
props_DFA = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 0.0]

# Reorientation block. H_Reo is a modulus in MPa, calibrated on multiaxial data: keep it
# in the stable range of this law (H_Reo <= 5000 with these transformation parameters).
Y_Reo = 200.0
H_Reo = 5000.0
ETR_max = 0.05
props_Reo = [Y_Reo, H_Reo, ETR_max, 1.0e-6, 1.0e-3, 1.0, 1.0e8]

props = np.array(props_T + props_DFA + props_Reo)

###################################################################################
# Loading path
# ------------
#
# Reorientation needs a change of loading direction, so the path is a square in the
# :math:`(\varepsilon_{11}, \gamma_{12})` plane: tension first, then a full
# counter-clockwise box at constant amplitude, back to the origin. The strain
# components 11 and 12 are controlled, the four others are stress-free.

eps = 0.015
control = ["strain", "stress", "stress", "strain", "stress", "stress"]
corners = [(eps, 0), (eps, eps), (-eps, eps), (-eps, -eps), (eps, -eps), (eps, 0), (0, 0)]
steps = [
    StepMeca(control=control, value=[e11, 0, 0, g12, 0, 0], time=1.0, ninc=300, Dn_mini=0.01)
    for e11, g12 in corners
]
blocks = [Block(steps=steps)]


def run(name, props, nstatev):
    return sim.solver.solve(blocks, name, props, nstatev, T_init=T_init, corate=3)


res = run(umat_name, props, nstatev)

###################################################################################
# The transformation-only response on the same path is the reference to read the
# effect of reorientation against: the same law with reorientation switched off
# (:math:`Y^{Reo} = 10^{10}` MPa), which reproduces ``SMADI``.

props_off = np.array(props_T + props_DFA + [1.0e10] + props_Reo[1:])
res_off = run(umat_name, props_off, nstatev)

###################################################################################
# Plotting the results
# --------------------
#
# The stress path shows the distorted square typical of the multiaxial SMA
# experiments of Grabe and Bruhns: reorientation softens the corners where the
# loading direction turns. The state-variable plot shows the martensite fraction
# :math:`\xi`, the reorientation multiplier :math:`p^{TR}` and the norm of the
# reorientation strain :math:`\mathbf{E}^{Reo}`. The energy plot checks the split
# :math:`W_m = W_m^r + W_m^d` (:math:`W_m^{ir} = 0` for the SMA laws).

e11, _, _, g12, _, _ = res["Strain"]
s11, _, _, s12, _, _ = res["Stress"]
s11_off, _, _, s12_off, _, _ = res_off["Stress"]
time = res["Time"]
xi = res["Statev"][1]
pTR = res["Statev"][17]
EReo = np.linalg.norm(res["Statev"][24:30], axis=0)
Wm, Wm_r, Wm_ir, Wm_d = res["Wm"]

fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)
ax.grid(True)
ax.set_xlabel(r"$\sigma_{11}$ (MPa)", size=15)
ax.set_ylabel(r"$\sigma_{12}$ (MPa)", size=15)
ax.plot(s11_off, s12_off, c="gray", ls="--", label="transformation only")
ax.plot(s11, s12, c="blue", label="transformation + reorientation")
ax.legend(loc="best")

ax = fig.add_subplot(2, 2, 2)
ax.grid(True)
ax.set_xlabel(r"$\varepsilon_{11}$", size=15)
ax.set_ylabel(r"$\sigma_{11}$ (MPa)", size=15)
ax.plot(e11, s11_off, c="gray", ls="--", label="transformation only")
ax.plot(e11, s11, c="blue", label="transformation + reorientation")
ax.legend(loc="best")

ax = fig.add_subplot(2, 2, 3)
ax.grid(True)
ax.set_xlabel("time (s)", size=15)
ax.plot(time, xi, c="black", label=r"$\xi$")
ax.plot(time, pTR, c="red", label=r"$p^{TR}$")
ax.plot(time, EReo / Hmax, c="green", label=r"$\|\mathbf{E}^{Reo}\| / H_{max}$")
for k in range(1, len(corners)):
    ax.axvline(k, c="lightgray", lw=0.8)
ax.legend(loc="best")

ax = fig.add_subplot(2, 2, 4)
ax.grid(True)
ax.set_xlabel("time (s)", size=15)
ax.set_ylabel(r"$W_m$", size=15)
ax.plot(time, Wm, c="black", label=r"$W_m$")
ax.plot(time, Wm_r, c="green", label=r"$W_m^r$")
ax.plot(time, Wm_ir, c="blue", label=r"$W_m^{ir}$")
ax.plot(time, Wm_d, c="red", label=r"$W_m^d$")
ax.legend(loc="best")

plt.show()
