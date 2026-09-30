"""
Shape Memory Alloy - Transformation and Reorientation
=====================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import simcoon as sim
from simcoon.solver import Block, StepMeca

plt.rcParams["figure.figsize"] = (18, 8)

###################################################################################
# The ``SMRAI`` constitutive law adds a third rate-independent mechanism to the
# superelastic transformation model of :ref:`SMA_T <sphx_glr_gallery_mechanical_SMA_T.py>`:
# the reorientation of the martensite variants. It is what lets the model describe
# the whole thermomechanical map of a shape memory alloy, from the superelastic loop
# above :math:`A_f` to the shape memory effect below :math:`M_f`, where the
# self-accommodated martensite formed on cooling is oriented by the applied stress
# and the strain is recovered on heating. Its four variants are
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
# operator is the von Mises one and ``SMRAI`` coincides with ``SMRDI``.
#
# The thirty state variables hold the transformation block of ``SMADI`` (17 values),
# the cumulative reorientation multiplier :math:`p^{TR}`, the reorientation back-strain
# :math:`\mathbf{v}^{re}` (6) and the macroscopic reorientation strain
# :math:`\mathbf{E}^{Reo}` (6).

umat_name = "SMRAI"
nstatev = 30

# Transformation block: NiTi-like parameters (same order as in SMA_T.py)
flagT = 0
E_A, E_M = 70000.0, 70000.0
nu_A, nu_M = 0.3, 0.3
alphaA, alphaM = 1.0e-6, 1.0e-6
Hmin, Hmax, k1, sigmacrit = 0.0, 0.05, 0.021, 0.0
C_A, C_M = 6.0, 5.0  # Clausius-Clapeyron slopes (MPa/K)
Ms0, Mf0, As0, Af0 = 293.15, 273.15, 313.15, 333.15  # 20, 0, 40, 60 degrees C
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
# Loading protocol
# ----------------
#
# Every run starts in the austenitic state at 80 degrees C, above :math:`A_f`, and
# reaches its test temperature by a stress-free cooling step, so that the martensite
# present below :math:`M_s` is the self-accommodated one formed on cooling, with no
# macroscopic strain. Loading is uniaxial and strain-controlled, unloading is
# stress-controlled.

T_start = 273.15 + 80.0
zero = [0.0] * 6
stress_free = ["stress"] * 6
uniaxial = ["strain"] + ["stress"] * 5
eps_max = 0.04


def stress_free_step(T_final=None, ninc=200):
    """Stress-free step: a temperature ramp to T_final, or an unloading at constant temperature."""
    return StepMeca(control=stress_free, value=zero, time=1.0, ninc=ninc, T_final=T_final, Dn_mini=0.01)


def load_step(eps, ninc=200):
    return StepMeca(control=uniaxial, value=[eps, 0, 0, 0, 0, 0], time=1.0, ninc=ninc, Dn_mini=0.01)


def run(steps):
    return sim.solver.solve([Block(steps=steps)], umat_name, props, nstatev, T_init=T_start, corate=3)


###################################################################################
# From shape memory to superelasticity
# ------------------------------------
#
# Five isothermal tension-unloading cycles at temperatures spanning the transformation
# temperatures. Below :math:`M_f` the material is fully martensitic when loaded: the
# stress orients the variants (reorientation) and the strain remains after unloading.
# Between :math:`M_s` and :math:`A_f` the stress-induced transformation is only partly
# reversed on unloading and a residual strain remains. Above :math:`A_f` the
# transformation reverses completely: the superelastic loop closes.

temperatures = [-20.0, 10.0, 30.0, 50.0, 80.0]  # degrees C: < Mf, Mf-Ms, Ms-As, As-Af, > Af
sweep = {}
for Tc in temperatures:
    res = run([stress_free_step(273.15 + Tc), load_step(eps_max), stress_free_step()])
    e11 = res["Strain"][0]
    s11 = res["Stress"][0]
    loaded = res["Time"] >= 1.0  # drop the cooling step
    sweep[Tc] = (e11[loaded], s11[loaded], res["Statev"][1][loaded])

# One temperature colormap for the whole example: blue below Mf, red above Af
T_cmap = plt.cm.turbo
T_norm = plt.Normalize(-40.0, 110.0)

fig = plt.figure()
grid = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0])
colors = [T_cmap(T_norm(Tc)) for Tc in temperatures]

# stress-strain loops stacked along the temperature axis
ax = fig.add_subplot(grid[0], projection="3d")
for c, Tc in zip(colors, temperatures):
    e11, s11, _ = sweep[Tc]
    ax.plot(e11, np.full_like(e11, Tc), s11, c=c, lw=2, label=f"T = {Tc:.0f} °C")
ax.set_xlabel(r"$\varepsilon_{11}$", size=13)
ax.set_ylabel("Temperature (°C)", size=13)
ax.set_zlabel(r"$\sigma_{11}$ (MPa)", size=13)
ax.set_xlim(0.0, eps_max)
ax.set_ylim(temperatures[0] - 10.0, temperatures[-1] + 10.0)
ax.set_box_aspect((1.0, 1.6, 0.9))
ax.view_init(elev=22, azim=-40)
ax.legend(loc="upper left")

ax = fig.add_subplot(grid[1])
ax.grid(True)
ax.set_xlabel(r"Strain $\varepsilon_{11}$", size=15)
ax.set_ylabel(r"Martensite volume fraction $\xi$", size=15)
for c, Tc in zip(colors, temperatures):
    e11, _, xi = sweep[Tc]
    ax.plot(e11, xi, c=c, lw=2, label=f"T = {Tc:.0f} °C")
ax.legend(loc="best")
plt.show()

###################################################################################
# Shape memory effect in the stress-temperature-strain space
# ----------------------------------------------------------
#
# One cycle: cooling from 80 to -20 degrees C under no stress (self-accommodated
# martensite, no strain), tension to 4 % (variant reorientation), unloading (the
# oriented martensite keeps most of the strain), then heating to 100 degrees C
# under no stress: the reverse transformation recovers the strain.
#
# The path is drawn in the :math:`(T, \varepsilon, \sigma)` space, stress vertical,
# temperature horizontal and strain towards the viewer. Its projection on the
# strain-free back wall is the classical Clausius-Clapeyron diagram, with the four
# transformation lines :math:`\sigma = C_M (T - M_s)`, :math:`C_M (T - M_f)`,
# :math:`C_A (T - A_s)` and :math:`C_A (T - A_f)` drawn for reference: the loading
# branch crosses none of them (it is reorientation, not transformation), the heating
# branch crosses the :math:`A_s` and :math:`A_f` lines at zero stress.

res = run([stress_free_step(273.15 - 20.0), load_step(eps_max), stress_free_step(),
           stress_free_step(273.15 + 100.0, ninc=400)])
T = res["Temp"] - 273.15
s11 = res["Stress"][0]
e11 = res["Strain"][0]
xi = res["Statev"][1]
time = res["Time"]

from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.collections import LineCollection


def coloured_segments(*coords):
    """Consecutive-point segments of a polyline, to colour each by its temperature."""
    pts = np.stack(coords, axis=-1)
    return np.stack([pts[:-1], pts[1:]], axis=1)


fig = plt.figure(figsize=(18, 9))
grid = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0])
ax = fig.add_subplot(grid[0], projection="3d")
T_min, T_max, s_max = -40.0, 110.0, 450.0
T_line = np.linspace(T_min, T_max, 100)
for T0, C, label in [(Mf0, C_M, r"$M_f$"), (Ms0, C_M, r"$M_s$"), (As0, C_A, r"$A_s$"), (Af0, C_A, r"$A_f$")]:
    sig_line = C * (T_line - (T0 - 273.15))
    keep = (sig_line >= 0.0) & (sig_line <= s_max)
    ax.plot(T_line[keep], np.zeros(keep.sum()), sig_line[keep], c="black", lw=1, ls="--")
    ax.text(T_line[keep][-1], 0.0, sig_line[keep][-1], label, color="black", size=12)
# projection of the cycle on the strain-free back wall: the classical stress-temperature view
ax.plot(T, np.zeros_like(T), s11, c="lightgray", lw=1)
path = Line3DCollection(coloured_segments(T, e11, s11), cmap=T_cmap, norm=T_norm, lw=2.5)
path.set_array(0.5 * (T[:-1] + T[1:]))
ax.add_collection(path)
for label, t_mid in [("cooling", 0.5), ("loading", 1.5), ("unloading", 2.5), ("heating", 3.6)]:
    k = np.argmin(np.abs(time - t_mid))
    ax.text(T[k], e11[k], s11[k], "  " + label, size=11)
ax.set_xlim(T_min, T_max)
ax.set_ylim(eps_max, 0.0)  # strain grows towards the viewer
ax.set_zlim(0.0, s_max)
ax.set_xlabel("Temperature (°C)", size=13)
ax.set_ylabel(r"$\varepsilon_{11}$", size=13)
ax.set_zlabel(r"$\sigma_{11}$ (MPa)", size=13)
ax.set_box_aspect((1.7, 1.0, 1.0))
ax.view_init(elev=18, azim=-70)
fig.colorbar(path, ax=ax, orientation="horizontal", shrink=0.5, pad=0.08, label="Temperature (°C)")

ax = fig.add_subplot(grid[1])
ax.grid(True)
ax.set_xlabel("time (s)", size=15)
strain_line = LineCollection(coloured_segments(time, e11 / eps_max), cmap=T_cmap, norm=T_norm, lw=2.5)
strain_line.set_array(0.5 * (T[:-1] + T[1:]))
ax.add_collection(strain_line)
ax.plot([], [], c="black", lw=2.5, label=r"$\varepsilon_{11} / \varepsilon_{max}$ (coloured by T)")
ax.plot(time, s11 / s11.max(), c="black", ls="--", label=r"$\sigma_{11} / \sigma_{max}$")
ax.plot(time, xi, c="gray", label=r"$\xi$")
for t in (1.0, 2.0, 3.0):
    ax.axvline(t, c="lightgray", lw=0.8)
ax.set_xlim(time[0], time[-1])
ax.set_ylim(-0.02, 1.05)
ax.legend(loc="best")
plt.show()
