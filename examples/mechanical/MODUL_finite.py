"""
Modular UMAT under Finite Strain — Hencky Hyperelasto-Plasticity
=================================================================

Runs the same composable ``"MODUL"`` material as ``MODUL.py``, but under
finite strain (NLGEOM): control type ``"logarithmic"`` drives the logarithmic
strain / Kirchhoff stress conjugate pair, and the solver kinematics hand the
model the log strain of the actual deformation gradient. The elasticity block
then acts as a Hencky stored-energy function of ln V and the plasticity
mechanism rides additively on that measure — a genuine hyperelasto-plastic
model.

This holds only when the accumulated corotational strain is exactly ln V,
which is the log_R corate: ``corate="logarithmic_R"`` is REQUIRED for MODUL
under NLGEOM and any other corate raises a ``RuntimeError`` (a Jaumann- or
Green-Naghdi-integrated model would be hypoelastic and dissipate spuriously
in closed cycles).

The loading is a log-strain cycle +15% / -15% / 0 — genuinely finite
stretches (lambda from 0.86 to 1.16) — exposing the elasto-plastic
hysteresis loop in the (ln V, tau) work-conjugate plane.

Both the material and the loading path are built in Python: no ``path.txt``,
no ``material.dat``, no result file on disk.
"""

import matplotlib.pyplot as plt
from simcoon import solver
from simcoon.modular import (
    ModularMaterial,
    IsotropicElasticity,
    Plasticity,
    VonMisesYield,
    VoceHardening,
)

plt.rcParams["figure.figsize"] = (14, 6)

###################################################################################
# 1. Compose the constitutive model
# ----------------------------------
# Same composition as MODUL.py: isotropic elasticity ("Enu": C1 = E, C2 = nu)
# + von Mises yield + Voce isotropic hardening.

mat = ModularMaterial(
    elasticity=IsotropicElasticity(
        C1=210000.0, C2=0.3, alpha=1.2e-5, convention="Enu"
    ),
    mechanisms=[
        Plasticity(
            sigma_Y=300.0,
            yield_criterion=VonMisesYield(),
            isotropic_hardening=VoceHardening(Q=200.0, b=10.0),
        ),
    ],
)

print(mat.summary())

###################################################################################
# 2. Build the loading path
# --------------------------
# Three steps driving the axial log strain to +15%, back to -15% and finally
# to 0, the five other components held stress-free (uniaxial tension /
# compression). Each step covers 1 s in 100 increments, with adaptive
# sub-stepping down to Dn_mini = 1e-3 of an increment.

steps = [
    solver.StepMeca(
        control=["strain"] + ["stress"] * 5,
        value=[target, 0, 0, 0, 0, 0],
        time=1.0,
        ninc=100,
        Dn_mini=1.0e-3,
    )
    for target in (0.15, -0.15, 0.0)
]

###################################################################################
# 3. Run the solver under NLGEOM
# -------------------------------
# ``control_type="logarithmic"`` is the finite-strain (log strain / Kirchhoff
# stress) control, and log_R is the only hyper/hypo-consistent corate for
# MODUL. It is also the solver default, but spelled out here because the
# model depends on it.

res = solver.solve(
    solver.Block(steps=steps, control_type="logarithmic"),
    mat.umat_name,
    mat.props,
    mat.nstatev,
    T_init=293.0,
    corate="logarithmic_R",
)

###################################################################################
# 4. Plot the response
# ---------------------
# The in-memory results carry every stress and strain measure the solver
# integrated, so the model's own work-conjugate pair is read directly:
# ``LogStrain`` is ln V and ``Kirchhoff`` is tau. (The legacy file output
# reported Green-Lagrange strain and Cauchy stress, which had to be converted
# by hand — and only exactly so on a rotation-free path such as this one.)

e11 = res["LogStrain"][0]
tau11 = res["Kirchhoff"][0]
Wm, Wm_r, _, Wm_d = res["Wm"]

fig = plt.figure()

ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=13)
plt.xlabel(r"Log strain $(\ln \mathbf{V})_{11}$", size=14)
plt.ylabel(r"Kirchhoff stress $\tau_{11}$ (MPa)", size=14)
plt.plot(e11, tau11, c="royalblue", lw=1.5,
         label="MODUL, NLGEOM logarithmic, corate log_R")
plt.axhline(y=300.0, color="0.6", linestyle="--", lw=0.8,
            label=r"initial $\sigma_Y$")
plt.axhline(y=-300.0, color="0.6", linestyle="--", lw=0.8)
plt.legend(loc="best")
plt.title("Finite-strain hysteresis loop")

ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=13)
plt.xlabel("time (s)", size=14)
plt.ylabel("Work (MPa)", size=14)
plt.plot(res["Time"], Wm, c="black", label=r"$W_m$ (total)")
plt.plot(res["Time"], Wm_r, c="green", label=r"$W_m^r$ (recoverable)")
plt.plot(res["Time"], Wm_d, c="red", label=r"$W_m^d$ (dissipated)")
plt.legend(loc="best")
plt.title("Energy decomposition")

plt.tight_layout()
plt.savefig("MODUL_finite_stress_strain.png", dpi=120)
plt.show()
