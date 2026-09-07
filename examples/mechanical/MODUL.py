"""
Modular UMAT Example — Composable Elasto-Plasticity
====================================================

Demonstrates the ``simcoon.modular`` high-level Python interface that composes a
constitutive model declaratively and runs it through the in-memory
:func:`simcoon.solver.solve` driver. The C++ ``ModularUMAT`` infrastructure
(ElasticityModule, YieldCriterion, hardening, ...) is internal — the user only
builds a ``ModularMaterial`` and hands its ``.props`` / ``.nstatev`` to the
``"MODUL"`` UMAT code registered in simcoon's UMAT table.

The loading applies a monotonic tensile ramp to 2% strain, then two
strain-controlled cycles between -2% and +2%, exposing isotropic-hardening
growth and the initial yield plateau.

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
# Isotropic elasticity + von Mises yield + Voce isotropic hardening.
# Parameters: E=210 GPa, nu=0.3, sigma_Y=300 MPa, Q=200 MPa, b=10.
#
# The elastic constants C1/C2 are ordinal slots whose meaning is fixed by the
# ``convention`` argument — here ``"Enu"`` (C1 = E, C2 = nu). Other
# parameterizations ("Kmu", "lambdamu", ...) are accepted as well.

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
# Two blocks, both small-strain and uniaxial (axial strain driven, the five
# other components stress-free): a monotonic ramp to +2%, then a two-cycle
# block alternating -2% / +2%. A ``Block`` repeats its step sequence
# ``ncycle`` times, which is how the cycling is expressed.

uniaxial = ["strain"] + ["stress"] * 5


def ramp(target, ninc):
    return solver.StepMeca(
        control=uniaxial, value=[target, 0, 0, 0, 0, 0],
        time=1.0, ninc=ninc, Dn_mini=1.0,
    )


path = [
    solver.Block(steps=[ramp(0.02, 200)]),
    solver.Block(steps=[ramp(-0.02, 100), ramp(0.02, 100)], ncycle=2),
]

###################################################################################
# 3. Run the solver
# ------------------
# ``mat.umat_name`` is ``"MODUL"``, the UMAT code registered at
# ``umat_smart.cpp:316`` (id 200). ``mat.props`` serializes the composition into
# the flat array that the C++ ``umat_modular`` deserializes.

res = solver.solve(path, mat.umat_name, mat.props, mat.nstatev, T_init=293.0)

###################################################################################
# 4. Plot the stress-strain curve
# --------------------------------
# Results come back as numpy arrays in a components-first layout: ``[0]`` is
# the 11 component of each tensor history.

e11 = res["Strain"][0]
s11 = res["Stress"][0]
Wm, Wm_r, Wm_ir, _ = res["Wm"]

fig = plt.figure()

# Stress-strain curve
ax1 = fig.add_subplot(1, 2, 1)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=13)
plt.xlabel(r"Strain $\varepsilon_{11}$", size=14)
plt.ylabel(r"Stress $\sigma_{11}$ (MPa)", size=14)
plt.plot(e11, s11, c="royalblue", lw=1.5, label="MODUL: iso-elastic + VM + Voce")
plt.axhline(y=300.0, color="0.6", linestyle="--", lw=0.8, label=r"initial $\sigma_Y$")
plt.axhline(y=-300.0, color="0.6", linestyle="--", lw=0.8)
plt.legend(loc="best")
plt.title("Stress-strain response")

# Work terms vs time
ax2 = fig.add_subplot(1, 2, 2)
plt.grid(True)
plt.tick_params(axis="both", which="major", labelsize=13)
plt.xlabel("time (s)", size=14)
plt.ylabel("Work (MPa)", size=14)
plt.plot(res["Time"], Wm, c="black", label=r"$W_m$ (total)")
plt.plot(res["Time"], Wm_r, c="green", label=r"$W_m^r$ (recoverable)")
plt.plot(res["Time"], Wm_ir, c="blue", label=r"$W_m^{ir}$ (irreversible)")
plt.legend(loc="best")
plt.title("Energy decomposition")

plt.tight_layout()
plt.savefig("MODUL_stress_strain.png", dpi=120)
plt.show()
