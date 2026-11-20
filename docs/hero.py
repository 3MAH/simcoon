"""Hero image of the documentation and the README, computed with simcoon.

Run ``python docs/hero.py`` to regenerate ``docs/_static/hero.png``.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import simcoon as sim
from simcoon import solver
from simcoon.solver.micromechanics import Ellipsoid, MoriTanaka, SelfConsistent

RED = "#a61b29"
NAVY = "#24445c"


def elastoplastic_cycles():
    props = np.array([200000.0, 0.3, 1.0e-5, 200.0, 100.0, 10.0, 60000.0, 400.0, 5000.0, 20.0])
    mixed = ["strain"] + ["stress"] * 5
    steps = []
    for amplitude in (0.004, 0.008, 0.012):
        steps.append(solver.StepMeca(control=mixed, value=[amplitude, 0, 0, 0, 0, 0], ninc=150))
        steps.append(solver.StepMeca(control=mixed, value=[-amplitude, 0, 0, 0, 0, 0], ninc=300))
    res = solver.solve(solver.Block(steps=steps), "EPCHA", props, nstatev=33)
    return 100 * res["Strain"][0], res["Stress"][0]


def superelastic_cycle():
    props = np.array([
        5.5, 0.5, 0.5, 0,
        67538.0, 67538.0, 0.349, 0.349,
        1.0e-6, 1.0e-6,
        0.0, 0.0418, 0.021, 0.0,
        10.0, 10.0,
        300.0, 290.0, 295.0, 305.0,
        0.2, 0.2, 0.2, 0.2,
        300.0, 0.2, 2.0,
        1.0e-6, 1.0e-5, 1.0, 0.0,
    ])
    mixed = ["strain"] + ["stress"] * 5
    loading = solver.StepThermomeca(control=mixed, value=[0.06, 0, 0, 0, 0, 0], ninc=300,
                                    time=30.0, thermal_control="heat_flux", Q=0.0)
    unloading = solver.StepThermomeca(control=mixed, value=[0.0, 0, 0, 0, 0, 0], ninc=300,
                                      time=30.0, thermal_control="heat_flux", Q=0.0)
    res = solver.solve(solver.Block(steps=[loading, unloading]), "SMADI", props,
                       nstatev=17, T_init=323.15)
    return 100 * res["Strain"][0], res["Stress"][0], res["Temp"]


def effective_modulus():
    matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, concentration=1.0, nstatev=1,
                       props=np.array([2250.0, 0.19, 8.8e-5]))
    fibre = Ellipsoid(number=1, umat_name="ELISO", save=1, concentration=0.0, nstatev=1,
                      props=np.array([73000.0, 0.19, 0.5e-6]))
    c = np.linspace(0.0, 0.5, 26)
    curves = {}
    for name, scheme in (("Mori-Tanaka", MoriTanaka(50, 50, 0)),
                         ("self-consistent", SelfConsistent(50, 50, 0))):
        E = np.zeros_like(c)
        for i, x in enumerate(c):
            fibre.concentration, matrix.concentration = x, 1.0 - x
            L = sim.L_eff(scheme, nstatev=0, orientation=(0.0, 0.0, 0.0),
                          phases=[matrix, fibre])
            E[i] = sim.L_iso_props(L).flatten()[0]
        curves[name] = E
    return c, curves


def main(path):
    with plt.rc_context({"text.usetex": True, "font.family": "serif", "font.size": 13}):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

        e, s = elastoplastic_cycles()
        axes[0].plot(e, s, color=RED, lw=1.8)
        axes[0].set_xlabel(r"$\varepsilon_{11}$ (\%)")
        axes[0].set_ylabel(r"$\sigma_{11}$ (MPa)")
        axes[0].set_title("Cyclic plasticity (Chaboche)")

        e, s, T = superelastic_cycle()
        points = axes[1].scatter(e, s, c=T, s=6, cmap="coolwarm")
        axes[1].set_xlabel(r"$\varepsilon_{11}$ (\%)")
        axes[1].set_ylabel(r"$\sigma_{11}$ (MPa)")
        axes[1].set_title("Adiabatic superelasticity")
        fig.colorbar(points, ax=axes[1], label=r"$T$ (K)")

        c, curves = effective_modulus()
        for (name, E), color in zip(curves.items(), (NAVY, RED)):
            axes[2].plot(c, E / 1000.0, color=color, lw=2, label=name)
        axes[2].set_xlabel(r"glass volume fraction $c$")
        axes[2].set_ylabel(r"$E^{\mathrm{eff}}$ (GPa)")
        axes[2].set_title("Mean-field homogenization")
        axes[2].legend(frameon=False)

        for ax in axes:
            ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(path, dpi=110, facecolor="white")
        plt.close(fig)


if __name__ == "__main__":
    main(os.path.join(os.path.dirname(os.path.abspath(__file__)), "_static", "hero.png"))
