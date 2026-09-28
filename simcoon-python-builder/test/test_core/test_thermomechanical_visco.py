"""Thermomechanical viscoelasticity (PRONK_T): thermal strain in the branch driving force and
the heat source during relaxation."""

import numpy as np

import simcoon as sim
from simcoon.solver import StepThermomeca

E0, NU0, E1, NU1, ETA_B, ETA_S = 3000.0, 0.35, 1500.0, 0.35, 3000.0, 1200.0
RHO, C_P, ALPHA, T0 = 4.4, 0.656, 1e-4, 293.15
PROPS = np.array([RHO, C_P, E0, NU0, ALPHA, 1.0, E1, NU1, ETA_B, ETA_S])   # one Prony branch


def test_free_thermal_expansion_creates_no_viscous_strain():
    """Stress-free heating: the strain is alpha dT and the branches do not flow. The branch
    driving force used to omit the thermal strain, so it 'relaxed' the thermal expansion
    (13 % extra strain for this set)."""
    st = StepThermomeca(control="stress", value=[0.0] * 6, ninc=50, time=1.0, T_final=T0 + 50.0)
    r = sim.solver.solve([st], "PRONK", PROPS, 14, T_init=T0)
    np.testing.assert_allclose(r["Strain"][:3, -1], ALPHA * 50.0, rtol=1e-10)
    assert np.abs(r["Statev"][1:7, -1]).max() < 1e-14


def test_heat_released_during_a_strain_hold_balances_the_dissipation():
    """Strain ramp then hold at constant temperature: over the hold, the integrated heat source
    equals the dissipation minus the thermoelastic term T alpha : dsigma. It used to be a
    linearisation in (dE, dT), zero during the hold (and of the wrong sign here)."""
    uni = ["strain"] + ["stress"] * 5
    ramp = StepThermomeca(control=uni, value=[0.01, 0, 0, 0, 0, 0], ninc=20, time=0.1, T_final=T0)
    hold = StepThermomeca(control=uni, value=[0.01, 0, 0, 0, 0, 0], ninc=200, time=5.0, T_final=T0)
    r = sim.solver.solve([ramp, hold], "PRONK", PROPS, 14, T_init=T0)
    t = np.asarray(r["Time"])
    rate = np.asarray(r["r"]).ravel()
    i1 = int(np.argmin(np.abs(t - 0.1)))
    heat = np.sum(rate[i1 + 1:] * np.diff(t)[i1:])
    dWd = r["Wm"][3, -1] - r["Wm"][3, i1]
    dsig = r["Stress"][:, -1] - r["Stress"][:, i1]
    assert dWd > 0.0
    np.testing.assert_allclose(heat, dWd - T0 * ALPHA * np.sum(dsig[:3]), rtol=1e-8)
