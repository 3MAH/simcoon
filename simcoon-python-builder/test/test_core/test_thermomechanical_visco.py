"""Thermomechanical viscoelasticity (PRONK_T): thermal strain in the branch driving force and
the heat source during relaxation."""

import numpy as np
import pytest

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


# ----- dissipation of the Kelvin branches in creep and recovery ------------------------------

_T0 = 293.15
_CREEP = {  # rho, c_p, E0, nu0, alpha, (N,) E1, nu1, etaB, etaS -- alpha = 0: isothermal balance
    "ZENER": (np.array([4.4, 0.656, 3000., 0.35, 0., 1500., 0.35, 3000., 1200.]), 8),
    "ZENNK": (np.array([4.4, 0.656, 3000., 0.35, 0., 1., 1500., 0.35, 3000., 1200.]), 14),
    "PRONK": (np.array([4.4, 0.656, 3000., 0.35, 0., 1., 1500., 0.35, 3000., 1200.]), 14),
}


@pytest.mark.parametrize("umat", sorted(_CREEP))
def test_creep_recovery_dissipation_is_nonnegative(umat):
    """Stress held, released, held again: the branches relax both ways, Wm_d never decreases
    (the Zener_T kernels used the elastic L EV instead of the viscous stress and lost
    dissipation in recovery), and the heat source integrates to Wm_d at constant T."""
    props, nstatev = _CREEP[umat]

    def step(v, t, n):
        return StepThermomeca(control="stress", value=[v, 0, 0, 0, 0, 0], ninc=n, time=t,
                              T_final=_T0)

    r = sim.solver.solve([step(20., 0.1, 10), step(20., 5., 100), step(0., 0.1, 10),
                          step(0., 5., 100)], umat, props, nstatev, T_init=_T0)
    wd = r["Wm"][3]
    assert wd[-1] > 0.05
    assert np.diff(wd).min() > -1e-6 * wd[-1]
    heat = np.sum(r["r"][1:] * np.diff(r["Time"]))
    np.testing.assert_allclose(heat, wd[-1] - wd[0], rtol=1e-4)


@pytest.mark.parametrize("umat", sorted(_CREEP))
def test_heat_source_derivatives_follow_the_returned_r(umat):
    """drdE/drdT differentiate the r the kernel returns (flow directions frozen, continuum
    dDs/dE): approximate, but no worse than the continuum dSdE they are built on."""
    props = _CREEP[umat][0].copy()
    props[4] = 1e-4                                    # thermal coupling on
    nstatev = _CREEP[umat][1]
    uni = ["strain"] + ["stress"] * 5
    st = StepThermomeca(control=uni, value=[0.01, 0, 0, 0, 0, 0], ninc=20, time=0.2,
                        T_final=_T0 + 5.)
    res = sim.solver.solve([st], umat, props, nstatev, T_init=_T0)
    col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
    etot, sig, sv = col(res["Strain"][:, -1]), col(res["Stress"][:, -1]), col(res["Statev"][:, -1])
    T = np.array([res["Temp"][-1]])
    De0, DT0 = np.array([2e-4, -5e-5, -5e-5, 1e-4, 0., 0.]), 0.5

    def call(De, DT):
        out = sim.umat_T(umat, etot, col(De), sig, np.asfortranarray(np.eye(3)[:, :, None]),
                         col(props), sv.copy(order="F"), 1.0, 0.05, np.zeros((4, 1), order="F"),
                         np.zeros((3, 1), order="F"), T, np.array([DT]))
        return out[0][:, 0], out[5][:, :, 0], out[4].ravel()[0], out[7][:, 0], out[8].ravel()[0]

    _, dSdE, _, drdE, drdT = call(De0, DT0)
    h, hT = 1e-7, 1e-5
    fdS = np.column_stack([(call(De0 + h * e, DT0)[0] - call(De0 - h * e, DT0)[0]) / (2 * h)
                           for e in np.eye(6)])
    fdr = np.array([(call(De0 + h * e, DT0)[2] - call(De0 - h * e, DT0)[2]) / (2 * h)
                    for e in np.eye(6)])
    fdrT = (call(De0, DT0 + hT)[2] - call(De0, DT0 - hT)[2]) / (2 * hT)
    err_S = np.abs(dSdE - fdS).max() / np.abs(fdS).max()
    err_r = np.abs(drdE - fdr).max() / np.abs(fdr).max()
    assert err_r < 1.5 * err_S + 1e-3
    assert abs(drdT - fdrT) < 5e-3 * abs(fdrT)
