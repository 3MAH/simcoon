"""Heat contract of the thermomechanical UMATs (sim.umat_T, the thermomechanical solver).

r is the heat of the increment per unit reference volume, storage included:
r = Delta gamma - T_{n+1} (eta_{n+1} - eta_n), with Delta gamma the dissipation added to Wm_d and
eta = c_0 ln(T/T_init) + alpha:sigma + ... the entropy. r = 0 is the adiabatic condition, at any
time increment (zero included), and drdE / drdT are its derivatives.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import StepThermomeca

RHO, C_P = 7.8e-9, 4.6e8
C0 = RHO * C_P
ALPHA = 1.2e-5
UNI = ["strain"] + ["stress"] * 5
SMA = [5.5, 0.5, 0.5, 0, 67538., 67538., 0.349, 0.349, 1e-6, 1e-6, 0., 0.0418, 0.021, 0., 10., 10.,
       300., 290., 295., 305., 0.2, 0.2, 0.2, 0.2, 300., 0.2, 2., 1e-6, 1e-3, 1., 1e8]

# name: (props, nstatev, strain reached before the probe, T_init, tolerance on drdE (drdT: a tenth))
# Exact for the elastic and viscoelastic kernels; the plastic and SMA ones hold the flow /
# transformation direction fixed over the increment in the derivative of the dissipation.
CASES = {
    "ELISO": ([RHO, C_P, 210000., 0.3, ALPHA], 1, 1e-3, 293.15, 1e-8),
    "ELIST": ([RHO, C_P, 1, 140000., 9000., 0.3, 0.4, 4500., 1e-6, 3e-5], 1, 1e-3, 293.15, 1e-8),
    "ELORT": ([RHO, C_P, 140000., 9000., 9500., 0.3, 0.32, 0.4, 4500., 4600., 3000., 1e-6, 3e-5, 2.5e-5],
              1, 1e-3, 293.15, 1e-8),
    "EPICP": ([RHO, C_P, 210000., 0.3, ALPHA, 300., 500., 0.4], 8, 0.01, 293.15, 1e-3),
    "EPKCP": ([RHO, C_P, 210000., 0.3, ALPHA, 300., 500., 0.4, 20000.], 14, 0.01, 293.15, 1e-3),
    "EPJCK": ([7.85e-9, 4.75e8, 200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0],
              9, 0.02, 293.15, 1e-3),
    "ZENER": ([4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.], 8, 0.01, 293.15, 1e-8),
    "PRONK": ([4.4, 0.656, 3000., 0.35, 1e-4, 1., 1500., 0.35, 3000., 1200.], 14, 0.01, 293.15, 1e-8),
    "SMADI": (SMA, 17, 0.02, 320., 1e-3),
}


def _col(a):
    return np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))


def _loaded_state(name):
    props, nstatev, eps, T0, _ = CASES[name]
    st = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=40, time=0.4, T_final=T0 + 2.)
    res = sim.solver.solve([st], name, np.array(props), nstatev, T_init=T0)
    return (np.array(props), _col(res["Strain"][:, -1]), _col(res["Stress"][:, -1]),
            _col(res["Statev"][:, -1]), np.array([res["Temp"][-1]]))


def _heat(name, state, De, DT, dtime=0.05):
    props, etot, sigma, statev, T = state
    out = sim.umat_T(name, etot, _col(De), sigma, np.asfortranarray(np.eye(3)[:, :, None]), _col(props),
                     statev.copy(order="F"), 1.0, dtime, np.zeros((4, 1), order="F"),
                     np.zeros((3, 1), order="F"), T, np.array([DT]), n_threads=1)
    return out[4].ravel()[0], out[7][:, 0], out[8].ravel()[0]


@pytest.mark.parametrize("name", list(CASES))
def test_heat_derivatives_match_finite_differences(name):
    state = _loaded_state(name)
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    _, drdE, drdT = _heat(name, state, De0, DT0)
    h, hT = 1e-7, 1e-5
    fdE = np.array([(_heat(name, state, De0 + h * e, DT0)[0] - _heat(name, state, De0 - h * e, DT0)[0]) / (2 * h)
                    for e in np.eye(6)])
    fdT = (_heat(name, state, De0, DT0 + hT)[0] - _heat(name, state, De0, DT0 - hT)[0]) / (2 * hT)
    assert np.abs(drdE - fdE).max() < CASES[name][4] * np.abs(fdE).max()
    assert abs(drdT - fdT) < max(1e-6, 0.1 * CASES[name][4]) * abs(fdT)


@pytest.mark.parametrize("name", ["ELISO", "EPICP", "SMADI"])
def test_heat_of_an_increment_is_defined_at_zero_time(name):
    """Rate-independent kernels: the heat of an increment does not depend on its duration, and
    is not zero at DTime = 0 (the instantaneous adiabatic jump is resolvable)."""
    state = _loaded_state(name)
    De, DT = np.array([2e-4, -6e-5, -6e-5, 0., 0., 0.]), 0.0
    r0 = _heat(name, state, De, DT, dtime=0.0)[0]
    assert abs(r0) > 1e-6
    np.testing.assert_allclose(r0, _heat(name, state, De, DT, dtime=0.05)[0], rtol=1e-12)


def test_adiabatic_elastic_cycle_follows_the_isentrope():
    """ELISO under Q = 0: c_0 ln(T/T_init) + alpha tr(sigma) = 0 at every increment (the isentrope,
    exact for the discrete balance), and the cycle returns to T_init."""
    props = np.array([RHO, C_P, 210000., 0.3, ALPHA])
    T0 = 293.15
    load = StepThermomeca(control=UNI, value=[2e-3, 0, 0, 0, 0, 0], ninc=20, time=1.0,
                          thermal_control="heat_flux", Q=0.0)
    unload = StepThermomeca(control=UNI, value=[0., 0, 0, 0, 0, 0], ninc=20, time=1.0,
                            thermal_control="heat_flux", Q=0.0)
    res = sim.solver.solve([load, unload], "ELISO", props, 1, T_init=T0)
    T = np.asarray(res["Temp"]).ravel()
    trs = np.asarray(res["Stress"])[:3].sum(axis=0)
    assert T.min() < T0 - 0.1                       # uniaxial tension cools
    np.testing.assert_allclose(T, T0 * np.exp(-ALPHA * trs / C0), rtol=0, atol=1e-5)
    np.testing.assert_allclose(T[-1], T0, rtol=0, atol=1e-5)


def test_adiabatic_plastic_heating_balances_the_dissipation():
    """EPICP under Q = 0: each increment satisfies T_{n+1} (eta_{n+1} - eta_n) = Delta Wm_d, with
    eta = c_0 ln(T/T_init) + alpha tr(sigma)."""
    props = np.array([RHO, C_P, 210000., 0.3, ALPHA, 300., 500., 0.4])
    T0 = 293.15
    st = StepThermomeca(control=UNI, value=[0.05, 0, 0, 0, 0, 0], ninc=50, time=1.0,
                        thermal_control="heat_flux", Q=0.0)
    res = sim.solver.solve([st], "EPICP", props, 8, T_init=T0)
    T = np.asarray(res["Temp"]).ravel()
    eta = C0 * np.log(T / T0) + ALPHA * np.asarray(res["Stress"])[:3].sum(axis=0)
    dWd = np.diff(np.asarray(res["Wm"])[3])
    assert T[-1] > T0 + 1.0
    np.testing.assert_allclose(T[1:] * np.diff(eta), dWd, rtol=0, atol=1e-5 * dWd.max())
    np.testing.assert_allclose(np.asarray(res["r"]).ravel()[1:], 0.0, atol=1e-6)
