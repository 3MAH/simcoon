"""Thermomechanical blocks under finite strain (control types 2-6), rate form (log_R) or total form.

The kernels are secant in the strain they are fed, so each reads as a law in the log strain: the
entropy eta = c_0 ln(T/T_init) + alpha tr(tau) + ..., the heat r = Delta gamma - T_{n+1} Delta eta
and the works per reference volume, tau the Kirchhoff stress.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, StepThermomeca, solve

RHO, C_P, E, NU, ALPHA = 7.8e-9, 4.6e8, 210000., 0.3, 1.2e-5
C0 = RHO * C_P
T0 = 293.15
ELISO = np.array([RHO, C_P, E, NU, ALPHA])
EPICP = np.array([RHO, C_P, E, NU, ALPHA, 300., 500., 0.4])
UNI = ["strain"] + ["stress"] * 5
# (control type, strain target of the 11 component, kinematic forms it serves)
CONTROLS = [("green_lagrange", 0.25, ("rate", "total")), ("logarithmic", 0.2, ("rate",)),
            ("biot", 1.2, ("rate", "total"))]
CASES = [(ct, eps, form) for ct, eps, forms in CONTROLS for form in forms]


def _tau(res):
    return np.asarray(res["Kirchhoff"])


def _adiabatic(umat, props, nstatev, ct, eps, form, ninc=20):
    st = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=ninc, thermal_control="heat_flux", Q=0.)
    return solve(Block(steps=[st], control_type=ct), umat, props, nstatev, T_init=T0, kinematic_form=form)


@pytest.mark.parametrize("ct, eps, form", CASES)
def test_adiabatic_elastic_stretch_follows_the_isentrope(ct, eps, form):
    """ELISO under Q = 0 and uniaxial stress: c_0 ln(T/T_init) + alpha tr(tau) = 0 at every
    increment, the lateral faces stress-free."""
    res = _adiabatic("ELISO", ELISO, 1, ct, eps, form)
    T = np.asarray(res["Temp"]).ravel()
    tau = _tau(res)
    assert T[-1] < T0 - 5.
    # to the Newton tolerance on the heat balance
    np.testing.assert_allclose(T, T0 * np.exp(-ALPHA * tau[:3].sum(axis=0) / C0), rtol=1e-8)
    np.testing.assert_allclose(tau[1:], 0., atol=1e-6 * np.abs(tau[0]).max())


def test_logarithmic_control_is_the_small_strain_one_in_ln_V():
    """Uniaxial ln V = 0.2 under log_R: the isotropic law is linear in ln V, so tau and T are those
    of the small-strain block driven to eps = 0.2."""
    fin = _adiabatic("ELISO", ELISO, 1, "logarithmic", 0.2, "rate")
    small = _adiabatic("ELISO", ELISO, 1, "small_strain", 0.2, "rate")
    np.testing.assert_allclose(_tau(fin), _tau(small), rtol=1e-9, atol=1e-9 * np.abs(_tau(small)).max())
    np.testing.assert_allclose(np.asarray(fin["Temp"]), np.asarray(small["Temp"]), rtol=1e-12)


@pytest.mark.parametrize("ct, eps, form", CASES)
def test_adiabatic_plastic_heating_balances_the_dissipation(ct, eps, form):
    """EPICP under Q = 0: T_{n+1} (eta_{n+1} - eta_n) = Delta Wm_d at every increment, and the
    heat of every increment vanishes."""
    res = _adiabatic("EPICP", EPICP, 8, ct, eps, form, ninc=40)
    T = np.asarray(res["Temp"]).ravel()
    eta = C0 * np.log(T / T0) + ALPHA * _tau(res)[:3].sum(axis=0)
    dWd = np.diff(np.asarray(res["Wm"])[3])
    assert T[-1] > T0 + 10.
    np.testing.assert_allclose(T[1:] * np.diff(eta), dWd, rtol=0, atol=1e-6 * dWd.max())
    np.testing.assert_allclose(np.asarray(res["r"]).ravel(), 0., atol=1e-6 * dWd.max())


@pytest.mark.parametrize("umat, props, nstatev", [("ELISO", ELISO, 1), ("EPICP", EPICP, 8)])
@pytest.mark.parametrize("ct, eps, form", CASES)
def test_isothermal_block_is_the_mechanical_one(umat, props, nstatev, ct, eps, form):
    """At imposed constant temperature the thermomechanical block is the mechanical one: same
    Kirchhoff stress and works along the path."""
    th = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=20, T_final=T0)
    me = StepMeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=20)
    r_T = solve(Block(steps=[th], control_type=ct), umat, props, nstatev, T_init=T0, kinematic_form=form)
    r_M = solve(Block(steps=[me], control_type=ct), umat, props[2:], nstatev, T_init=T0, kinematic_form=form)
    scale = np.abs(_tau(r_M)).max()
    np.testing.assert_allclose(_tau(r_T), _tau(r_M), rtol=0, atol=1e-8 * scale)
    np.testing.assert_allclose(np.asarray(r_T["Wm"]), np.asarray(r_M["Wm"]), rtol=0,
                               atol=1e-8 * np.abs(np.asarray(r_M["Wm"])).max())


def _shear_F(gamma):
    return np.array([[1., gamma, 0.], [0., 1., 0.], [0., 0., 1.]]).ravel()


@pytest.mark.parametrize("form", ["rate", "total"])
def test_adiabatic_simple_shear_keeps_the_temperature(form):
    """F prescribed (simple shear, isochoric), heat flux Q = 0: tr(ln V) = ln J = 0, so the
    thermoelastic heat vanishes and T stays at T_init, while the shear stress builds up."""
    st = StepThermomeca(control="F", value=_shear_F(1.0), ninc=40, thermal_control="heat_flux", Q=0.)
    res = solve(Block(steps=[st], control_type="F"), "ELISO", ELISO, 1, T_init=T0, kinematic_form=form)
    np.testing.assert_allclose(np.asarray(res["Temp"]), T0, rtol=0, atol=1e-8)
    assert np.abs(_tau(res)[3]).max() > 1e4


@pytest.mark.parametrize("form", ["rate", "total"])
def test_F_control_with_convection_balances_the_heat(form):
    """F prescribed (stretch with shear), 0D convection: every increment satisfies the discrete
    balance r = DTime q_conv (T_{n+1} - T_init) the Newton on DT solves: the heat produced leaves
    through the convection."""
    F = np.array([[1.15, 0.3, 0.], [0., 0.95, 0.], [0., 0., 1.02]]).ravel()
    q_conv = 0.5
    st = StepThermomeca(control="F", value=F, ninc=30, time=3.0, thermal_control="convection", q_conv=q_conv)
    res = solve(Block(steps=[st], control_type="F"), "EPICP", EPICP, 8, T_init=T0, kinematic_form=form)
    T, r = np.asarray(res["Temp"]).ravel(), np.asarray(res["r"]).ravel()
    assert np.abs(T - T0).max() > 1.
    np.testing.assert_allclose(r, 0.1 * q_conv * (T - T0), rtol=0, atol=1e-6)   # Newton tolerance


def test_F_control_with_imposed_temperature_is_the_mechanical_block():
    F = np.array([[1.15, 0.3, 0.], [0., 0.95, 0.], [0., 0., 1.02]]).ravel()
    th = StepThermomeca(control="F", value=F, ninc=30, T_final=T0)
    me = StepMeca(control="F", value=F, ninc=30)
    r_T = solve(Block(steps=[th], control_type="F"), "EPICP", EPICP, 8, T_init=T0)
    r_M = solve(Block(steps=[me], control_type="F"), "EPICP", EPICP[2:], 8, T_init=T0)
    np.testing.assert_allclose(_tau(r_T), _tau(r_M), rtol=0, atol=1e-10 * np.abs(_tau(r_M)).max())


@pytest.mark.parametrize("form", ["jaumann", "green_naghdi"])
def test_thermomechanical_finite_strain_refuses_other_rates(form):
    st = StepThermomeca(control=UNI, value=[0.1, 0, 0, 0, 0, 0], ninc=5, T_final=T0)
    with pytest.raises(RuntimeError, match="corate_type 3"):
        solve(Block(steps=[st], control_type="green_lagrange"), "ELISO", ELISO, 1, T_init=T0, kinematic_form=form)


SMA = [5.5, 0.5, 0.5, 0, 67538., 67538., 0.349, 0.349, 1e-6, 1e-6, 0., 0.0418, 0.021, 0., 10., 10.,
       300., 290., 295., 305., 0.2, 0.2, 0.2, 0.2, 300., 0.2, 2., 1e-6, 1e-3, 1., 1e8]
# name: (props, nstatev, T_init, uniaxial Green-Lagrange target), as test_thermo_heat_contract
LAWS = {
    "ELISO": (ELISO, 1, T0, 0.1),
    "ELIST": ([RHO, C_P, 1, 140000., 9000., 0.3, 0.4, 4500., 1e-6, 3e-5], 1, T0, 0.02),
    "ELORT": ([RHO, C_P, 140000., 9000., 9500., 0.3, 0.32, 0.4, 4500., 4600., 3000., 1e-6, 3e-5, 2.5e-5], 1, T0, 0.02),
    "EPICP": (EPICP, 8, T0, 0.1),
    "EPKCP": ([RHO, C_P, E, NU, ALPHA, 300., 500., 0.4, 20000.], 14, T0, 0.1),
    "EPJCK": ([7.85e-9, 4.75e8, 200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0], 9, T0, 0.1),
    "ZENER": ([4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.], 8, T0, 0.05),
    "PRONK": ([4.4, 0.656, 3000., 0.35, 1e-4, 1., 1500., 0.35, 3000., 1200.], 14, T0, 0.05),
    "SMADI": (SMA, 17, 320., 0.06),
}
SPIN = np.array([[0., -2., 0.], [2., 0., 0.], [0., 0., 0.]])


def _spun(name, form, ninc=40, spin=None):
    props, nstatev, T_init, eps = LAWS[name]
    st = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=ninc, time=1.0,
                        thermal_control="heat_flux", Q=0., BC_w=spin)
    return solve(Block(steps=[st], control_type="green_lagrange"), name, np.array(props), nstatev,
                 T_init=T_init, kinematic_form=form)


@pytest.mark.parametrize("form", ["rate", "total"])
@pytest.mark.parametrize("name", list(LAWS))
def test_every_law_balances_its_heat_in_the_solver(name, form):
    """Adiabatic finite stretch: the heat of the increments sums to Wm_d - Wt (zero here, the
    Newton closing r = 0), and the temperature moves."""
    res = _spun(name, form)
    r, Wm, Wt = np.asarray(res["r"]).ravel(), np.asarray(res["Wm"]), np.asarray(res["Wt"])
    assert abs(res["Temp"][-1] - LAWS[name][2]) > 0.1
    np.testing.assert_allclose(r.sum(), Wm[3, -1] - Wt[0, -1], rtol=0, atol=1e-10 * max(abs(Wt[0, -1]), 1.))


@pytest.mark.parametrize("name", list(LAWS))
def test_total_form_is_objective_under_a_superposed_spin(name):
    """A spin of 2 rad over the step leaves PKII and the temperature unchanged to round-off, for
    every law: the kernels run in the material frame, anisotropy axes included."""
    a, b = _spun(name, "total"), _spun(name, "total", spin=SPIN)
    S = np.asarray(a["PKII"])
    np.testing.assert_allclose(np.asarray(b["PKII"]), S, rtol=0, atol=1e-11 * np.abs(S).max())
    np.testing.assert_allclose(np.asarray(b["Temp"]), np.asarray(a["Temp"]), rtol=0, atol=1e-10)


@pytest.mark.parametrize("name", ["ELIST", "EPICP", "SMADI"])
def test_rate_form_is_objective_to_second_order_under_a_spin(name):
    """Under log_R the discrete increment of the spin is objective asymptotically: the PKII and
    temperature gaps fall as the square of the step, anisotropic law included."""
    def gap(n):
        a, b = _spun(name, "rate", ninc=n), _spun(name, "rate", ninc=n, spin=SPIN)
        S = np.asarray(a["PKII"])
        return (np.abs(np.asarray(b["PKII"])[:, -1] - S[:, -1]).max() / np.abs(S).max(),
                abs(b["Temp"][-1] - a["Temp"][-1]))
    (dS40, dT40), (dS160, dT160) = gap(40), gap(160)
    assert dS40 < 5e-3 and dT40 < 5e-3
    assert dS160 < dS40 / 10. and dT160 < dT40 / 10.


def test_block_start_probe_relaxes_nothing():
    """A second finite block primes its tangent with a zero-time call: PRONK must not relax its
    branches there, so two blocks give the one-block path."""
    props, nstatev, T_init, eps = LAWS["PRONK"]
    s1 = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=20, time=1.0, thermal_control="heat_flux", Q=0.)
    s2 = StepThermomeca(control=UNI, value=[0., 0, 0, 0, 0, 0], ninc=20, time=1.0, thermal_control="heat_flux", Q=0.)
    one = solve(Block(steps=[s1, s2], control_type="green_lagrange"), "PRONK", np.array(props), nstatev, T_init=T_init)
    two = solve([Block(steps=[s1], control_type="green_lagrange"), Block(steps=[s2], control_type="green_lagrange")],
                "PRONK", np.array(props), nstatev, T_init=T_init)
    np.testing.assert_allclose(np.asarray(two["PKII"]), np.asarray(one["PKII"]), rtol=0,
                               atol=1e-12 * np.abs(np.asarray(one["PKII"])).max())
    np.testing.assert_allclose(np.asarray(two["Temp"]), np.asarray(one["Temp"]), rtol=0, atol=1e-10)
