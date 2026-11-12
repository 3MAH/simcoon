"""The thermomechanical face of MODUL (umat_modular_T, dispatch id 200): props [rho, c_p, stream].

The heat r = Delta gamma - T_{n+1} (eta_{n+1} - eta_n) and its derivatives are written once on the
mechanisms (ModularUMAT::run_T, PlasticityMechanism::dissipation_sensitivity), so a composition
equivalent to a dedicated thermomechanical kernel (EPICP_T, EPKCP_T) is that kernel. Every parity
test plays the legacy name twice: as served by default (the modular engine through the legacy
adapter) and as the dedicated kernel (the _closed name).
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import (ArmstrongFrederickHardening, Damage, IsotropicElasticity, ModularMaterial,
                             NeoHookeanElasticity, Plasticity, PowerLawHardening, PragerHardening,
                             Viscoelasticity, VoceHardening, VonMisesYield)
from simcoon.solver import Block, StepMeca, StepThermomeca, solve

E, NU, ALPHA, SY, K, M, KX = 210000., 0.3, 1.2e-5, 300., 500., 0.4, 20000.
RHO, C_P, T0 = 7.8e-9, 4.6e8, 293.15
UNI = ["strain"] + ["stress"] * 5

ISO = ModularMaterial(IsotropicElasticity(E, NU, ALPHA), [Plasticity(SY, VonMisesYield(), PowerLawHardening(K, M))])
PRAGER = ModularMaterial(IsotropicElasticity(E, NU, ALPHA),
                         [Plasticity(SY, VonMisesYield(), PowerLawHardening(K, M), PragerHardening(1.5 * KX))])
AF = ModularMaterial(IsotropicElasticity(E, NU, ALPHA),
                     [Plasticity(SY, VonMisesYield(), VoceHardening(200., 20.), ArmstrongFrederickHardening(30000., 50.))])
EPICP = np.array([RHO, C_P, E, NU, ALPHA, SY, K, M])
EPKCP = np.array([RHO, C_P, E, NU, ALPHA, SY, K, M, KX])


def _props(mat):
    return np.r_[RHO, C_P, mat.props]


def _run(name, props, nstatev, thermal, ct="small_strain", form=None, eps=0.05, ninc=40, tangent_mode=2,
         integrator="cutting_plane", ncp="fischer_burmeister", **step_kw):
    kw = dict(thermal_control="heat_flux", Q=0.) if thermal == "adiabatic" else dict(T_final=T0 + 20.)
    st = StepThermomeca(control=UNI, value=[eps, 0, 0, 0, 0, 0], ninc=ninc, **kw, **step_kw)
    extra = dict(kinematic_form=form) if form else {}
    return solve(Block(steps=[st], control_type=ct), name, np.asarray(props, float), nstatev, T_init=T0,
                 tangent_mode=tangent_mode, integrator=integrator, ncp=ncp, **extra)


def _same(a, b, keys=("Stress", "Temp", "Wm", "Wt"), tol=1e-10):
    for k in keys:
        x, y = np.asarray(a[k]), np.asarray(b[k])
        np.testing.assert_allclose(x, y, rtol=0, atol=tol * np.abs(y).max(), err_msg=k)


ENGINES = ["", "_closed"]
LEGACY = [(ISO, "EPICP", EPICP, 8), (PRAGER, "EPKCP", EPKCP, 14)]


@pytest.mark.parametrize("thermal", ["adiabatic", "imposed"])
@pytest.mark.parametrize("tangent_mode", [1, 2])
@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("mat, legacy, props, nstatev", LEGACY)
def test_modular_face_is_the_dedicated_kernel(mat, legacy, props, nstatev, engine, thermal, tangent_mode):
    """Stress, temperature and every work (Wm, Wm_r, Wm_ir, Wm_d; Wt) of the composition are those of
    EPICP_T / EPKCP_T, to round-off, under both tangent modes."""
    a = _run("MODUL", _props(mat), nstatev, thermal, tangent_mode=tangent_mode)
    b = _run(legacy + engine, props, nstatev, thermal, tangent_mode=tangent_mode)
    _same(a, b)
    if thermal == "imposed":
        _same(a, b, keys=("r",))


@pytest.mark.parametrize("engine", ENGINES)
def test_modular_face_closest_point_is_the_dedicated_kernel(engine):
    a = _run("MODUL", _props(ISO), 8, "adiabatic", integrator="closest_point")
    b = _run("EPICP" + engine, EPICP, 8, "adiabatic", integrator="closest_point")
    _same(a, b)


def test_modular_face_reads_ncp_min():
    a = _run("MODUL", _props(ISO), 8, "adiabatic", ncp="min")
    b = _run("MODUL", _props(ISO), 8, "adiabatic")
    _same(a, b, tol=1e-6)   # the same return mapping, another complementarity function


@pytest.mark.parametrize("form", ["rate", "total"])
@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("mat, legacy, props, nstatev", LEGACY)
def test_modular_face_under_finite_strain_is_the_dedicated_kernel(mat, legacy, props, nstatev, engine, form):
    """Under both finite-strain forms the modular face runs in the same box as the dedicated kernels."""
    a = _run("MODUL", _props(mat), nstatev, "adiabatic", ct="green_lagrange", form=form, eps=0.2)
    b = _run(legacy + engine, props, nstatev, "adiabatic", ct="green_lagrange", form=form, eps=0.2)
    _same(a, b)


# --- heat contract: drdE / drdT against finite differences of r (as test_thermo_heat_contract) -----

def _col(a):
    return np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))


def _loaded_state(props, nstatev):
    st = StepThermomeca(control=UNI, value=[0.01, 0, 0, 0, 0, 0], ninc=40, time=0.4, T_final=T0 + 2.)
    res = solve([st], "MODUL", props, nstatev, T_init=T0)
    return (props, _col(res["Strain"][:, -1]), _col(res["Stress"][:, -1]), _col(res["Statev"][:, -1]),
            np.array([res["Temp"][-1]]))


def _call(state, De, DT, dtime=0.05, integrator="cutting_plane"):
    """One sim.umat_T call on a loaded state: (sigma, r, dSdE, dSdT, drdE, drdT)."""
    props, etot, sigma, statev, T = state
    out = sim.umat_T("MODUL", etot, _col(De), sigma, np.asfortranarray(np.eye(3)[:, :, None]), _col(props),
                     statev.copy(order="F"), 1.0, dtime, np.zeros((4, 1), order="F"),
                     np.zeros((3, 1), order="F"), T, np.array([DT]), n_threads=1, integrator=integrator)
    return out[0][:, 0], out[4].ravel()[0], out[5][:, :, 0], out[6][:, 0], out[7][:, 0], out[8].ravel()[0]


def _heat(state, De, DT, dtime=0.05, integrator="cutting_plane"):
    _, r, _, _, drdE, drdT = _call(state, De, DT, dtime, integrator)
    return r, drdE, drdT


# exact for isotropic hardening (the turning of the flow direction is in the derivative); Prager
# holds the direction of its back stress fixed as EPKCP_T (2e-6, tested at 1e-5); Armstrong-Frederick:
# the secant of the explicit cutting-plane iterates (4e-5), the consistent derivative of the
# backward-Euler update under the closest point (3e-7)
@pytest.mark.parametrize("mat, nstatev, tol, integrator", [(ISO, 8, 1e-8, "cutting_plane"), (PRAGER, 14, 1e-5, "cutting_plane"),
                                                           (AF, 14, 1e-4, "cutting_plane"), (AF, 14, 1e-6, "closest_point")])
def test_heat_derivatives_match_finite_differences(mat, nstatev, tol, integrator):
    state = _loaded_state(_props(mat), nstatev)
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    heat = lambda De, DT: _heat(state, De, DT, integrator=integrator)
    _, drdE, drdT = heat(De0, DT0)
    h, hT = 1e-7, 1e-5
    fdE = np.array([(heat(De0 + h * e, DT0)[0] - heat(De0 - h * e, DT0)[0]) / (2 * h) for e in np.eye(6)])
    fdT = (heat(De0, DT0 + hT)[0] - heat(De0, DT0 - hT)[0]) / (2 * hT)
    assert np.abs(drdE - fdE).max() < tol * np.abs(fdE).max()
    assert abs(drdT - fdT) < max(1e-6, 0.1 * tol) * abs(fdT)


def test_adiabatic_plastic_heating_balances_the_dissipation():
    """EPICP-like composition under Q = 0: T_{n+1} Delta eta = Delta Wm_d at every increment, with the
    thermodynamic split (isotropic hardening work in Wm_ir, back-stress energy in Wm_r, not dissipated)."""
    res = _run("MODUL", _props(AF), 14, "adiabatic")
    T, Wm = np.asarray(res["Temp"]).ravel(), np.asarray(res["Wm"])
    eta = RHO * C_P * np.log(T / T0) + ALPHA * np.asarray(res["Stress"])[:3].sum(axis=0)
    assert T[-1] > T0 + 5. and Wm[2, -1] > 0.
    np.testing.assert_allclose(T[1:] * np.diff(eta), np.diff(Wm[3]), rtol=0, atol=1e-6 * np.diff(Wm[3]).max())
    np.testing.assert_allclose(Wm[0], Wm[1:].sum(axis=0), rtol=1e-12)


# --- entropic hyperelastic block (omega = 1) composed with mechanisms: the entropy of the potential on
# the elastic strain H_e = eps - eps_in - alpha dtheta, its derivatives through dH_e = L^-1 (dsigma_eff - ...)

MU_E, KAPPA_E = E / (2 * (1 + NU)), E / (3 * (1 - 2 * NU))


def _entropic(*mechanisms, mu=MU_E, kappa=KAPPA_E, alpha=ALPHA):
    return ModularMaterial(NeoHookeanElasticity(mu=mu, kappa=kappa, alpha=alpha, omega=1.), list(mechanisms))


ENTROPIC_PLAST = _entropic(Plasticity(SY, VonMisesYield(), PowerLawHardening(K, M)))
ENTROPIC_AF = _entropic(Plasticity(SY, VonMisesYield(), VoceHardening(200., 20.), ArmstrongFrederickHardening(30000., 50.)))
ENTROPIC_VISCO = _entropic(Viscoelasticity(((1500., 0.35, 3000., 1200.),)), mu=3000. / 2.7, kappa=3000. / 0.9, alpha=1e-4)
ENTROPIC_DAMAGE = _entropic(Damage(0.05, 5.))


# dS/dT and the heat derivatives against finite differences (h = 1e-6, hT = 1e-4: the round-off of the
# isochoric potential at these moduli is 1e-11 on r, 1e-7 relative on the quotient at h = 1e-7); the thermal
# rows go through M^-1 of the assembly as the tangent, so dS/dT is as exact as dS/dE: 1e-8 (closest point,
# J2 and Armstrong-Frederick: the entropic stress rate is deviatoric, 2e-4 without that mapping), the
# cutting plane with Armstrong-Frederick 5e-5 (its dS/dE 4e-4, as without the entropic block)
@pytest.mark.parametrize("mat, tol, integrator", [(ENTROPIC_PLAST, 1e-7, "cutting_plane"), (ENTROPIC_PLAST, 1e-7, "closest_point"),
                                                  (ENTROPIC_AF, 1e-6, "closest_point"), (ENTROPIC_VISCO, 1e-8, "cutting_plane"),
                                                  (ENTROPIC_DAMAGE, 1e-7, "cutting_plane")])
def test_entropic_block_with_a_mechanism_heat_derivatives_match_finite_differences(mat, tol, integrator):
    state = _loaded_state(_props(mat), mat.nstatev)
    assert np.asarray(state[3]).ravel()[1] > 0.   # the mechanism has moved (p, v_1 or D)
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    call = lambda De, DT: _call(state, De, DT, integrator=integrator)
    _, _, _, dSdT, drdE, drdT = call(De0, DT0)
    h, hT = 1e-6, 1e-4
    fdE = np.array([(call(De0 + h * e, DT0)[1] - call(De0 - h * e, DT0)[1]) / (2 * h) for e in np.eye(6)])
    fdT = (call(De0, DT0 + hT)[1] - call(De0, DT0 - hT)[1]) / (2 * hT)
    fdST = (call(De0, DT0 + hT)[0] - call(De0, DT0 - hT)[0]) / (2 * hT)
    np.testing.assert_allclose(dSdT, fdST, rtol=0, atol=tol * np.abs(fdST).max())
    np.testing.assert_allclose(drdE, fdE, rtol=0, atol=tol * np.abs(fdE).max())
    np.testing.assert_allclose(drdT, fdT, rtol=0, atol=tol * abs(fdT))


def test_entropic_block_with_plasticity_heat_is_the_dissipation_minus_the_entropy():
    """omega = 1 + J2 plasticity, imposed temperature ramp, total form: r_n = dWm_d - T_{n+1} d eta with
    eta = c_0 ln(T/T_0) + alpha tr tau - W_iso(H_e)/T_0, H_e = ln V - EP - alpha (T - T_0) rebuilt from
    the outputs and the stored plastic strain (statev [T_init, p, EP])."""
    res = _run("MODUL", _props(ENTROPIC_PLAST), ENTROPIC_PLAST.nstatev, "imposed", ct="green_lagrange", form="total",
               eps=0.03, ninc=60)
    T, r, Wm = np.asarray(res["Temp"]).ravel(), np.asarray(res["r"]).ravel(), np.asarray(res["Wm"])
    EP = np.asarray(res["Statev"])[2:5]
    h = np.asarray(res["LogStrain"])[:3] - EP - ALPHA * (T - T0)[None, :]
    lam2 = np.exp(2 * h)
    I1b = np.exp(-2 * h.sum(axis=0) / 3) * lam2.sum(axis=0)
    W_iso = 0.5 * MU_E * (I1b - 3)
    eta = RHO * C_P * np.log(T / T0) + ALPHA * np.asarray(res["Kirchhoff"])[:3].sum(axis=0) - W_iso / T0
    assert Wm[3, -1] > 0. and np.abs(EP).max() > 1e-3
    np.testing.assert_allclose(r[1:], np.diff(Wm[3]) - T[1:] * np.diff(eta), rtol=0, atol=1e-8 * np.abs(r).max())


# --- viscoelastic face: PRONK_T / ZENER_T parity, the heat contract, the thermal strain of the branches

E_V, NU_V, ALPHA_V, RHO_V, C_P_V = 3000., 0.35, 1e-4, 4.4, 0.656
BRANCH = (1500., 0.35, 3000., 1200.)   # E_1, nu_1, etaB_1, etaS_1
VISCO = ModularMaterial(IsotropicElasticity(E_V, NU_V, ALPHA_V), [Viscoelasticity((BRANCH,))])
VISCO_PRAGER = ModularMaterial(IsotropicElasticity(E_V, NU_V, ALPHA_V),
                               [Viscoelasticity((BRANCH,)),
                                Plasticity(30., VonMisesYield(), PowerLawHardening(50., 0.4), PragerHardening(300.))])
PRONK = np.array([RHO_V, C_P_V, E_V, NU_V, ALPHA_V, 1., *BRANCH])
ZENER = np.array([RHO_V, C_P_V, E_V, NU_V, ALPHA_V, *BRANCH])


def _props_v(mat):
    return np.r_[RHO_V, C_P_V, mat.props]


def _zener_as_maxwell(E0, nu0, E1, nu1, etaB, etaS):
    """The Kelvin branch (L_1, H_1) in series with L_0 as the Maxwell branch in parallel with it, per
    isotropic mode X in (K, G): X_1^M = X_0^2 / (X_0 + X_1), eta^M = eta^K X_0^2 / (X_0 + X_1)^2 (same
    instantaneous and long-term moduli, same relaxation time)."""
    K0, G0 = E0 / (3 * (1 - 2 * nu0)), E0 / (2 * (1 + nu0))
    K1, G1 = E1 / (3 * (1 - 2 * nu1)), E1 / (2 * (1 + nu1))
    K1m, G1m = K0 ** 2 / (K0 + K1), G0 ** 2 / (G0 + G1)
    return (9 * K1m * G1m / (3 * K1m + G1m), (3 * K1m - 2 * G1m) / (2 * (3 * K1m + G1m)),
            etaB * K0 ** 2 / (K0 + K1) ** 2, etaS * G0 ** 2 / (G0 + G1) ** 2)


@pytest.mark.parametrize("thermal", ["adiabatic", "imposed"])
@pytest.mark.parametrize("tangent_mode", [1, 2])
@pytest.mark.parametrize("engine", ENGINES)
def test_viscoelastic_face_is_pronk_T(engine, thermal, tangent_mode):
    """A Maxwell branch on the elastic block is PRONK_T: stress, temperature, the works and the heat,
    with the thermal strain in the branch driving force (imposed temperature, alpha != 0). Under
    tangent_mode 0 both return the elastic operator and the solver's Newton stops within its own
    tolerance (temperatures 1e-7 K apart), so that mode is not a parity statement."""
    a = _run("MODUL", _props_v(VISCO), 8, thermal, tangent_mode=tangent_mode, time=2.)
    b = _run("PRONK" + engine, PRONK, 14, thermal, tangent_mode=tangent_mode, time=2.)
    _same(a, b)
    if thermal == "imposed":
        _same(a, b, keys=("r",))


@pytest.mark.parametrize("form", ["rate", "total"])
@pytest.mark.parametrize("engine", ENGINES)
def test_viscoelastic_face_under_finite_strain_is_pronk_T(engine, form):
    a = _run("MODUL", _props_v(VISCO), 8, "adiabatic", ct="green_lagrange", form=form, eps=0.2, time=2.)
    b = _run("PRONK" + engine, PRONK, 14, "adiabatic", ct="green_lagrange", form=form, eps=0.2, time=2.)
    _same(a, b)


@pytest.mark.parametrize("thermal", ["adiabatic", "imposed"])
def test_zener_T_is_a_one_branch_maxwell_composition(thermal):
    """ZENER_T (a Kelvin branch in series with the spring) is the Maxwell composition with the mapped
    branch: same state ODE under a linear change of the viscous strain, so the same stress, works,
    entropy and heat (to the round-off of the two inverses)."""
    mat = ModularMaterial(IsotropicElasticity(E_V, NU_V, ALPHA_V), [Viscoelasticity((_zener_as_maxwell(E_V, NU_V, *BRANCH),))])
    a = _run("MODUL", _props_v(mat), 8, thermal, time=2.)
    b = _run("ZENER", ZENER, 8, thermal, time=2.)
    _same(a, b, tol=1e-10)
    if thermal == "imposed":
        _same(a, b, keys=("r",), tol=1e-10)


@pytest.mark.parametrize("mat, nstatev, tol", [(VISCO, 8, 1e-8), (VISCO_PRAGER, 21, 1e-5)])
def test_viscoelastic_heat_derivatives_match_finite_differences(mat, nstatev, tol):
    """drdE / drdT of the viscous dissipation are exact (the branches are linear, nothing held fixed);
    with a Prager mechanism on top the back-stress direction is held fixed as EPKCP_T."""
    state = _loaded_state(_props(mat), nstatev)
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    _, drdE, drdT = _heat(state, De0, DT0)
    h, hT = 1e-7, 1e-5
    fdE = np.array([(_heat(state, De0 + h * e, DT0)[0] - _heat(state, De0 - h * e, DT0)[0]) / (2 * h)
                    for e in np.eye(6)])
    fdT = (_heat(state, De0, DT0 + hT)[0] - _heat(state, De0, DT0 - hT)[0]) / (2 * hT)
    np.testing.assert_allclose(drdE, fdE, rtol=0, atol=tol * np.abs(fdE).max())
    np.testing.assert_allclose(drdT, fdT, rtol=0, atol=tol * abs(fdT))


def test_viscous_heat_over_a_strain_hold_is_the_dissipation():
    """Ramp then hold at constant temperature (as the PRONK_T test): the heat released over the hold
    is the dissipation minus the thermoelastic term T alpha : dsigma; nothing stored (Wm_ir = 0)."""
    ramp = StepThermomeca(control=UNI, value=[0.01, 0, 0, 0, 0, 0], ninc=20, time=0.1, T_final=T0)
    hold = StepThermomeca(control=UNI, value=[0.01, 0, 0, 0, 0, 0], ninc=200, time=5.0, T_final=T0)
    res = solve([ramp, hold], "MODUL", _props_v(VISCO), 8, T_init=T0)
    t, heat_inc, Wm = np.asarray(res["Time"]), np.asarray(res["r"]).ravel(), np.asarray(res["Wm"])
    i1 = int(np.argmin(np.abs(t - 0.1)))
    dWd = Wm[3, -1] - Wm[3, i1]
    dsig = np.asarray(res["Stress"])[:, -1] - np.asarray(res["Stress"])[:, i1]
    assert dWd > 0. and np.abs(Wm[2]).max() == 0.
    np.testing.assert_allclose(np.sum(heat_inc[i1 + 1:]), dWd - T0 * ALPHA_V * np.sum(dsig[:3]), rtol=1e-8)
    np.testing.assert_allclose(Wm[0], Wm[1:].sum(axis=0), rtol=1e-12)


@pytest.mark.parametrize("engine", ENGINES)
def test_mechanical_block_drives_the_branches_with_the_thermal_strain_removed(engine):
    """A mechanical block with a temperature change: the branches see eps - alpha (T - T_init), as
    PRONK; free thermal expansion creates no viscous strain and dissipates nothing."""
    props_m = np.asarray(VISCO.props, float)
    st = StepMeca(control=UNI, value=[0.02, 0, 0, 0, 0, 0], ninc=40, time=1., T_final=T0 + 30.)
    a = solve(Block(steps=[st]), "MODUL", props_m, 8, T_init=T0)
    b = solve(Block(steps=[st]), "PRONK" + engine, PRONK[2:], 14, T_init=T0)
    _same(a, b, keys=("Stress", "Wm"), tol=1e-12)
    free = StepMeca(control=["stress"] * 6, value=[0] * 6, ninc=20, time=1., T_final=T0 + 50.)
    c = solve(Block(steps=[free]), "MODUL", props_m, 8, T_init=T0)
    assert np.abs(np.asarray(c["Wm"])).max() < 1e-14 and np.abs(np.asarray(c["Statev"])[1:7]).max() < 1e-14


# --- damage face: the trapezoid 1/2 (Y_n + Y_n+1) dD on the driving force, its heat and derivatives

Y0, YC = 10., 2000.
DAMAGE_0 = ModularMaterial(IsotropicElasticity(E, NU, 0.), [Damage(Y0, YC)])
DAMAGE = ModularMaterial(IsotropicElasticity(E, NU, ALPHA), [Damage(Y0, YC)])
DAMAGE_PLAST = ModularMaterial(IsotropicElasticity(E, NU, ALPHA),
                               [Plasticity(SY, VonMisesYield(), PowerLawHardening(K, M)), Damage(0.05, 5.)])
DAMAGE_PRAGER = ModularMaterial(IsotropicElasticity(E, NU, ALPHA),
                                [Plasticity(SY, VonMisesYield(), PowerLawHardening(K, M), PragerHardening(1.5 * KX)),
                                 Damage(0.05, 5.)])


def test_damage_dissipation_is_the_trapezoid_on_the_driving_force():
    """Uniaxial, alpha = 0, imposed temperature ramp: the heat of each increment is the dissipation
    minus T_{n+1} c_0 ln(T_{n+1}/T_n); the trapezoid is exact for the linear law (it telescopes to
    (Y_{n+1}^2 - Y_n^2) / (2 (Y_c - Y_0)) once damage grows at both ends); nothing stored."""
    res = _run("MODUL", _props(DAMAGE_0), DAMAGE_0.nstatev, "imposed", eps=0.05, ninc=40)
    Wm, r, S, T = (np.asarray(res[k]) for k in ("Wm", "r", "Stress", "Temp"))
    D = np.asarray(res["Statev"])[1]
    T = T.ravel()
    assert D[-1] > 0.1 and np.abs(Wm[2]).max() == 0.
    np.testing.assert_allclose(r.ravel()[1:], np.diff(Wm[3]) - T[1:] * RHO * C_P * np.log(T[1:] / T[:-1]),
                               rtol=0, atol=1e-10 * np.abs(r).max())
    Y = 0.5 * (S[0] / (1 - D)) ** 2 / E
    on = np.nonzero(D[:-1] > 0.)[0]   # increments whose both ends are on the damage branch
    np.testing.assert_allclose(np.diff(Wm[3])[on].sum(), (Y[on[-1] + 1] ** 2 - Y[on[0]] ** 2) / (2 * (YC - Y0)), rtol=1e-12)
    np.testing.assert_allclose(Wm[0], Wm[1:].sum(axis=0), rtol=1e-12)


@pytest.mark.parametrize("mat, tol", [(DAMAGE, 1e-8), (DAMAGE_PLAST, 1e-8), (DAMAGE_PRAGER, 1e-5)])
def test_damage_heat_derivatives_match_finite_differences(mat, tol):
    """drdE / drdT on the effective stress (D'(Y) and the compliance of the committed state), alone
    and on top of a plasticity mechanism, damage growing in the probed increment; with a Prager
    back stress the turning term is on the effective stress (its direction held fixed, as without
    damage: 2e-6 at this probe)."""
    state = _loaded_state(_props(mat), mat.nstatev)
    assert np.asarray(state[3]).ravel()[-2] > 0.   # D of the loaded state
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    _, drdE, drdT = _heat(state, De0, DT0)
    h, hT = 1e-7, 1e-5
    fdE = np.array([(_heat(state, De0 + h * e, DT0)[0] - _heat(state, De0 - h * e, DT0)[0]) / (2 * h)
                    for e in np.eye(6)])
    fdT = (_heat(state, De0, DT0 + hT)[0] - _heat(state, De0, DT0 - hT)[0]) / (2 * hT)
    np.testing.assert_allclose(drdE, fdE, rtol=0, atol=tol * np.abs(fdE).max())
    np.testing.assert_allclose(drdT, fdT, rtol=0, atol=tol * abs(fdT))


# --- hyperelastic block: omega = 0 is the linear face at small strain, omega = 1 the entropic rubber

MU_R, KAPPA_R, ALPHA_R = 0.5, 1000., 2e-4
RUBBER = ModularMaterial(NeoHookeanElasticity(mu=MU_R, kappa=KAPPA_R, alpha=ALPHA_R, omega=1.))


def test_energetic_hyperelastic_block_is_the_linear_face_at_small_strain():
    """omega = 0, neo-Hookean (mu, kappa) against ELISO_T with the same ground-state moduli, under the
    total form at 1e-4 strain: stress, temperature, works and heat agree to the order of the strain."""
    mu, kappa = E / (2 * (1 + NU)), E / (3 * (1 - 2 * NU))
    hyper = ModularMaterial(NeoHookeanElasticity(mu=mu, kappa=kappa, alpha=ALPHA))
    for thermal in ("adiabatic", "imposed"):
        a = _run("MODUL", _props(hyper), hyper.nstatev, thermal, ct="green_lagrange", form="total", eps=1e-4)
        b = _run("ELISO", np.array([RHO, C_P, E, NU, ALPHA]), 1, thermal, ct="green_lagrange", form="total", eps=1e-4)
        _same(a, b, keys=("Stress", "Wm"), tol=1e-3)
        np.testing.assert_allclose(np.asarray(a["Temp"]), np.asarray(b["Temp"]), rtol=0, atol=1e-6)
        np.testing.assert_allclose(np.asarray(a["Wt"]), np.asarray(b["Wt"]), rtol=0,   # Wt is round-off under Q = 0
                                   atol=1e-6 * max(np.abs(b["Wt"]).max(), np.abs(b["Wm"]).max()))
        if thermal == "imposed":
            _same(a, b, keys=("r",), tol=1e-3)


def test_entropic_rubber_heat_is_the_entropy_of_the_potential():
    """omega = 1, imposed temperature ramp: r_n = -T_{n+1} [c_0 ln(T_{n+1}/T_n) + alpha d(tr tau)
    - d W_iso / theta_0] exactly (no dissipation), W_iso = mu/2 (I1_bar - 3) rebuilt from ln V
    - alpha (T - T_0) I; nothing dissipated nor stored."""
    res = _run("MODUL", _props(RUBBER), RUBBER.nstatev, "imposed", ct="green_lagrange", form="total", eps=0.3, ninc=60)
    T, r, Wm = np.asarray(res["Temp"]).ravel(), np.asarray(res["r"]).ravel(), np.asarray(res["Wm"])
    h = np.asarray(res["LogStrain"])[:3] - ALPHA_R * (T - T0)[None, :]
    lam2 = np.exp(2 * h)
    I1b = np.exp(-2 * h.sum(axis=0) / 3) * lam2.sum(axis=0)
    W_iso = 0.5 * MU_R * (I1b - 3)
    tr_tau = np.asarray(res["Kirchhoff"])[:3].sum(axis=0)
    eta = RHO * C_P * np.log(T / T0) + ALPHA_R * tr_tau - W_iso / T0
    np.testing.assert_allclose(r[1:], -T[1:] * np.diff(eta), rtol=0, atol=1e-8 * np.abs(r).max())
    assert np.abs(Wm[2:]).max() == 0. and np.abs(W_iso).max() > 1e-2


def test_entropic_rubber_heat_derivatives_match_finite_differences():
    state = _loaded_state(_props(RUBBER), RUBBER.nstatev)
    De0, DT0 = np.array([2e-4, -6e-5, -6e-5, 1e-4, 0., 0.]), 0.5
    _, drdE, drdT = _heat(state, De0, DT0)
    h, hT = 1e-7, 1e-5
    fdE = np.array([(_heat(state, De0 + h * e, DT0)[0] - _heat(state, De0 - h * e, DT0)[0]) / (2 * h)
                    for e in np.eye(6)])
    fdT = (_heat(state, De0, DT0 + hT)[0] - _heat(state, De0, DT0 - hT)[0]) / (2 * hT)
    np.testing.assert_allclose(drdE, fdE, rtol=0, atol=1e-8 * np.abs(fdE).max())
    np.testing.assert_allclose(drdT, fdT, rtol=0, atol=1e-8 * abs(fdT))


def test_entropic_rubber_shows_the_thermoelastic_inversion():
    """Adiabatic uniaxial stretch of the entropic rubber: the temperature first drops (thermal
    expansion) then rises (Gough-Joule), the minimum at lambda - 1 = alpha theta_0 (Joule 1859)."""
    res = _run("MODUL", _props(RUBBER), RUBBER.nstatev, "adiabatic", ct="green_lagrange", form="total", eps=0.2, ninc=200)
    T, e = np.asarray(res["Temp"]).ravel(), np.asarray(res["LogStrain"])[0]
    i = int(np.argmin(T))
    assert T[i] < T0 - 1e-4 and T[-1] > T0 + 1e-3, (T[i] - T0, T[-1] - T0)
    np.testing.assert_allclose(e[i], ALPHA_R * T0, rtol=0.2, err_msg=f"inversion at {e[i]}")


# --- closed cycles: what is recovered and what is not

def _cycle(name, props, nstatev, ends, ninc=200, ct="small_strain", form=None, **thermal):
    steps = [StepThermomeca(control=UNI, value=[e, 0, 0, 0, 0, 0], ninc=ninc, **thermal) for e in ends]
    extra = dict(kinematic_form=form) if form else {}
    return solve(Block(steps=steps, control_type=ct), name, np.asarray(props, float), nstatev, T_init=T0,
                 tangent_mode=2, **extra)


@pytest.mark.parametrize("name, props, nstatev", [("MODUL", _props(PRAGER), PRAGER.nstatev), ("EPKCP", EPKCP, 14),
                                                  ("EPKCP_closed", EPKCP, 14)])
def test_back_stress_energy_is_recoverable(name, props, nstatev):
    """Tension then compression with a Prager back stress X = kX a, isothermal: Wm_r is the elastic energy
    plus the energy of the back stress, sigma^2 / 2E + 3/4 kX ep^2, a state function given back when the
    plastic strain returns to zero; Wm_ir only holds the isotropic hardening work and never decreases."""
    res = _cycle(name, props, nstatev, (0.02, -0.02), T_final=T0)
    s, e, Wm = np.asarray(res["Stress"])[0], np.asarray(res["Strain"])[0], np.asarray(res["Wm"])
    ep = e - s / E
    W_el, W_kin = s**2 / (2 * E), 0.75 * KX * ep**2
    np.testing.assert_allclose(Wm[1], W_el + W_kin, rtol=0, atol=1e-9 * Wm[0].max())
    np.testing.assert_allclose(Wm[0], Wm[1:].sum(axis=0), rtol=1e-12, atol=1e-12)
    i_peak = len(s) // 2
    assert W_kin[i_peak] > 1. and W_kin[i_peak:].min() < 1e-3 * W_kin[i_peak]   # stored, then recovered
    assert np.diff(Wm[2]).min() >= 0. and Wm[2, -1] > 0. and np.diff(Wm[3]).min() >= -1e-12


@pytest.mark.parametrize("thermal", [dict(T_final=T0), dict(thermal_control="heat_flux", Q=0.)])
def test_entropic_rubber_cycle_closes_without_hysteresis(thermal):
    """omega = 1, stretch and release: no dissipation (Wm_ir = Wm_d = 0), so the cycle closes under
    isothermal and adiabatic conditions, the work and the heat of the loading given back on unloading."""
    res = _cycle("MODUL", _props(RUBBER), RUBBER.nstatev, (0.5, 0.), ct="green_lagrange", form="total", **thermal)
    T, Wm, r = np.asarray(res["Temp"]).ravel(), np.asarray(res["Wm"]), np.asarray(res["r"]).ravel()
    assert np.abs(Wm[2:]).max() == 0. and Wm[0].max() > 0.1
    assert abs(Wm[0, -1]) < 1e-9 * Wm[0].max() and abs(T[-1] - T0) < 1e-6
    assert abs(r.sum()) < 1e-5 * max(np.abs(r).sum(), Wm[0].max())


def test_entropic_rubber_hysteresis_is_the_heat_exchanged():
    """The same cycle with a convective exchange: the loop that opens is thermoelastic damping, not a
    dissipation of the law (Wm_d = 0); it is booked in Wm_r and is the net heat given to the surroundings."""
    res = _cycle("MODUL", _props(RUBBER), RUBBER.nstatev, (0.5, 0.), ct="green_lagrange", form="total",
                 thermal_control="convection", q_conv=10.)
    Wm = np.asarray(res["Wm"])
    assert np.abs(Wm[2:]).max() == 0. and 1e-7 < Wm[0, -1] < 1e-3 * Wm[0].max()
    assert Wm[1, -1] == Wm[0, -1]
