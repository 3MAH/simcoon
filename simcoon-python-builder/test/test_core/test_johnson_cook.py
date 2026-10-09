"""EPJCK (Johnson-Cook) through the Python API: the block solver at several strain rates,
the thermomechanical twin under adiabatic self-heating, and the point-wise binding in the
tangent modes (0 = elastic L for explicit integration).

The C++ kernel tests (Ttangent_EPJCK) check the law against its closed form; here the checks
are the ones a user of ``sim.solver.solve`` relies on.
"""

import functools

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import StepMeca, StepThermomeca, solve

# AISI 4340 (Johnson & Cook, 1983): E, nu, alpha, A, B, n, C, edot0, m, T_ref, T_melt
JC = np.array([200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0])
NSTATEV = 9
UNIAXIAL = ["strain"] + ["stress"] * 5


@functools.lru_cache(maxsize=None)
def _tension(time, eps=0.10, ninc=100, T_init=293.15):
    # cached: several tests read the same uniaxial runs; the results are not mutated
    step = StepMeca(control=UNIAXIAL, value=[eps, 0, 0, 0, 0, 0], time=time, ninc=ninc)
    return solve(step, "EPJCK", JC, NSTATEV, T_init=T_init)


def _power_hardening(p, B, n, p_reg=1.0e-6):
    """B p^n, exact for p >= p_reg and the C1 quadratic onset blend of PowerLawHardening below."""
    if n >= 1.0 or p >= p_reg:
        return B * p**n if p > 0.0 else 0.0
    a = B * p_reg ** (n - 1.0) * (2.0 - n)
    b = B * p_reg ** (n - 2.0) * (n - 1.0)
    return a * p + b * p * p if p > 0.0 else a * p


def _jc_yield(p, pdot, T):
    A, B, n, C, edot0, m, T_ref, T_melt = JC[3:11]
    rate = 1.0 + C * max(np.log(pdot / edot0), 0.0) if pdot > 0.0 else 1.0
    Tstar = min(max((T - T_ref) / (T_melt - T_ref), 0.0), 1.0)
    return (A + _power_hardening(p, B, n)) * rate * (1.0 - Tstar**m)


@pytest.mark.parametrize("time, pdot_expected", [(10.0, 0.01), (1.0, 0.1), (1.0e-3, 100.0)])
def test_uniaxial_tension_sits_on_the_johnson_cook_surface(time, pdot_expected):
    res = _tension(time)
    s11, p, pdot = res["Stress"][0], res["Statev"][1], res["Statev"][8]
    # uniaxial: Mises = s11; the plastic strain rate is the one of the loading
    assert pdot[-1] == pytest.approx(pdot_expected, rel=0.05)
    # on the surface at EVERY plastic point, the onset increment included (the un-regularized
    # power hardening used to leave it 2-3 % below the surface)
    plastic = p > 0.0
    sigma_Y = np.array([_jc_yield(pk, pdk, 293.15) for pk, pdk in zip(p[plastic], pdot[plastic])])
    np.testing.assert_allclose(s11[plastic], sigma_Y, rtol=1e-6)


def test_fine_time_steps_stay_on_the_surface_and_refine_consistently():
    """Below DTime ~ C (A + B p^n) / (3 G edot0) (~5e-5 s for 4340) the clamped rate factor used
    to make the cutting-plane Newton cycle between its two branches: silently off-surface states
    and a final stress 30 % too high at 6400 increments. The branch selection at the kink keeps
    every increment on the surface and the response independent of the step."""
    coarse = _tension(1.0e-3, ninc=100)
    fine = _tension(1.0e-3, ninc=3200)          # DTime = 3e-7 s
    s11, p, pdot = fine["Stress"][0], fine["Statev"][1], fine["Statev"][8]
    plastic = p > 0.0
    sigma_Y = np.array([_jc_yield(pk, pdk, 293.15) for pk, pdk in zip(p[plastic], pdot[plastic])])
    np.testing.assert_allclose(s11[plastic], sigma_Y, rtol=1e-7)
    assert fine["Stress"][0, -1] == pytest.approx(coarse["Stress"][0, -1], rel=1e-4)


def test_stress_increases_with_strain_rate_and_is_clamped_below_edot0():
    slow, ref, fast = (_tension(t)["Stress"][0, -1] for t in (10.0, 0.1, 1.0e-3))
    # below the reference rate the rate factor is 1: slow and reference-rate runs coincide
    assert slow == pytest.approx(ref, rel=1e-8)
    assert fast > 1.05 * ref


def test_thermal_softening_lowers_the_flow_stress():
    cold = _tension(1.0, T_init=293.15)["Stress"][0, -1]
    hot = _tension(1.0, T_init=793.15)["Stress"][0, -1]
    assert hot < 0.75 * cold


def test_thermomechanical_twin_heats_up_adiabatically_and_softens():
    rho, c_p = 7.85e-9, 4.75e8   # t/mm^3, mJ/(t K)
    props_T = np.concatenate([[rho, c_p], JC])
    step = StepThermomeca(control=UNIAXIAL, value=[0.10, 0, 0, 0, 0, 0], time=1.0e-3, ninc=100,
                          thermal_control="heat_flux", Q=0.0)
    res = solve(step, "EPJCK", props_T, NSTATEV, T_init=293.15)
    T_end, Wm_d, Wt = res["Temp"][-1], res["Wm"][3, -1], res["Wt"][0, -1]
    assert T_end > 293.15 + 10.0, "plastic dissipation must heat the material"
    # adiabatic first law: the thermal work (heat capacity + thermoelastic + hardening-entropy
    # couplings) equals the dissipated work, to the step error
    assert Wt == pytest.approx(Wm_d, rel=0.02)
    assert rho * c_p * (T_end - 293.15) == pytest.approx(Wm_d, rel=0.15)
    # the isothermal mechanical run at the same rate is stiffer than the self-heated one
    iso = _tension(1.0e-3)["Stress"][0, -1]
    assert res["Stress"][0, -1] < iso


def _point(tangent_mode, DTime=1.0e-3):
    De = np.array([[0.01, -0.005, -0.005, 0.0, 0.0, 0.0]]).T
    col = lambda v: np.asfortranarray(np.asarray(v, dtype=float).reshape(-1, 1))
    DR = np.asfortranarray(np.eye(3)[:, :, None])
    return sim.umat("EPJCK", col(np.zeros(6)), col(De), np.empty(0), np.empty(0), col(np.zeros(6)), DR,
                    col(JC), col(np.zeros(NSTATEV)), 0.0, DTime, col(np.zeros(4)),
                    temp=np.array([293.15]), n_threads=1, tangent_mode=tangent_mode, start=True)


def test_tangent_mode_none_returns_elastic_operator_for_explicit_integration():
    sigma0, statev0, _, Lt0 = _point(0)
    sigma2, _, _, Lt2 = _point(2)
    assert statev0[1, 0] > 0.0, "the point must yield"
    L = sim.L_iso([JC[0], JC[1]], "Enu")
    np.testing.assert_allclose(Lt0[:, :, 0], L, rtol=1e-12)
    np.testing.assert_allclose(sigma0, sigma2, rtol=1e-12)
    assert not np.allclose(Lt2[:, :, 0], L, rtol=1e-3)
