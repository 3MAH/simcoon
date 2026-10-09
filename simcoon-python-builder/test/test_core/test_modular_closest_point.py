"""tangent_mode = 3 (closest-point projection) on the modular UMAT.

The cutting-plane loop (modes 1/2) accumulates the plastic strain along the flow direction of
each iterate; the closest-point branch solves stress, multipliers and backward-Euler state
together, so the consistent operator becomes the exact Jacobian of the discrete update for
every criterion with a flow Hessian and every hardening law.
"""
import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import (ModularMaterial, IsotropicElasticity, NeoHookeanElasticity, Plasticity,
                             VonMisesYield, HillYield, TrescaYield, VoceHardening,
                             LinearIsotropicHardening, PragerHardening, ArmstrongFrederickHardening,
                             ChabocheHardening, Viscoelasticity, Damage)
from simcoon.solver import StepMeca, solve

_UNIAXIAL = ["strain"] + ["stress"] * 5
HILL = dict(F=0.5, G=0.7, H=0.4, L=1.6, M=1.3, N=1.8)
col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
I3 = np.eye(3).reshape(3, 3, 1).copy(order="F")

CASES = {
    # name: (plasticity, exact operator symmetric?)
    "J2+Voce": (Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                           isotropic_hardening=VoceHardening(Q=10., b=50.)), True),
    "J2+Prager": (Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                             isotropic_hardening=LinearIsotropicHardening(H=100.),
                             kinematic_hardening=PragerHardening(C=2000.)), True),
    "J2+AF": (Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                         isotropic_hardening=LinearIsotropicHardening(H=100.),
                         kinematic_hardening=ArmstrongFrederickHardening(C=2000., D=50.)), False),
    "J2+Chaboche2": (Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                                isotropic_hardening=VoceHardening(Q=10., b=50.),
                                kinematic_hardening=ChabocheHardening(terms=((2000., 50.), (500., 5.)))),
                     False),
    "Hill+Voce": (Plasticity(sigma_Y=20., yield_criterion=HillYield(**HILL),
                             isotropic_hardening=VoceHardening(Q=10., b=50.)), True),
    "Hill+AF": (Plasticity(sigma_Y=20., yield_criterion=HillYield(**HILL),
                           isotropic_hardening=LinearIsotropicHardening(H=100.),
                           kinematic_hardening=ArmstrongFrederickHardening(C=2000., D=50.)), False),
}


def _material(*mechs):
    return ModularMaterial(elasticity=IsotropicElasticity(C1=3000., C2=0.35), mechanisms=list(mechs))


def _plastic_state(mat, mode, ninc=10):
    r = solve(StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=ninc, time=0.5),
              "MODUL", mat.props, mat.nstatev, T_init=290., tangent_mode=mode)
    assert r.status == 0
    return r


def _umat(mat, e, s, sv, De, mode, time=0.5, dtime=0.3):
    """One increment of the modular UMAT from state (e, s, sv): (sigma, statev, Wm, Lt)."""
    return sim.umat("MODUL", col(e), col(De), I3, I3, col(s), I3, col(mat.props), col(sv),
                    time, dtime, np.zeros((4, 1), order="F"), n_threads=1, tangent_mode=mode)


def _tangent_and_fd(mat, e, s, sv, mode, De0=None, h=1e-8):
    De0 = np.array([4e-4, 1e-4, -1e-4, 1.5e-4, 0., 0.]) if De0 is None else De0
    Lt = _umat(mat, e, s, sv, De0, mode)[3][:, :, 0]
    fd = np.column_stack([(_umat(mat, e, s, sv, De0 + h * u, mode)[0][:, 0]
                           - _umat(mat, e, s, sv, De0 - h * u, mode)[0][:, 0]) / (2 * h)
                          for u in np.eye(6)])
    return Lt, fd


@pytest.mark.parametrize("name", list(CASES))
def test_closest_point_tangent_is_the_exact_jacobian(name):
    """Mode 3: Lt matches central differences of the discrete map to FD noise (1e-9) for
    kinematic hardening and anisotropic criteria (the mode-2 operator is 1e-4..1e-3 off there,
    recorded in the PR). The exact operator is symmetric for isotropic / Prager hardening and
    NOT for Armstrong-Frederick / Chaboche (dynamic recovery is not generalised-standard)."""
    pl, symmetric = CASES[name]
    mat = _material(pl)
    r = _plastic_state(mat, 3)
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    Lt3, fd3 = _tangent_and_fd(mat, e, s, sv, 3)
    err3 = np.linalg.norm(Lt3 - fd3) / np.linalg.norm(fd3)
    assert err3 < 1e-8, err3
    sym = np.linalg.norm(Lt3 - Lt3.T) / np.linalg.norm(Lt3)
    if symmetric:
        assert sym < 1e-12, sym
    else:
        assert sym > 1e-7, sym   # recorded asymmetry of the exact AF/Chaboche operator


def test_closest_point_equals_cutting_plane_for_radial_return():
    """J2 + isotropic hardening: the normal does not rotate, the two integrators solve the
    same equations; each stops at the local tolerance (|Phi|/sigma_Y < 1e-9), so they agree to
    that order, not to the bit."""
    mat = _material(CASES["J2+Voce"][0])
    r2, r3 = _plastic_state(mat, 2), _plastic_state(mat, 3)
    np.testing.assert_allclose(r3["Stress"], r2["Stress"], rtol=1e-7, atol=1e-7)
    np.testing.assert_allclose(r3["Statev"], r2["Statev"], rtol=1e-7, atol=1e-10)


def test_closest_point_differs_from_cutting_plane_by_second_order():
    """Where the normal rotates (Hill + AF) the two integrators differ by O(h^2) on ONE increment
    from the same state: halving the increment divides the gap by about 4 (over a whole path
    both are first order, and their global gap only halves)."""
    mat = _material(CASES["Hill+AF"][0])
    r = _plastic_state(mat, 3)
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    De0 = np.array([4e-3, 1e-3, -1e-3, 1.5e-3, 0., 0.])

    stress = lambda d, mode: _umat(mat, e, s, sv, d, mode)[0][:, 0]
    gaps = [np.linalg.norm(stress(De0 / k, 2) - stress(De0 / k, 3)) for k in (1, 2, 4, 8)]
    ratios = np.array(gaps[:-1]) / np.array(gaps[1:])
    assert gaps[0] > 1e-3 and (ratios > 3.0).all(), (gaps, ratios)


def test_elastic_increment_is_bitwise_identical_across_modes():
    mat = _material(CASES["Hill+AF"][0])
    e0, s0, sv0 = np.zeros(6), np.zeros(6), np.r_[290., np.zeros(mat.nstatev - 1)]
    De = np.array([1e-4, 0, 0, 0, 0, 0])
    outs = [_umat(mat, e0, s0, sv0, De, m, time=0.0, dtime=1.0) for m in (1, 2, 3)]
    for o in outs[1:]:
        assert np.array_equal(o[0], outs[0][0]) and np.array_equal(o[1], outs[0][1])
        assert np.array_equal(o[3], outs[0][3])


@pytest.mark.parametrize("mechs", ["visco+plast", "plast+damage", "visco+plast+damage"])
def test_closest_point_composite_tangent_is_exact(mechs):
    """Viscoelastic branches (closed-form predict) and damage (strain equivalence) are rows
    without multipliers: they sit around the closest-point solve and keep the composite
    operator exact (same construction as the mode-2 composite test, with the plastic row now
    kinematic so that mode 2 would NOT be exact)."""
    terms = [(1500., 0.35, 3000., 1200.), (800., 0.35, 30000., 12000.)]
    parts = {"visco": Viscoelasticity(terms=terms), "plast": CASES["J2+AF"][0],
             "damage": Damage(Y_0=0.001, Y_c=0.5)}
    mat = _material(*[parts[k] for k in mechs.split("+")])
    r = _plastic_state(mat, 3)
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    Lt, fd = _tangent_and_fd(mat, e, s, sv, 3)
    assert np.linalg.norm(Lt - fd) < 1e-7 * np.linalg.norm(fd)


def test_closest_point_with_hyperelastic_block_is_exact():
    """A hyperelastic block is handed to the helper as its elastic response and evaluated at
    every iterate (the tangent moves with the strain): the mode-3 operator stays the exact
    Jacobian with a rotating normal on top (J2 + AF)."""
    mat = ModularMaterial(elasticity=NeoHookeanElasticity(mu=1100., kappa=5000.),
                          mechanisms=[CASES["J2+AF"][0]])
    r = _plastic_state(mat, 3)
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    Lt3, fd3 = _tangent_and_fd(mat, e, s, sv, 3)
    assert np.linalg.norm(Lt3 - fd3) < 1e-7 * np.linalg.norm(fd3)
    assert sv[1] > 1e-4   # the plastic row is active (p grew)


def test_tresca_degrades_to_the_cutting_plane_loop():
    """Tresca has no flow Hessian: under mode 3 the UMAT keeps the cutting-plane loop and the
    algorithmic operator, i.e. mode 3 == mode 2 to the bit."""
    pl = Plasticity(sigma_Y=20., yield_criterion=TrescaYield(), isotropic_hardening=VoceHardening(Q=10., b=50.))
    mat = _material(pl)
    r2, r3 = _plastic_state(mat, 2), _plastic_state(mat, 3)
    assert np.array_equal(r2["Stress"], r3["Stress"])
    assert np.array_equal(r2["TangentMatrix"], r3["TangentMatrix"])


def test_two_mechanisms_reload_converges_without_step_cut():
    """Two von Mises rows with different yield stresses under a large reversed reload: the
    cutting-plane loop needs the commit-time drift guard here (Aug 2026); the closest-point
    solve satisfies both flow rules by construction and takes every increment at full size."""
    p1 = Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                    isotropic_hardening=LinearIsotropicHardening(H=300.),
                    kinematic_hardening=ArmstrongFrederickHardening(C=1500., D=30.))
    p2 = Plasticity(sigma_Y=35., yield_criterion=VonMisesYield(),
                    isotropic_hardening=LinearIsotropicHardening(H=100.))
    mat = _material(p1, p2)
    steps = [StepMeca(control=_UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=10, time=1.),
             StepMeca(control=_UNIAXIAL, value=[-0.02, 0, 0, 0, 0, 0], ninc=4, time=1.),
             StepMeca(control=_UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=4, time=1.)]
    r3 = solve(steps, "MODUL", mat.props, mat.nstatev, T_init=290., tangent_mode=3)
    assert r3.status == 0 and len(r3) == 18
    s = r3["Stress"][0]
    assert s[9] > 20. and s[13] < -20. and s[-1] > 20.      # loads, reverses, reloads
    # no step cut: the recorded increments are exactly the prescribed ones
    np.testing.assert_allclose(np.diff(r3["Time"]), np.r_[np.full(9, 0.1), np.full(8, 0.25)])
    # against a fine cutting-plane reference, the coarse closest-point history is first-order close
    ref = solve([StepMeca(control=_UNIAXIAL, value=v, ninc=200, time=1.) for v in ([0.02, 0, 0, 0, 0, 0], [-0.02, 0, 0, 0, 0, 0], [0.02, 0, 0, 0, 0, 0])],
                "MODUL", mat.props, mat.nstatev, T_init=290., tangent_mode=2)
    assert abs(s[-1] - ref["Stress"][0, -1]) < 0.05 * abs(ref["Stress"][0, -1])


def test_two_active_rows_closest_point_tangent_is_exact():
    """Two plasticity rows active together: the exact Jacobian couples them through a
    non-symmetric local Jacobian, and the multiplier elimination of the tangent assembly must
    use its transpose (regression of the elimination fix): mode 3 vs central differences."""
    p1 = Plasticity(sigma_Y=20., yield_criterion=VonMisesYield(),
                    isotropic_hardening=LinearIsotropicHardening(H=300.),
                    kinematic_hardening=ArmstrongFrederickHardening(C=1500., D=30.))
    p2 = Plasticity(sigma_Y=35., yield_criterion=VonMisesYield(),
                    isotropic_hardening=LinearIsotropicHardening(H=100.))
    mat = _material(p1, p2)
    r = solve(StepMeca(control=_UNIAXIAL, value=[0.06, 0, 0, 0, 0, 0], ninc=20, time=1.),
              "MODUL", mat.props, mat.nstatev, T_init=290., tangent_mode=3)
    assert r.status == 0
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    assert sv[1] > 1e-3 and sv[14] > 1e-3     # both rows carried multipliers
    Lt, fd = _tangent_and_fd(mat, e, s, sv, 3)
    assert np.linalg.norm(Lt - fd) < 1e-8 * np.linalg.norm(fd)
