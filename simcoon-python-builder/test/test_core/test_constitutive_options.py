"""Numerical options of the local and global Newton loops (all off by default).

- solve(ncp=...), sim.umat(ncp=...): complementarity function of the local cutting plane,
  'fischer_burmeister' (default) or 'min' (min-NCP with a Fischer-Burmeister fallback; single
  surface, isotropic or Prager hardening), read by EPICP, EPKCP and MODUL;
- solve(precision_rel=...): stress-scaled tolerance of the global residual.
The default is the historical behaviour, bit for bit.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import (
    IsotropicElasticity,
    ModularMaterial,
    Plasticity,
    PowerLawHardening,
    PragerHardening,
)
from simcoon.solver import Block, StepMeca

EPICP = [200000.0, 0.3, 1e-5, 300.0, 1000.0, 0.5]
EPICP_LIN = [200000.0, 0.3, 1e-5, 300.0, 2000.0, 1.0]   # linear hardening: one return step
EPKCP = [200000.0, 0.3, 1e-5, 300.0, 1000.0, 0.5, 10000.0]
EPCHA = [140000.0, 0.3, 1e-6, 62.859017, 416.004456, 4.788635, 30382.293921, 172.425687, 195142.490843,
         3012.614659]


def cyclic(emax=0.02, ninc=50, ncycle=2):
    uni = ["strain"] + ["stress"] * 5
    steps = []
    for _ in range(ncycle):
        steps.append(StepMeca(control=uni, value=[emax, 0, 0, 0, 0, 0], time=1.0, ninc=ninc))
        steps.append(StepMeca(control=uni, value=[-emax, 0, 0, 0, 0, 0], time=1.0, ninc=2 * ninc))
    return Block(steps=steps)


def run(name, props, nstatev, **kw):
    res = sim.solver.solve(cyclic(), name, np.asarray(props, float), nstatev, T_init=293.15, **kw)
    assert res.status == 0
    return res


def rel(a, b):
    return np.max(np.abs(np.asarray(a) - np.asarray(b))) / np.max(np.abs(np.asarray(b)))


def test_default_ncp_is_the_default():
    ref = run("EPKCP", EPKCP, 14)
    res = run("EPKCP", EPKCP, 14, ncp="fischer_burmeister")
    for key in ("Stress", "Strain", "Statev", "TangentMatrix"):
        np.testing.assert_array_equal(res[key], ref[key])


def test_min_ncp_epicp():
    # linear hardening: the min-NCP step is the exact return, so the result is the FB one to the
    # local tolerance; power law: within the local tolerance of the return
    for props, tol in ((EPICP_LIN, 1e-10), (EPICP, 1e-7)):
        ref = run("EPICP", props, 8)
        res = run("EPICP", props, 8, ncp="min")
        assert rel(res["Stress"], ref["Stress"]) < tol
        assert rel(res["Statev"][1], ref["Statev"][1]) < 10 * tol


def test_closed_form_tangent_epicp():
    # tangent_mode 2 takes the radial-return tangent in closed form: the exact operator of the
    # closest-point update (integrator="closest_point", generic assembly) at the same converged state
    ref = run("EPICP", EPICP, 8, integrator="closest_point")
    res = run("EPICP", EPICP, 8)
    assert rel(res["Stress"], ref["Stress"]) < 1e-8
    assert rel(res["TangentMatrix"], ref["TangentMatrix"]) < 1e-7


def test_closest_point_serves_the_continuum_and_elastic_operators():
    """integrator="closest_point" with tangent_mode "continuum" or "none": the same converged state as
    its algorithmic run; the continuum operator is the rate operator of the converged state, equal to
    the cutting plane's on a radial return (same state), away from the exact operator on a plastic
    step; "none" hands the elastic operator."""
    cp = {m: run("EPICP", EPICP, 8, integrator="closest_point", tangent_mode=m) for m in ("none", "continuum", "algorithmic")}
    cc = run("EPICP", EPICP, 8, tangent_mode="continuum")
    for m in ("none", "continuum"):   # the same state to the global Newton tolerance
        assert rel(cp[m]["Stress"], cp["algorithmic"]["Stress"]) < 1e-6
    assert rel(cp["continuum"]["TangentMatrix"], cc["TangentMatrix"]) < 1e-6
    assert rel(cp["continuum"]["TangentMatrix"], cp["algorithmic"]["TangentMatrix"]) > 1e-3
    L = run("EPICP", EPICP, 8, tangent_mode="none")["TangentMatrix"]
    np.testing.assert_array_equal(cp["none"]["TangentMatrix"], L)
    # the modular engine and a kinematic law alike
    mat = ModularMaterial(IsotropicElasticity(C1=70000.0, C2=0.3, alpha=0.0),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5))])
    for name, props, nstatev in (("MODUL", mat.props, mat.nstatev), ("EPCHA", EPCHA, 33)):
        a = run(name, props, nstatev, integrator="closest_point")
        c = run(name, props, nstatev, integrator="closest_point", tangent_mode="continuum")
        assert rel(c["Stress"], a["Stress"]) < 1e-6
        assert rel(c["TangentMatrix"], a["TangentMatrix"]) > 1e-3
    # sim.umat takes the same pair
    col = lambda v: np.asfortranarray(np.asarray(v, float).reshape(-1, 1))
    out = sim.umat("EPICP", col(np.zeros(6)), col([4e-3, 0, 0, 0, 0, 0]), np.empty(0), np.empty(0), col(np.zeros(6)),
                   np.asfortranarray(np.eye(3)[:, :, None]), col(EPICP), col(np.zeros(8)), 0.0, 1.0, col(np.zeros(4)),
                   temp=np.full(1, 293.15), n_threads=1, integrator="closest_point", tangent_mode=1)
    assert np.isfinite(out[3]).all()


@pytest.mark.parametrize("method", ["active_set"])
@pytest.mark.parametrize("integrator", ["cutting_plane", "closest_point"])
def test_other_complementarity_functions_find_the_same_return(method, integrator):
    """active_set (the min step on every row, coupled surfaces together): the same return as
    Fischer-Burmeister to the local tolerance on a single surface and on a coupled two-surface
    composition under the closest point (one root); under the cutting plane, a different iteration
    path on coupled surfaces, it lands within the path dependence of that integrator."""
    iso = IsotropicElasticity(C1=200000.0, C2=0.3, alpha=1e-5)
    one = ModularMaterial(iso, [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5))])
    two = ModularMaterial(iso, [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5)),
                                Plasticity(sigma_Y=320.0, isotropic_hardening=PowerLawHardening(k=3000.0, m=0.7))])
    for name, props, nstatev, tol in (("MODUL", one.props, one.nstatev, 1e-8), ("EPICP_closed", EPICP, 8, 1e-8),
                                      ("MODUL", two.props, two.nstatev, 1e-8 if (integrator == "closest_point" or method != "active_set") else 1e-3)):
        ref = run(name, props, nstatev, integrator=integrator)
        res = run(name, props, nstatev, integrator=integrator, ncp=method)
        assert rel(res["Stress"], ref["Stress"]) < tol, (name, method, integrator)


def test_sma_laws_read_the_complementarity_function():
    """The unified SMA laws take their local step through ncp_step: under the closest point every
    function lands on the Fischer-Burmeister root (one backward-Euler system); min falls back to
    Fischer-Burmeister on their coupled rows, bit for bit."""
    sma = [0, 67538.0, 67538.0, 0.349, 0.349, 1e-6, 1e-6, 0.0, 0.0418, 0.021, 0.0, 10.0, 10.0, 250.0, 230.0,
           260.0, 280.0, 0.2, 0.2, 0.2, 0.2, 300.0, 1.4, 2.0, 1e-6, 1e-3, 1.0, 1e8]
    uni = ["strain"] + ["stress"] * 5
    path = Block(steps=[StepMeca(control=uni, value=[0.05, 0, 0, 0, 0, 0], time=1.0, ninc=50, Dn_mini=0.01),
                        StepMeca(control=["stress"] * 6, value=[0] * 6, time=1.0, ninc=50, Dn_mini=0.01)])
    go = lambda **kw: sim.solver.solve(path, "SMADI", np.asarray(sma, float), 50, T_init=320.0, **kw)
    for integrator in ("cutting_plane", "closest_point"):
        ref = go(integrator=integrator)
        np.testing.assert_array_equal(go(integrator=integrator, ncp="min")["Stress"], ref["Stress"])
        tol = 1e-8 if integrator == "closest_point" else 1e-5   # the cutting plane: another iteration path
        assert rel(go(integrator=integrator, ncp="active_set")["Stress"], ref["Stress"]) < tol


def test_complementarity_function_admission():
    """The semi-smooth functions (min, active_set) are refused with an Armstrong-Frederick back stress
    under the cutting plane only; an unknown name (the penalized Fischer-Burmeister function, measured
    and removed in 2.5) and a law that does not read ncp are refused."""
    from simcoon.modular import ArmstrongFrederickHardening
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3, alpha=1e-5),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5),
                                      kinematic_hardening=ArmstrongFrederickHardening(C=30000.0, D=100.0))])
    with pytest.raises(Exception, match="Armstrong-Frederick"):
        run("MODUL", mat.props, mat.nstatev, ncp="active_set")
    assert run("MODUL", mat.props, mat.nstatev, ncp="active_set", integrator="closest_point").status == 0
    for unknown in ("smoothed", "penalized_fischer_burmeister"):
        with pytest.raises(ValueError):
            run("MODUL", mat.props, mat.nstatev, ncp=unknown)
    with pytest.raises(ValueError, match="does not read"):   # a dedicated kernel without the option
        run("ZENER_closed", [3000.0, 0.35, 0.0, 1500.0, 0.35, 3000.0, 1200.0], 8, ncp="active_set")


def test_min_ncp_batch_umat():
    # the same option through sim.umat
    N = 200
    rng = np.random.default_rng(2)
    F = np.asfortranarray
    Detot = F(4e-3 * rng.standard_normal((6, N)))
    props = F(np.asarray(EPICP, float).reshape(-1, 1))
    out = []
    for ncp in ("fischer_burmeister", "min"):
        sig, sv, Wm, Lt = sim.umat("EPICP", F(np.zeros((6, N))), Detot, np.array([]), np.array([]), F(np.zeros((6, N))),
                                   F(np.repeat(np.eye(3)[:, :, None], N, axis=2)), props, F(np.zeros((8, N))),
                                   0.0, 1.0, F(np.zeros((4, N))), temp=np.array([293.15]), n_threads=1, ncp=ncp)
        out.append((sig, Lt))
    assert rel(out[1][0], out[0][0]) < 1e-7
    assert rel(out[1][1], out[0][1]) < 1e-6


def test_min_ncp_modular_and_epkcp_adapter():
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3, alpha=1e-5),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5),
                                      kinematic_hardening=PragerHardening(C=15000.0))])
    ref = run("EPKCP", EPKCP, 14)
    res = run("EPKCP", EPKCP, 14, ncp="min")
    mod = run("MODUL", mat.props, mat.nstatev, ncp="min")
    assert rel(res["Stress"], ref["Stress"]) < 1e-6
    assert rel(mod["Stress"], res["Stress"]) < 1e-12   # the adapter forwards the option


@pytest.mark.parametrize("extra", ["prony", "damage"])
def test_min_ncp_modular_with_inert_rows(extra):
    # viscoelastic and damage rows are always satisfied and uncoupled: min-NCP applies to the
    # yield surface next to them (the option has an effect) and gives the same return
    from simcoon.modular import Damage, Viscoelasticity
    other = (Viscoelasticity(terms=((20000.0, 0.3, 2000.0, 2000.0),) * 2) if extra == "prony"
             else Damage(Y_0=0.1, Y_c=5.0))
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3, alpha=1e-5),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5)), other])
    ref = run("MODUL", mat.props, mat.nstatev)
    res = run("MODUL", mat.props, mat.nstatev, ncp="min")
    assert not np.array_equal(res["Stress"], ref["Stress"])
    assert rel(res["Stress"], ref["Stress"]) < 1e-6


def test_min_ncp_refuses_armstrong_frederick_under_the_cutting_plane_only():
    # the cutting plane integrates AF back stresses along its iterates: min-NCP would change the
    # return far beyond the tolerance on coarse increments, so MODUL refuses the combination there;
    # the closest point solves one backward-Euler system whose root does not depend on the step
    # function, so it takes min-NCP (measured 1.5e-10 against Fischer-Burmeister on a coarse
    # non-proportional path, where the cutting plane sits 4e-2 away)
    from simcoon.modular import ArmstrongFrederickHardening
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3),
                          [Plasticity(sigma_Y=300.0, kinematic_hardening=ArmstrongFrederickHardening(C=30000.0, D=200.0))])
    run("MODUL", mat.props, mat.nstatev)
    with pytest.raises(Exception, match="Armstrong-Frederick"):
        run("MODUL", mat.props, mat.nstatev, ncp="min")
    fb = run("MODUL", mat.props, mat.nstatev, integrator="closest_point")
    mn = run("MODUL", mat.props, mat.nstatev, integrator="closest_point", ncp="min")
    assert rel(mn["Stress"], fb["Stress"]) < 1e-8


def test_ncp_values_and_laws_are_checked():
    with pytest.raises(ValueError):
        run("EPICP", EPICP, 8, ncp="smooth")
    with pytest.raises(ValueError):
        run("EPICP", EPICP, 8, ncp=1)
    with pytest.raises(ValueError, match="EPCHA"):
        run("EPCHA", EPCHA, 33, ncp="min")
    # the closest point takes the min step on its reduced system: the same root
    for name in ("EPICP", "EPICP_closed"):
        fb = run(name, EPICP, 8, integrator="closest_point")
        mn = run(name, EPICP, 8, ncp="min", integrator="closest_point")
        assert not np.array_equal(mn["Stress"], fb["Stress"]) or name == "EPICP_closed"
        assert rel(mn["Stress"], fb["Stress"]) < 1e-8
    with pytest.raises(ValueError, match="ncp"):
        sim.umat("EPCHA", np.zeros((6, 1)), np.zeros((6, 1)), np.array([]), np.array([]), np.zeros((6, 1)),
                 np.eye(3).reshape(3, 3, 1), np.asarray(EPCHA, float), np.zeros((33, 1)), 0.0, 1.0,
                 np.zeros((4, 1)), ncp="min")


def test_relative_global_tolerance():
    ref = run("EPICP", EPICP, 8)
    same = run("EPICP", EPICP, 8, precision_rel=0.0)
    for key in ("Stress", "Strain", "Statev"):
        np.testing.assert_array_equal(same[key], ref[key])
    res = run("EPICP", EPICP, 8, precision_rel=1e-7)
    assert rel(res["Stress"], ref["Stress"]) < 1e-5
    with pytest.raises(Exception):
        run("EPICP", EPICP, 8, precision_rel=-1.0)
