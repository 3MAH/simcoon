"""Numerical options of the local and global Newton loops (lot 6, all off by default).

- local complementarity (trailing props: EPICP props[6], EPKCP props[7], MODUL through
  ModularMaterial(ncp=...)): 0 the Fischer-Burmeister function, 1 min-NCP with a
  Fischer-Burmeister fallback (single surface, isotropic or Prager hardening);
- EPICP props[7] = 1: algorithmic tangent of the radial return in closed form;
- solve(precision_rel=...): stress-scaled tolerance of the global residual.
The default (absent option, or 0) is the historical behaviour, bit for bit.
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


@pytest.mark.parametrize("name, props, nstatev, slots", [
    ("EPICP", EPICP, 8, [0.0, 0.0]),
    ("EPKCP", EPKCP, 14, [0.0]),
])
def test_zero_options_are_the_default(name, props, nstatev, slots):
    ref = run(name, props, nstatev)
    res = run(name, props + slots, nstatev)
    for key in ("Stress", "Strain", "Statev", "TangentMatrix"):
        np.testing.assert_array_equal(res[key], ref[key])


def test_modular_default_props_unchanged():
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5))])
    mat_min = ModularMaterial(mat.elasticity, mat.mechanisms, ncp="min")
    np.testing.assert_array_equal(mat_min.props[:-1], mat.props)
    assert mat_min.props[-1] == 1.0
    with pytest.raises(ValueError):
        ModularMaterial(mat.elasticity, mat.mechanisms, ncp="smooth")


def test_min_ncp_epicp():
    # linear hardening: the min-NCP step is the exact return, so the result is the FB one to the
    # local tolerance; power law: within the local tolerance of the return
    for props, tol in ((EPICP_LIN, 1e-10), (EPICP, 1e-7)):
        ref = run("EPICP", props, 8)
        res = run("EPICP", props + [1.0], 8)
        assert rel(res["Stress"], ref["Stress"]) < tol
        assert rel(res["Statev"][1], ref["Statev"][1]) < 10 * tol


def test_closed_form_tangent_epicp():
    ref = run("EPICP", EPICP, 8)
    res = run("EPICP", EPICP + [0.0, 1.0], 8)
    assert rel(res["Stress"], ref["Stress"]) < 1e-8
    assert rel(res["TangentMatrix"], ref["TangentMatrix"]) < 1e-7
    both = run("EPICP", EPICP + [1.0, 1.0], 8)
    assert rel(both["Stress"], ref["Stress"]) < 1e-7


def test_min_ncp_batch_umat():
    # the same option through sim.umat (props per call)
    N = 200
    rng = np.random.default_rng(2)
    F = np.asfortranarray
    Detot = F(4e-3 * rng.standard_normal((6, N)))
    out = []
    for slots in ([], [1.0, 1.0]):
        props = F(np.asarray(EPICP + slots, float).reshape(-1, 1))
        sig, sv, Wm, Lt = sim.umat("EPICP", F(np.zeros((6, N))), Detot, np.array([]), np.array([]), F(np.zeros((6, N))),
                                   F(np.repeat(np.eye(3)[:, :, None], N, axis=2)), props, F(np.zeros((8, N))),
                                   0.0, 1.0, F(np.zeros((4, N))), temp=np.array([293.15]), n_threads=1)
        out.append((sig, Lt))
    assert rel(out[1][0], out[0][0]) < 1e-7
    assert rel(out[1][1], out[0][1]) < 1e-6


def test_min_ncp_modular_and_epkcp_adapter():
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3, alpha=1e-5),
                          [Plasticity(sigma_Y=300.0, isotropic_hardening=PowerLawHardening(k=1000.0, m=0.5),
                                      kinematic_hardening=PragerHardening(C=15000.0))], ncp="min")
    ref = run("EPKCP", EPKCP, 14)
    res = run("EPKCP", EPKCP + [1.0], 14)
    mod = run("MODUL", mat.props, mat.nstatev)
    assert rel(res["Stress"], ref["Stress"]) < 1e-6
    assert rel(mod["Stress"], res["Stress"]) < 1e-12   # the adapter forwards the option


def test_min_ncp_refuses_armstrong_frederick():
    # the cutting plane integrates AF back stresses along its iterates: min-NCP would change the
    # return far beyond the tolerance on coarse increments, so MODUL refuses the combination
    from simcoon.modular import ArmstrongFrederickHardening
    mat = ModularMaterial(IsotropicElasticity(C1=200000.0, C2=0.3),
                          [Plasticity(sigma_Y=300.0, kinematic_hardening=ArmstrongFrederickHardening(C=30000.0, D=200.0))],
                          ncp="min")
    with pytest.raises(Exception):
        run("MODUL", mat.props, mat.nstatev)


def test_option_values_are_checked():
    with pytest.raises(Exception):
        run("EPICP", EPICP + [2.0], 8)
    with pytest.raises(Exception):
        run("EPICP", EPICP + [0.0, 0.5], 8)
    with pytest.raises(Exception):
        run("EPKCP", EPKCP + [3.0], 14)


def test_relative_global_tolerance():
    ref = run("EPICP", EPICP, 8)
    same = run("EPICP", EPICP, 8, precision_rel=0.0)
    for key in ("Stress", "Strain", "Statev"):
        np.testing.assert_array_equal(same[key], ref[key])
    res = run("EPICP", EPICP, 8, precision_rel=1e-7)
    assert rel(res["Stress"], ref["Stress"]) < 1e-5
    with pytest.raises(Exception):
        run("EPICP", EPICP, 8, precision_rel=-1.0)
