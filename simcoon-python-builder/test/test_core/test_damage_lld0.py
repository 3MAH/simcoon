"""LLDM0 (Ladeveze-Le Dantec ply damage, plasticity in effective stress) from Python: the solver
under mixed boundary conditions and sim.umat, 10 state variables."""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, solve

T0 = 293.15
# carbon/epoxy, nu_TL the MINOR Poisson ratio
PROPS = np.array([1, 138000., 11000., 0.0223, 0.40, 5500., -0.3e-6, 30e-6,
                  0.10, 0.30, 0.15, 0.40, 1.20, 2.5, 2.0, 50., 0.5, 300.])


@pytest.mark.parametrize("axis", [1, 2, 3])
def test_off_axis_tension_damages_and_balances_the_work(axis):
    """Tension along the transverse direction of the ply, other stresses free: both damage
    variables grow, and Wm = Wm_r + Wm_ir + Wm_d along the path, whatever the fibre axis."""
    props = PROPS.copy()
    props[0] = axis
    k = {1: 1, 2: 0, 3: 0}[axis]   # a direction transverse to the fibres
    ctl = ["stress"] * 6
    ctl[k] = "strain"
    val = [0.] * 6
    val[k] = 0.008
    res = solve(Block(steps=[StepMeca(control=ctl, value=val, ninc=80)], control_type="small_strain"),
                "LLDM0", props, 10, T_init=T0)
    sv, wm, s = np.asarray(res["Statev"]), np.asarray(res["Wm"]), np.asarray(res["Stress"])
    assert sv.shape[0] == 10
    assert sv[1, -1] > 0.1 and sv[2, -1] > 0.1
    assert np.all(np.diff(sv[1]) >= -1e-14) and np.all(np.diff(sv[2]) >= -1e-14)
    np.testing.assert_allclose(wm[0], wm[1] + wm[2] + wm[3], rtol=0, atol=1e-10 * np.abs(wm[0]).max())
    free = [i for i in range(6) if i != k]
    assert np.abs(s[free, -1]).max() < 1e-6 * np.abs(s[k, -1])


def test_fibre_axes_give_the_rotated_response():
    """The same ply with its fibres along 1 and along 3, loaded along its own transverse and shear
    directions: identical stresses and damage."""
    out = {}
    for axis, k_t, k_s in [(1, 1, 3), (3, 1, 5)]:   # fibres 1: transverse 2, shear 12; fibres 3: transverse 2, shear 23 (= 3-2)
        props = PROPS.copy()
        props[0] = axis
        ctl, val = ["stress"] * 6, [0.] * 6
        ctl[k_t], ctl[k_s] = "strain", "strain"
        val[k_t], val[k_s] = 0.004, 0.008
        res = solve(Block(steps=[StepMeca(control=ctl, value=val, ninc=60)], control_type="small_strain"),
                    "LLDM0", props, 10, T_init=T0)
        s = np.asarray(res["Stress"])[:, -1]
        out[axis] = (s[k_t], abs(s[k_s]), np.asarray(res["Statev"])[1:4, -1])
    np.testing.assert_allclose(out[3][0], out[1][0], rtol=1e-8)
    np.testing.assert_allclose(out[3][1], out[1][1], rtol=1e-8)
    np.testing.assert_allclose(out[3][2], out[1][2], rtol=1e-8, atol=1e-12)


def test_umat_with_ten_state_variables():
    """sim.umat runs the law with its documented 10 state variables."""
    col = lambda v: np.asfortranarray(np.asarray(v, float).reshape(-1, 1))
    e, s, sv, wm = np.zeros(6), np.zeros(6), np.zeros(10), np.zeros(4)
    for i in range(12):
        de = np.array([0., 2.5e-4, 0., 5e-4, 0., 0.])
        s, sv, wm, Lt = sim.umat("LLDM0", col(e), col(de), np.array([]), np.array([]), col(s),
                                 np.asfortranarray(np.eye(3).reshape(3, 3, 1)), col(PROPS),
                                 np.asfortranarray(sv.reshape(-1, 1)), float(i), 0.01, col(wm), start=(i == 0))
        e = e + de
        s, sv, wm = np.asarray(s).ravel(), np.asarray(sv).ravel(), np.asarray(wm).ravel()
    assert sv[2] > 0.05
    assert np.all(np.isfinite(np.asarray(Lt)))
