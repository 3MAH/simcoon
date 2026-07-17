"""Tests for the point-wise thermomechanical UMAT binding sim.umat_T."""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import StepThermomeca, solve

ELISO_T_PROPS = [1.0e-9, 1.0, 70000.0, 0.3, 1.0e-5]  # rho c_p E nu alpha


def _batch(n, nstatev=1):
    return dict(
        etot=np.zeros((6, n), order="F"),
        Detot=np.zeros((6, n), order="F"),
        sigma=np.zeros((6, n), order="F"),
        DR=np.asfortranarray(np.tile(np.eye(3)[:, :, None], (1, 1, n))),
        statev=np.zeros((nstatev, n), order="F"),
        Wm=np.zeros((4, n), order="F"),
        Wt=np.zeros((3, n), order="F"),
    )


def test_eliso_T_constrained_thermal_stress():
    rho, c_p, E, nu, alpha = ELISO_T_PROPS
    n = 4
    b = _batch(n)
    props = np.array(ELISO_T_PROPS, dtype=float).reshape(-1, 1)
    T = np.full(n, 293.15)
    DT = np.full(n, 10.0)
    sigma, statev, Wm, Wt, r, dSdE, dSdT, drdE, drdT = sim.umat_T(
        "ELISO", b["etot"], b["Detot"], b["sigma"], b["DR"], props, b["statev"],
        0.0, 1.0, b["Wm"], b["Wt"], T, DT)
    # fully constrained heating: sigma = -E alpha DT / (1 - 2 nu) on the diagonal
    expected = -E * alpha * 10.0 / (1.0 - 2.0 * nu)
    np.testing.assert_allclose(sigma[0:3, :], expected, rtol=1e-10)
    np.testing.assert_allclose(sigma[3:6, :], 0.0, atol=1e-12)
    # shapes of the coupled tangents
    assert dSdE.shape == (6, 6, n)
    assert dSdT.shape == (6, n)
    assert drdE.shape == (6, n)
    assert drdT.shape in ((n,), (n, 1))
    assert r.shape in ((n,), (n, 1))
    # dSdT = -L : alpha
    L = sim.L_iso([E, nu], "Enu")
    a = np.array([alpha, alpha, alpha, 0, 0, 0])
    np.testing.assert_allclose(dSdT[:, 0], -L @ a, rtol=1e-10)


def test_umat_T_per_point_props():
    n = 2
    b = _batch(n)
    props = np.column_stack([
        np.array([1e-9, 1.0, 70000.0, 0.3, 1.0e-5]),
        np.array([1e-9, 1.0, 210000.0, 0.3, 1.0e-5]),
    ])
    T = np.full(n, 293.15)
    DT = np.full(n, 10.0)
    Detot = np.zeros((6, n), order="F")
    Detot[0, :] = 0.001  # uniaxial strain increment on both points
    sigma = sim.umat_T("ELISO", b["etot"], Detot, b["sigma"], b["DR"], props,
                       b["statev"], 0.0, 1.0, b["Wm"], b["Wt"], T, DT * 0.0)[0]
    # stiffer point develops proportionally larger stress
    np.testing.assert_allclose(sigma[0, 1] / sigma[0, 0], 3.0, rtol=1e-10)


def test_umat_T_consistent_with_block_solver():
    # drive the point-wise UMAT with the kinematics found by the block solver
    # (stress-free thermal expansion): it must reproduce the near-zero stress
    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=10, T_final=340.0)
    res = solve(step, "ELISO", ELISO_T_PROPS, 1, T_init=290.0)
    n = 1
    b = _batch(n)
    props = np.array(ELISO_T_PROPS, dtype=float).reshape(-1, 1)
    # last increment: state at inc -2, increments to inc -1
    etot = np.asfortranarray(res["Strain"][:, [-2]])
    Detot = np.asfortranarray(res["Strain"][:, [-1]] - res["Strain"][:, [-2]])
    statev = np.asfortranarray(res["Statev"][:, [-2]])
    sigma0 = np.asfortranarray(res["Stress"][:, [-2]])
    T = np.array([res["Temp"][-2]])
    DT = np.array([res["Temp"][-1] - res["Temp"][-2]])
    sigma_pt = sim.umat_T("ELISO", etot, Detot, sigma0, b["DR"], props, statev,
                          0.9, 0.1, b["Wm"], b["Wt"], T, DT)[0]
    np.testing.assert_allclose(sigma_pt[:, 0], res["Stress"][:, -1], atol=1e-6)


def test_umat_T_unknown_name_raises():
    b = _batch(1)
    props = np.array(ELISO_T_PROPS, dtype=float).reshape(-1, 1)
    with pytest.raises(Exception):
        sim.umat_T("XXXXX", b["etot"], b["Detot"], b["sigma"], b["DR"], props,
                   b["statev"], 0.0, 1.0, b["Wm"], b["Wt"],
                   np.array([290.0]), np.array([0.0]))


def test_umat_T_1d_props_accepted():
    # a 1-D (nprops,) props array is a valid shared-props input
    n = 3
    b = _batch(n)
    props_1d = np.array(ELISO_T_PROPS, dtype=float)
    T = np.full(n, 293.15)
    DT = np.full(n, 10.0)
    sigma = sim.umat_T("ELISO", b["etot"], b["Detot"], b["sigma"], b["DR"],
                       props_1d, b["statev"], 0.0, 1.0, b["Wm"], b["Wt"], T, DT)[0]
    expected = -70000.0 * 1.0e-5 * 10.0 / (1.0 - 2.0 * 0.3)
    np.testing.assert_allclose(sigma[0:3, :], expected, rtol=1e-10)


def test_umat_T_shape_mismatches_raise():
    n = 4
    b = _batch(n)
    T = np.full(n, 293.15)
    DT = np.zeros(n)
    # props with 1 < n_cols < nb_points
    bad_props = np.tile(np.array(ELISO_T_PROPS, dtype=float).reshape(-1, 1), (1, 2))
    with pytest.raises(ValueError):
        sim.umat_T("ELISO", b["etot"], b["Detot"], b["sigma"], b["DR"], bad_props,
                   b["statev"], 0.0, 1.0, b["Wm"], b["Wt"], T, DT)
    # T shorter than the batch
    props = np.array(ELISO_T_PROPS, dtype=float).reshape(-1, 1)
    with pytest.raises(ValueError):
        sim.umat_T("ELISO", b["etot"], b["Detot"], b["sigma"], b["DR"], props,
                   b["statev"], 0.0, 1.0, b["Wm"], b["Wt"], T[:2], DT)
    # statev column count mismatch
    with pytest.raises(ValueError):
        sim.umat_T("ELISO", b["etot"], b["Detot"], b["sigma"], b["DR"], props,
                   np.zeros((1, 2), order="F"), 0.0, 1.0, b["Wm"], b["Wt"], T, DT)
