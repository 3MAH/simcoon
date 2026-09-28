"""Objectivity of the finite-strain route under superposed rigid rotation.

A kernel must see its total strain, its start stress and its internal variables in the SAME
configuration. The route used to hand the total-form kernels a start strain not yet
transported by the current increment's rotation, while their internal variables were
rotated by it inside the kernel and their start stress by the previous increment's
rotation. Under a rigid rotation that produced spurious elastic strain -- 1.3 % stress
error and spurious plastic flow after 90 degrees in 100 increments -- and, through the
start stress, a ~60 % error in Biot stress control under spin.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca


def _v2t(v):
    return np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])


def _rz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


@pytest.mark.parametrize("corate", [0, 2, 3, 4])
def test_rigid_rotation_of_a_plastic_state_is_exact(corate):
    """EPICP stretched into plasticity, then rotated rigidly by 90 degrees: the stress just
    rotates, and no plastic strain or dissipation is created."""
    props = np.array([200000., 0.3, 0., 300., 1000., 0.5])
    F1 = np.diag([np.exp(0.05), 1., 1.])
    R = _rz(np.pi / 2)
    steps = [StepMeca(control="F", value=F1.ravel().tolist(), time=1., ninc=50,
                      Dn_init=1., Dn_mini=1e-4),
             StepMeca(control="F", value=(R @ F1).ravel().tolist(), time=1., ninc=100,
                      Dn_init=1., Dn_mini=1e-4)]
    r = sim.solver.solve(Block(steps=steps, control_type="F"), "EPICP", props, 8,
                         T_init=290., corate=corate)
    i1 = np.argmin(np.abs(r["Time"] - 1.0))
    assert r["Statev"][1, i1] > 1e-2, "the stretch must have yielded"
    tau1, tau2 = _v2t(r["Kirchhoff"][:, i1]), _v2t(r["Kirchhoff"][:, -1])
    np.testing.assert_allclose(tau2, R @ tau1 @ R.T, atol=1e-10 * np.abs(tau1).max())
    assert abs(r["Statev"][1, -1] - r["Statev"][1, i1]) < 1e-12
    assert abs(r["Wm"][3, -1] - r["Wm"][3, i1]) < 1e-8 * r["Wm"][3, i1]


@pytest.mark.parametrize("corate", [0, 2, 3])
@pytest.mark.parametrize("control_type", [2, 3, 4])
@pytest.mark.parametrize("umat", ["SNTVE", "ELISO"])
def test_stress_hold_under_spin(umat, control_type, corate):
    """Stress ramp, then every component held while a z-spin rotates the body: the material
    state does not change (PKII constant) and tau rotates with the spin, for every corate --
    the Kirchhoff increments of ct 3 are prescribed in the polar frame F = V R rebuilds. The tolerance on tau covers the Hughes-Winget (Cayley) integration of the
    spin, ~ theta * dtheta^2 / 12 = 4e-6 here."""
    th = 0.8
    target = 30000. if umat == "SNTVE" else 15000.
    w = [[0., -th, 0.], [th, 0., 0.], [0., 0., 0.]]
    steps = [StepMeca(control="stress", value=[target, 0, 0, 0, 0, 0], time=1., ninc=50,
                      Dn_init=1., Dn_mini=1e-3),
             StepMeca(control="stress", value=[target, 0, 0, 0, 0, 0], time=1., ninc=100,
                      Dn_init=1., Dn_mini=1e-3, BC_w=w)]
    r = sim.solver.solve(Block(steps=steps, control_type=control_type), umat,
                         np.array([100000., 0.3, 0.]), 1, T_init=290., corate=corate,
                         record_tangent=False)
    i1 = np.argmin(np.abs(r["Time"] - 1.0))
    S1, S2 = r["PKII"][:, i1], r["PKII"][:, -1]
    np.testing.assert_allclose(S2, S1, atol=1e-8 * np.abs(S1).max())
    R = _rz(th)
    tau1, tau2 = _v2t(r["Kirchhoff"][:, i1]), _v2t(r["Kirchhoff"][:, -1])
    np.testing.assert_allclose(tau2, R @ tau1 @ R.T, atol=1e-4 * np.abs(tau1).max())


# ----- corate 4: the genuine convected (Truesdell / Oldroyd) rate --------------------------

def _simple_shear(umat, props, ninc, corate=4, gamma=1.0):
    st = StepMeca(control="F", value=[1., gamma, 0., 0., 1., 0., 0., 0., 1.], time=1., ninc=ninc,
                  Dn_init=1., Dn_mini=1e-5)
    return sim.solver.solve(Block(steps=[st], control_type="F"), umat, np.asarray(props, float),
                            1, T_init=290., corate=corate, record_tangent=False)


def test_truesdell_total_form_is_the_almansi_law():
    """Corate 4 transports the strain lower-convected with the closed-form Almansi increment,
    so a total-form kernel gives tau = L : e_A(F) exactly, whatever the increment."""
    E, nu = 70000., 0.3
    L = np.asarray(sim.L_iso([E, nu], "Enu"))
    for ninc in (10, 100):
        r = _simple_shear("ELISO", [E, nu, 0.], ninc)
        F = r["F"][:, :, -1]
        e_A = 0.5 * (np.eye(3) - np.linalg.inv(F @ F.T))
        ref = L @ np.asarray(sim.t2v_strain(e_A)).ravel()
        np.testing.assert_allclose(r["Kirchhoff"][:, -1], ref, atol=1e-10 * np.abs(ref).max())


def test_truesdell_rate_form_converges_to_the_oldroyd_solution():
    """HYPOO at corate 4 integrates tau_dot = L tau + tau L^T + lambda tr(D) I + 2 mu D; in simple
    shear the exact solution is tau_12 = mu gamma, tau_11 = mu gamma^2 (the Truesdell normal
    stress). The explicit update converges at first order."""
    E, nu = 70000., 0.3
    mu = E / (2 * (1 + nu))
    iso = [E] * 3 + [nu] * 3 + [mu] * 3 + [0.] * 3
    err = []
    for ninc in (200, 800):
        tau = _v2t(_simple_shear("HYPOO", iso, ninc)["Kirchhoff"][:, -1])
        err.append(max(abs(tau[0, 1] - mu), abs(tau[0, 0] - mu)) / mu)
    assert err[1] < 5e-3
    assert err[1] < 0.3 * err[0], "first-order convergence to the Oldroyd solution"
