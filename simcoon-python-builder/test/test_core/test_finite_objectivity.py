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


_PLASTIC = {
    "EPICP": (np.array([200000., 0.3, 0., 300., 1000., 0.5]), 8, 0.05),
    # EPCHA carries stress-like back-stresses X_1, X_2: they used not to be rotated at all
    "EPCHA": (np.array([210000.0, 0.3, 0.0, 300.0, 200.0, 20.0, 30000.0, 172.0, 19500.0, 301.0]), 33, 0.02),
}


@pytest.mark.parametrize("umat, corate", [("EPICP", c) for c in (0, 2, 3, 4)]
                         + [("EPCHA", c) for c in (0, 2, 3)])
def test_rigid_rotation_of_a_plastic_state_is_exact(corate, umat):
    """A plastic prestretch, then a rigid 90 degree rotation: the stress just rotates, and no
    plastic strain or dissipation is created."""
    props, nstatev, eps = _PLASTIC[umat]
    F1 = np.diag([np.exp(eps), 1., 1.])
    R = _rz(np.pi / 2)
    steps = [StepMeca(control="F", value=F1.ravel().tolist(), time=1., ninc=50,
                      Dn_init=1., Dn_mini=1e-4),
             StepMeca(control="F", value=(R @ F1).ravel().tolist(), time=1., ninc=100,
                      Dn_init=1., Dn_mini=1e-4)]
    r = sim.solver.solve(Block(steps=steps, control_type="F"), umat, props, nstatev,
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

def test_epcha_refuses_truesdell():
    """EPCHA stores X_i = 2/3 C_i a_i: the Truesdell rate convects a (strain) and X (stress)
    differently and would pull the pair apart, so corate 4 is refused."""
    props, nstatev, eps = _PLASTIC["EPCHA"]
    st = StepMeca(control="F", value=np.diag([np.exp(eps), 1., 1.]).ravel().tolist(), ninc=5)
    with pytest.raises(Exception, match="Truesdell"):
        sim.solver.solve(Block(steps=[st], control_type="F"), "EPCHA", props, nstatev,
                         T_init=290., corate=4)


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


# ----- small-strain kernels driven with a rotation increment (sim.umat, Abaqus DROT) --------

_SMADI = np.array([0, 67538.0, 67538.0, 0.349, 0.349, 1.0e-6, 1.0e-6, 0.0, 0.0418, 0.021, 0.0,
                   10.0, 10.0, 300.0, 290.0, 295.0, 305.0, 0.2, 0.2, 0.2, 0.2, 300.0, 0.2, 2.0,
                   1.0e-6, 1.0e-5, 1.0, 0.0])


def test_sma_rotation_increment_rotates_the_transformation_strain():
    """SMADI transformed at 323.15 K, then one call with a 90 degree DR and the total strain
    rotated accordingly: the stress and the transformation strain just rotate, the martensite
    fraction is unchanged. The rotation of ET used to be computed and discarded, giving
    xi = 0.38 -> -2.5 on this very call."""
    col = lambda v: np.asfortranarray(np.asarray(v, float).reshape(-1, 1))
    temp = np.array([323.15])

    def call(etot, Detot, sigma, DR, statev, Wm, time):
        s, sv, wm, _ = sim.umat("SMADI", col(etot), col(Detot), np.array([]), np.array([]),
                                col(sigma), DR, col(_SMADI), col(statev), time, 0.01, col(Wm),
                                temp=temp, n_threads=1)
        return np.asarray(s).ravel(), np.asarray(sv).ravel(), np.asarray(wm).ravel()

    I3 = np.eye(3).reshape(3, 3, 1).copy(order="F")
    etot, sigma, statev, Wm = np.zeros(6), np.zeros(6), np.zeros(17), np.zeros(4)
    de = np.array([5e-4, 0, 0, 0, 0, 0])
    for i in range(60):
        sigma, statev, Wm = call(etot, de, sigma, I3, statev, Wm, i * 0.01)
        etot = etot + de
    assert statev[1] > 0.05, "the loading must have transformed the material"

    R = _rz(np.pi / 2)
    rot_e = lambda v: np.asarray(sim.t2v_strain(R @ sim.v2t_strain(v) @ R.T)).ravel()
    rot_s = lambda v: np.asarray(sim.t2v_stress(R @ sim.v2t_stress(v) @ R.T)).ravel()
    s2, sv2, _ = call(rot_e(etot), np.zeros(6), sigma, R.reshape(3, 3, 1).copy(order="F"),
                      statev, Wm, 0.6)
    np.testing.assert_allclose(s2, rot_s(sigma), atol=1e-9 * np.abs(sigma).max())
    np.testing.assert_allclose(sv2[2:8], rot_e(statev[2:8]), atol=1e-12)
    assert sv2[1] == pytest.approx(statev[1], abs=1e-12)


# ----- anisotropy axes follow the material (material frame of the box kernels) -------------

_ORTHO = np.array([70000., 30000., 15000., 0.3, 0.3, 0.3, 8000., 6000., 5000., 0., 0., 0.])


@pytest.mark.parametrize("orientation", [(0., 0., 0.), (30., 0., 0.)], ids=["aligned", "psi30"])
@pytest.mark.parametrize("corate", [0, 1, 2, 3, 4, 5])
@pytest.mark.parametrize("umat", ["ELORT", "HYPOO"])
def test_orthotropic_axes_follow_a_rigid_rotation(umat, corate, orientation):
    """An orthotropic body stretched along x, then rotated rigidly by 90 degrees: the stress
    just rotates. With lab-fixed axes ELORT used to read the rotated stretch against E_y
    (tau_22 = 664 instead of 1510, 56 %)."""
    F1 = np.diag([np.exp(0.02), 1., 1.])
    R = _rz(np.pi / 2)
    steps = [StepMeca(control="F", value=F1.ravel().tolist(), time=1., ninc=50,
                      Dn_init=1., Dn_mini=1e-4),
             StepMeca(control="F", value=(R @ F1).ravel().tolist(), time=1., ninc=100,
                      Dn_init=1., Dn_mini=1e-4)]
    r = sim.solver.solve(Block(steps=steps, control_type="F"), umat, _ORTHO, 1, T_init=290.,
                         corate=corate, orientation=orientation)
    i1 = np.argmin(np.abs(r["Time"] - 1.0))
    tau1, tau2 = _v2t(r["Kirchhoff"][:, i1]), _v2t(r["Kirchhoff"][:, -1])
    np.testing.assert_allclose(tau2, R @ tau1 @ R.T, atol=1e-10 * np.abs(tau1).max())


def test_modul_fibres_follow_the_material_like_standalone_holza():
    """A purely elastic MODUL-HOLZA material in simple shear must converge to the standalone
    HOLZA kernel, which pushes its fibres with F exactly. With lab-fixed axes the MODUL fibres
    stayed put: 190 % off at gamma = 0.5 and a factor 4e4 at gamma = 1."""
    from simcoon.modular import HolzapfelElasticity, ModularMaterial
    el = HolzapfelElasticity(C10=0.0354, k1=0.0107, k2=7.48, kappa_d=0.0,
                             fibres=sim.Rotation.from_euler("zxz", [[0.0, 0.0, 30.0]], degrees=True),
                             kappa=1000.0)
    mat = ModularMaterial(elasticity=el)

    def shear(name, props, nstatev, ninc):
        st = StepMeca(control="F", value=[1., 0.5, 0., 0., 1., 0., 0., 0., 1.], time=1., ninc=ninc,
                      Dn_init=1., Dn_mini=1e-5)
        return sim.solver.solve(Block(steps=[st], control_type="F"), name,
                                np.asarray(props, float), nstatev, T_init=290., corate=3,
                                record_tangent=False)["Kirchhoff"][:, -1]

    err = []
    for ninc in (100, 400):
        ref = shear("HOLZA", el.potential_params(), 1, ninc)
        err.append(np.abs(shear(mat.umat_name, mat.props, mat.nstatev, ninc) - ref).max()
                   / np.abs(ref).max())
    assert err[1] < 1e-4
    assert err[1] < 0.35 * err[0], "first-order convergence to the F-pushed fibres"


# ----- mechanical work: Wm is the true work per reference volume --------------------------------

@pytest.mark.parametrize("corate", [0, 2, 3, 5])
def test_wm_is_the_true_work_for_a_non_coaxial_state(corate):
    """Orthotropic ELORT in simple shear: tau is not coaxial with V, so tau : (corate strain rate)
    differs from the stress power tau : D. Wm must still converge to int P : dF."""
    props = np.array([150000., 10000., 10000., 0.3, 0.3, 0.45, 5000., 5000., 3500., 0., 0., 0.])
    rel = []
    for ninc in (100, 400):
        st = StepMeca(control="F", value=[1., 1., 0., 0., 1., 0., 0., 0., 1.], ninc=ninc,
                      Dn_init=1., Dn_mini=1e-5)
        r = sim.solver.solve(Block(steps=[st], control_type="F"), "ELORT", props, 1,
                             T_init=290., corate=corate, record_tangent=False)
        F = r["F"]
        P = [_v2t(r["Kirchhoff"][:, k]) @ np.linalg.inv(F[:, :, k]).T for k in range(len(r))]
        W = sum(0.5 * np.sum((P[k - 1] + P[k]) * (F[:, :, k] - F[:, :, k - 1]))
                for k in range(1, len(r)))
        rel.append(abs(r["Wm"][0, -1] - W) / W)
    assert rel[1] < 1e-3
    assert rel[1] < 0.3 * rel[0] or rel[1] < 1e-5
