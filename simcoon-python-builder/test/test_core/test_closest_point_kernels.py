"""tangent_mode = 3 on the dedicated plasticity kernels (EPICP, EPCHA): the closest-point
projection of return_mapping.hpp replaces their cutting-plane loops, with the exact operator."""
import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import StepMeca, solve

_UNIAXIAL = ["strain"] + ["stress"] * 5
EPICP_LIN = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 1.0]     # E nu alpha sigmaY k m (linear hardening)
EPICP_POW = [70000.0, 0.3, 1.0e-5, 300.0, 1000.0, 0.3]     # power law: singular slope at onset
EPCHA = [70000.0, 0.3, 1.0e-5, 300.0, 100.0, 20.0, 30000.0, 150.0, 5000.0, 10.0]  # + Q b C1 D1 C2 D2
col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
I3 = np.eye(3).reshape(3, 3, 1).copy(order="F")


def _run(name, props, nstatev, mode, ninc=10):
    step = StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=ninc, time=0.5)
    r = solve(step, name, props, nstatev, T_init=290.0, tangent_mode=mode)
    assert r.status == 0
    return r


def _umat(name, props, e, s, sv, De, mode):
    return sim.umat(name, col(e), col(De), I3, I3, col(s), I3, col(props), col(sv), 0.5, 0.3,
                    np.zeros((4, 1), order="F"), temp=np.full(1, 290.0), n_threads=1,
                    tangent_mode=mode, start=False)


def _tangent_vs_fd(name, props, e, s, sv, mode, h=1e-8):
    De0 = np.array([4e-4, 1e-4, -1e-4, 1.5e-4, 0.0, 0.0])
    Lt = _umat(name, props, e, s, sv, De0, mode)[3][:, :, 0]
    fd = np.column_stack([(_umat(name, props, e, s, sv, De0 + h * u, mode)[0][:, 0]
                           - _umat(name, props, e, s, sv, De0 - h * u, mode)[0][:, 0]) / (2 * h)
                          for u in np.eye(6)])
    return Lt, np.linalg.norm(Lt - fd) / np.linalg.norm(fd)


def test_epicp_closest_point_equals_cutting_plane_for_linear_hardening():
    """J2 + linear isotropic hardening: radial return, the two integrators solve the same
    equation and agree to the local tolerance."""
    r2, r3 = _run("EPICP", EPICP_LIN, 8, 2), _run("EPICP", EPICP_LIN, 8, 3)
    np.testing.assert_allclose(r3["Stress"], r2["Stress"], rtol=1e-7, atol=1e-7)
    np.testing.assert_allclose(r3["Statev"][1], r2["Statev"][1], rtol=1e-7, atol=1e-10)


def test_epicp_closest_point_stays_on_the_yield_surface_at_a_singular_onset():
    """Power-law hardening k p^m with m < 1 has an infinite slope at p = 0. The cutting-plane
    loop overshoots there and commits its unconverged-at-maxiter state (Phi tens of MPa
    inside the surface at the onset increment); the closest-point solve converges on the
    surface at every plastic increment."""
    r3 = _run("EPICP", EPICP_POW, 8, 3)
    k, m, sigma_Y = EPICP_POW[4], EPICP_POW[5], EPICP_POW[3]
    for i in range(len(r3)):
        p = r3["Statev"][1, i]
        if p > 1e-12:
            phi = sim.Mises_stress(r3["Stress"][:, i]) - k * p ** m - sigma_Y
            assert abs(phi) < 1e-6 * sigma_Y, (i, phi)


@pytest.mark.parametrize("name, props, nstatev, symmetric", [
    ("EPICP", EPICP_LIN, 8, True),
    ("EPCHA", EPCHA, 33, False),      # Armstrong-Frederick recovery: exact operator not symmetric
])
def test_dedicated_kernel_closest_point_tangent_is_exact(name, props, nstatev, symmetric):
    r = _run(name, props, nstatev, 3)
    e, s, sv = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1]
    Lt3, err3 = _tangent_vs_fd(name, props, e, s, sv, 3)
    assert err3 < 1e-8, err3
    sym = np.linalg.norm(Lt3 - Lt3.T) / np.linalg.norm(Lt3)
    assert (sym < 1e-12) if symmetric else (sym > 1e-7), sym
    if name == "EPCHA":
        _, err2 = _tangent_vs_fd(name, props, e, s, sv, 2)
        assert err2 > 1e-5   # the cutting-plane operator is approximate with kinematic hardening


def test_epcha_closest_point_cyclic_path_runs_at_full_increments():
    """Cyclic non-proportional loading of the Chaboche kernel under mode 3: no step cut."""
    steps = [StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10, time=1.0),
             StepMeca(control=["strain", "stress", "stress", "strain", "stress", "stress"],
                      value=[-0.01, 0, 0, 0.012, 0, 0], ninc=10, time=1.0),
             StepMeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10, time=1.0)]
    r = solve(steps, "EPCHA", EPCHA, 33, T_init=290.0, tangent_mode=3)
    assert r.status == 0 and len(r) == 30
    np.testing.assert_allclose(np.diff(r["Time"]), 0.1)
    assert r["Stress"][0, 9] > 300.0 and r["Stress"][0, 19] < -300.0


EPICP_T = [1e-9, 1.0, 70000.0, 0.3, 1e-5, 300.0, 1000.0, 1.0]           # rho c_p + EPICP props
EPKCP_T = [1e-9, 1.0, 70000.0, 0.3, 1e-5, 300.0, 1000.0, 1.0, 20000.0]  # + kX (Prager)


@pytest.mark.parametrize("name, props, nstatev", [("EPICP", EPICP_T, 8), ("EPKCP", EPKCP_T, 14)])
def test_thermomechanical_kernel_closest_point_tangent_is_exact(name, props, nstatev):
    """The thermomechanical twins take the same closest-point branch: dSdE is the exact
    Jacobian (the Prager backstress makes mode 2 approximate on EPKCP)."""
    from simcoon.solver import StepThermomeca
    step = StepThermomeca(control=_UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10, T_final=300.0)
    r = solve(step, name, props, nstatev, T_init=290.0, tangent_mode=3)
    assert r.status == 0
    e, s, sv, Tn = r["Strain"][:, -1], r["Stress"][:, -1], r["Statev"][:, -1], r["Temp"][-1]
    De0 = np.array([4e-4, 1e-4, -1e-4, 1.5e-4, 0.0, 0.0])

    def call(d, mode):
        out = sim.umat_T(name, col(e), col(d), col(s), I3, col(props), col(sv), 0.5, 0.3,
                         np.zeros((4, 1), order="F"), np.zeros((3, 1), order="F"),
                         np.array([Tn]), np.array([0.0]), tangent_mode=mode)
        return out[0][:, 0], out[5][:, :, 0]

    h = 1e-8
    fd = np.column_stack([(call(De0 + h * u, 3)[0] - call(De0 - h * u, 3)[0]) / (2 * h) for u in np.eye(6)])
    assert np.linalg.norm(call(De0, 3)[1] - fd) < 1e-8 * np.linalg.norm(fd)
    if name == "EPKCP":
        fd2 = np.column_stack([(call(De0 + h * u, 2)[0] - call(De0 - h * u, 2)[0]) / (2 * h) for u in np.eye(6)])
        assert np.linalg.norm(call(De0, 2)[1] - fd2) > 1e-5 * np.linalg.norm(fd2)


SMADI = [0, 70000.0, 70000.0, 0.3, 0.3, 1e-6, 1e-6, 0.0, 0.05, 0.021, 0.0, 6.0, 5.0,
         293.15, 273.15, 313.15, 333.15, 0.2, 0.2, 0.2, 0.2, 300.0, 1.4, 2.0, 1e-6, 1e-3, 1.0, 1e8]
T_SUPERELASTIC = 353.15


def test_sma_unified_T_closest_point_superelastic_loop_and_exact_tangent():
    """SMADI on the closest-point helper (frozen mixture stiffness, effective flows): the
    superelastic loop runs at full increments in mode 3 and the operator is the exact
    Jacobian in the transformation regime (mode 2 is clamped to the continuum operator on the
    SMA kernels)."""
    steps = [StepMeca(control=_UNIAXIAL, value=[0.04, 0, 0, 0, 0, 0], time=1.0, ninc=100, Dn_mini=0.01),
             StepMeca(control=_UNIAXIAL, value=[0.0] * 6, time=1.0, ninc=100, Dn_mini=0.01)]
    r2 = solve(steps, "SMADI", SMADI, 30, T_init=T_SUPERELASTIC, corate=3, tangent_mode=2)
    r3 = solve(steps, "SMADI", SMADI, 30, T_init=T_SUPERELASTIC, corate=3, tangent_mode=3)
    assert r3.status == 0 and len(r3) == 200
    np.testing.assert_allclose(np.diff(r3["Time"]), 0.01)
    np.testing.assert_allclose(r3["Stress"][0], r2["Stress"][0], atol=0.5)   # same loop, O(h^2) apart
    k = 60   # forward transformation
    e, s, sv = r3["Strain"][:, k], r3["Stress"][:, k], r3["Statev"][:, k]
    De0 = np.array([4e-4, -1e-4, -1e-4, 1.5e-4, 0.0, 0.0])

    def call(d, mode):
        return sim.umat("SMADI", col(e), col(d), I3, I3, col(s), I3, col(SMADI), col(sv), 0.5, 0.01,
                        np.zeros((4, 1), order="F"), temp=np.full(1, T_SUPERELASTIC), n_threads=1,
                        tangent_mode=mode, start=False)

    h = 1e-7
    Lt = call(De0, 3)[3][:, :, 0]
    fd = np.column_stack([(call(De0 + h * u, 3)[0][:, 0] - call(De0 - h * u, 3)[0][:, 0]) / (2 * h)
                          for u in np.eye(6)])
    assert np.linalg.norm(Lt - fd) < 1e-5 * np.linalg.norm(fd)


DFA_ISO = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 0.0]


@pytest.mark.parametrize("name, extra", [("SMRDI", []), ("SMRAI", DFA_ISO)])
def test_sma_unified_TR_closest_point_loop_and_tangent(name, extra):
    """unified_TR (transformation + reorientation, 3 rows) on the helper: the superelastic loop
    closes in mode 3 as in mode 2, and the operator is the exact Jacobian in the forward
    transformation (mode 2 is clamped to the continuum operator on the SMA kernels)."""
    props = np.asarray(SMADI + extra + [100.0, 5000.0, 0.05, 1e-6, 1e-3, 1.0, 1e8], dtype=float)
    steps = [StepMeca(control=_UNIAXIAL, value=[0.04, 0, 0, 0, 0, 0], time=1.0, ninc=200, Dn_mini=0.01),
             StepMeca(control=["stress"] * 6, value=[0.0] * 6, time=1.0, ninc=200, Dn_mini=0.01)]
    r2 = solve(steps, name, props, 30, T_init=T_SUPERELASTIC, corate=3, tangent_mode=2)
    r3 = solve(steps, name, props, 30, T_init=T_SUPERELASTIC, corate=3, tangent_mode=3)
    assert r2.status == 0 and r3.status == 0 and len(r3) == 400
    np.testing.assert_allclose(r3["Stress"][0], r2["Stress"][0], atol=0.5)
    assert abs(r3["Strain"][0, -1]) < 1e-3   # the loop closes (full reverse transformation)
    k = 120
    e, s, sv = r3["Strain"][:, k], r3["Stress"][:, k], r3["Statev"][:, k]
    De0 = np.array([4e-4, -1e-4, -1e-4, 1.5e-4, 0.0, 0.0])

    def call(d):
        return sim.umat(name, col(e), col(d), I3, I3, col(s), I3, col(props), col(sv), 0.5, 0.005,
                        np.zeros((4, 1), order="F"), temp=np.full(1, T_SUPERELASTIC), n_threads=1,
                        tangent_mode=3, start=False)

    h = 1e-7
    Lt = call(De0)[3][:, :, 0]
    fd = np.column_stack([(call(De0 + h * u)[0][:, 0] - call(De0 - h * u)[0][:, 0]) / (2 * h) for u in np.eye(6)])
    assert np.linalg.norm(Lt - fd) < 1e-5 * np.linalg.norm(fd)
