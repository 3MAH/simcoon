"""The total formulation (strain space), ``formulation="total"``.

The law works on ln U and on its dual stress T, where the rate formulation works on ln V and the
Kirchhoff stress. A log-strain box kernel is fed the material Hencky strain H = 1/2 ln C = ln U, a state function of
F in reference components, and its output is read as the dual stress T (T:dH = S:dE = tau:D), mapped
once to the Kirchhoff stress (Miehe, Apel and Lambrecht, CMAME 191, 2002). The box, its strain and
its internal variables are the ones the rate formulation works with under corate 3 (log_R); what
differs is

* the strain increment: a difference of a state function, exact at any step, against an integral;
* the reading of the box output: T against R^T tau R. They coincide on a state coaxial with U and
  differ otherwise, which moves the reported stress under plastic shear and makes an anisotropic
  elastic law a stored-energy one.
"""

import numpy as np
import pytest
from scipy.linalg import expm, logm

import simcoon as sim
from simcoon.solver import Block, StepMeca

EPICP = np.array([210000., 0.3, 0., 250., 1., 1.])      # sigma_Y = 250 MPa, linear hardening of 1 MPa
ELISO = np.array([210000., 0.3, 0.])
ORTHO = np.array([4500., 2300., 2700., 0.06, 0.08, 0.3, 2200., 2100., 2400., 0., 0., 0.])


def _v2t(v):
    return np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])


def _t2v(m):
    return np.array([m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[0, 2], m[1, 2]])


def _t2v_strain(m):
    return np.array([m[0, 0], m[1, 1], m[2, 2], 2 * m[0, 1], 2 * m[0, 2], 2 * m[1, 2]])


def _rz(th):
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


def _shear(gamma):
    return np.array([[1., gamma, 0.], [0., 1., 0.], [0., 0., 1.]])


def _polar(F):
    c, N = np.linalg.eigh(F.T @ F)
    U = N @ np.diag(np.sqrt(c)) @ N.T
    return F @ np.linalg.inv(U), U


def _lnU(F):
    c, N = np.linalg.eigh(F.T @ F)
    return N @ np.diag(0.5 * np.log(c)) @ N.T


def _lnV(F):
    b, N = np.linalg.eigh(F @ F.T)
    return N @ np.diag(0.5 * np.log(b)) @ N.T


def _mises(s):
    d = s - np.trace(s) / 3 * np.eye(3)
    return np.sqrt(1.5 * np.sum(d * d))


def _dual_stress(tau, F):
    """T from tau: the inverse of tau = A_R^T (R T R^T), reference components."""
    R, _ = _polar(F)
    return R.T @ _v2t(np.linalg.solve(np.asarray(sim.A_R(F)).T, _t2v(tau))) @ R


def _formulation(corate):
    """solve() keywords: an objective rate, or the total formulation."""
    return {"formulation": "total"} if corate == "total" else {"corate": corate}


def _solve(umat, props, nstatev, targets, ninc, corate, control_type="F"):
    steps = [StepMeca(control="F", value=np.asarray(F).ravel().tolist(), time=1., ninc=ninc,
                      Dn_init=1., Dn_mini=1e-6) for F in targets]
    return sim.solver.solve(Block(steps=steps, control_type=control_type), umat,
                            np.asarray(props, float), nstatev, T_init=290., **_formulation(corate))


def _umat6(umat, props, F0, F1, sigma, statev, Wm, start=False, **kw):
    """One sim.umat call under the total formulation: etot = H(F0), Detot = H(F1) - H(F0), DR = I."""
    col = lambda a: np.asfortranarray(np.asarray(a, float).reshape(-1, 1))
    cube = lambda a: np.asfortranarray(np.asarray(a, float).reshape(3, 3, 1))
    H0, H1 = _t2v_strain(_lnU(F0)), _t2v_strain(_lnU(F1))
    s, sv, wm, Lt = sim.umat(umat, col(H0), col(H1 - H0), cube(F0), cube(F1), col(sigma), cube(np.eye(3)),
                             col(props), col(statev), 0. if start else 1., 0. if start else 1., col(Wm),
                             np.array([290.]), n_threads=1, formulation="total", start=start, **kw)
    return s[:, 0], np.asarray(sv)[:, 0], np.asarray(wm)[:, 0], Lt[:, :, 0]


def test_total_is_a_formulation_not_a_corate():
    assert sim.solver.FORMULATIONS == ("rate", "total")
    assert sim.corate_type_total not in sim.solver.CORATE_TYPES.values()
    with pytest.raises(ValueError, match="corate"):
        _solve("ELISO", ELISO, 1, [_shear(0.1)], 1, sim.corate_type_total)
    step = StepMeca(control="F", value=_shear(0.1).ravel().tolist(), time=1., ninc=1)
    with pytest.raises(ValueError, match="no objective rate"):
        sim.solver.solve(Block(steps=[step], control_type="F"), "ELISO", ELISO, 1, corate=3, formulation="total")
    with pytest.raises(ValueError, match="formulation"):
        sim.solver.solve(Block(steps=[step], control_type="F"), "ELISO", ELISO, 1, formulation="lagrangian")


def test_elastic_branch_is_exact_in_one_increment():
    """Simple shear to gamma = 1 with an isotropic elastic law: tau = L : ln V in a single
    increment, where the accumulated log_R strain carries an O(dt) residue that halves with the
    step."""
    F = _shear(1.0)
    exact = _v2t(np.asarray(sim.L_iso([210000., 0.3], "Enu")) @ _t2v_strain(_lnV(F)))
    r = _solve("ELISO", ELISO, 1, [F], 1, "total")
    np.testing.assert_allclose(_v2t(r["Kirchhoff"][:, -1]), exact, atol=1e-12 * np.abs(exact).max())
    np.testing.assert_allclose(r["Strain"][:, -1], _t2v_strain(_lnU(F)), atol=1e-14)

    err = [np.linalg.norm(_v2t(_solve("ELISO", ELISO, 1, [F], n, 3)["Kirchhoff"][:, -1]) - exact) for n in (40, 80)]
    assert err[1] > 1e-5 * np.abs(exact).max()
    assert err[0] / err[1] == pytest.approx(2.0, abs=0.1)


def test_plastic_simple_shear_reads_the_box_as_the_dual_stress():
    """J2 plasticity in simple shear to gamma = 1. The criterion holds on T, so the reported
    stress sits inside it; the log_R route holds it on tau. Reference values of the strain-space
    framework (MFront LogarithmicStrainPlasticity, same parameters): s12 = 140.9, s22 = -15.7,
    von Mises 245.5 MPa. log_R gives 144.4 and -8.4: 10.7 MPa of difference that no step removes,
    while the two box histories are the same."""
    F = _shear(1.0)
    r6 = _solve("EPICP", EPICP, 8, [F], 400, "total")
    tau = _v2t(r6["Kirchhoff"][:, -1])
    assert tau[0, 1] == pytest.approx(140.9, abs=0.1)
    assert tau[1, 1] == pytest.approx(-15.7, abs=0.1)
    assert _mises(tau) == pytest.approx(245.5, abs=0.1)
    p = r6["Statev"][1, -1]
    T = _dual_stress(tau, F)
    assert _mises(T) == pytest.approx(250. + 1. * p, abs=1e-6)       # on the yield surface
    np.testing.assert_allclose(_v2t(r6["PKII"][:, -1]), np.linalg.inv(F) @ tau @ np.linalg.inv(F).T, atol=1e-9)

    r3 = _solve("EPICP", EPICP, 8, [F], 400, 3)
    tau3 = _v2t(r3["Kirchhoff"][:, -1])
    assert _mises(tau3) == pytest.approx(250. + r3["Statev"][1, -1], abs=1e-6)
    assert np.linalg.norm(tau3 - tau) > 5.
    R, _ = _polar(F)
    np.testing.assert_allclose(R.T @ tau3 @ R, T, atol=0.1)          # one box, two readings
    # ... and the reading is the A_R map
    np.testing.assert_allclose(_v2t(np.asarray(sim.A_R(F)).T @ _t2v(tau3)), tau, atol=0.1)


def test_orthotropic_closed_cycle_stores_no_net_work():
    """Stretch, shear and rotation along a closed path of F with an orthotropic elastic law:
    T = L : H derives from a stored energy, so the work returns to zero at any step. Read as
    tau = R (L : H) R^T (log_R) the same law is Cauchy-elastic and the cycle stays open."""
    Fa = _shear(0.5) @ np.diag([1.2, 0.9, 1.0])
    Fb = _rz(0.8) @ np.diag([0.9, 1.25, 1.05]) @ _shear(-0.3)
    r6 = _solve("ELORT", ORTHO, 1, [Fa, Fb, np.eye(3)], 20, "total")
    peak = np.abs(r6["Wm"][0]).max()
    assert peak > 100.
    assert abs(r6["Wm"][0, -1]) < 1e-10 * peak
    assert np.abs(r6["Kirchhoff"][:, -1]).max() < 1e-9

    r3 = _solve("ELORT", ORTHO, 1, [Fa, Fb, np.eye(3)], 80, 3)
    assert abs(r3["Wm"][0, -1]) > 1e-4 * peak


def test_work_is_the_stress_power():
    """Wm accumulates T:dH, which is tau:D exactly: no work correction on this route."""
    F = _rz(0.4) @ np.diag([1.15, 0.95, 1.0]) @ _shear(0.4)
    r = _solve("EPICP", EPICP, 8, [F], 200, "total")
    Fs = np.concatenate([np.eye(3)[:, :, None], r["F"]], axis=2)
    taus = np.concatenate([np.zeros((6, 1)), r["Kirchhoff"]], axis=1)
    power = 0.
    for i in range(Fs.shape[2] - 1):
        F0, F1 = Fs[:, :, i], Fs[:, :, i + 1]
        LDt = 2. * (F1 - F0) @ np.linalg.inv(F1 + F0)
        power += 0.5 * np.sum((_v2t(taus[:, i]) + _v2t(taus[:, i + 1])) * (0.5 * (LDt + LDt.T)))
    assert r["Wm"][0, -1] == pytest.approx(power, rel=2e-4)


def test_a_rate_law_is_refused():
    """HYPOO integrates from its start stress, which has no meaning in the strain space."""
    with pytest.raises(RuntimeError, match="total formulation"):
        _solve("HYPOO", ORTHO, 1, [_shear(0.1)], 5, "total")


@pytest.mark.parametrize("umat, props", [("NEOHC", [1000., 10000.]), ("SNTVE", [100000., 0.3, 0.])])
def test_a_law_built_on_F_is_the_same_law_in_both_formulations(umat, props):
    """A hyperelastic law built on F is a state function of F: stretch, shear, rotation, and back.
    Its stress does not depend on the formulation, its strain output is ln U, and its work, formed
    on the dual pair (T, ln U), returns to zero over the closed path."""
    Fa = _rz(0.6) @ np.diag([1.3, 0.85, 1.0]) @ _shear(0.4)
    path = [Fa, _rz(1.2) @ Fa, np.eye(3)]
    rt = _solve(umat, props, 1, path, 40, "total")
    rr = _solve(umat, props, 1, path, 40, 3)
    scale = np.abs(rr["Kirchhoff"]).max()
    np.testing.assert_allclose(rt["Kirchhoff"], rr["Kirchhoff"], atol=1e-9 * scale)
    np.testing.assert_allclose(rt["Strain"][:, 39], _t2v_strain(_lnU(Fa)), atol=1e-12)
    peak = np.abs(rt["Wm"][0]).max()
    assert peak > 1.
    # the work on the dual pair (T, ln U) is the stress power the rate formulation integrates
    assert rt["Wm"][0, 39] == pytest.approx(rr["Wm"][0, 39], rel=2e-3)
    # a rigid rotation does no work ...
    assert abs(rt["Wm"][0, 79] - rt["Wm"][0, 39]) < 1e-10 * peak
    # ... and the cycle closes at the accuracy of the trapezoidal rule
    assert abs(rt["Wm"][0, -1]) < 2e-3 * peak


def test_a_law_built_on_F_converges_under_mixed_control():
    """Uniaxial tension under Green-Lagrange / PKII control: the Newton loop runs on the dS/dE
    rebuilt from the dT/dH box of the law, and lands on the solution of the rate formulation."""
    st = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.3, 0, 0, 0, 0, 0], time=1., ninc=10,
                  Dn_init=1., Dn_mini=1e-3)
    out = {c: sim.solver.solve(Block(steps=[st], control_type="green_lagrange"), "NEOHC", [1000., 10000.], 1,
                               T_init=290., **_formulation(c)) for c in (3, "total")}
    assert len(out["total"]) == len(out[3]) == 10                 # no step cut
    np.testing.assert_allclose(out["total"]["Kirchhoff"][:, -1], out[3]["Kirchhoff"][:, -1], atol=1e-6)


def test_small_strain_control_is_refused():
    """A small-strain block accumulates a strain that F does not carry, and the total
    formulation re-reads its strain from F."""
    st = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.01, 0, 0, 0, 0, 0], time=1., ninc=2)
    with pytest.raises(RuntimeError, match="control_type 1"):
        sim.solver.solve(Block(steps=[st]), "EPICP", EPICP, 8, T_init=290., formulation="total")


def test_logarithmic_control_is_refused():
    st = StepMeca(control=["strain"] + ["stress"] * 5, value=[0.1, 0, 0, 0, 0, 0], time=1., ninc=10)
    with pytest.raises(RuntimeError, match="control_type 3"):
        sim.solver.solve(Block(steps=[st], control_type="logarithmic"), "EPICP", EPICP, 8, T_init=290.,
                         formulation="total")


@pytest.mark.parametrize("control_type, target", [("green_lagrange", 0.15), ("biot", 1.15)])
def test_mixed_control_converges_like_log_R_on_a_coaxial_path(control_type, target):
    """Uniaxial tension under mixed control: the Newton loop runs on the dS/dE of the total
    formulation. The state is coaxial, so both formulations describe the same law and agree to
    the log_R quadrature."""
    st = StepMeca(control=["strain"] + ["stress"] * 5, value=[target, 0, 0, 0, 0, 0], time=1., ninc=50,
                  Dn_init=1., Dn_mini=1e-3)
    out = {c: sim.solver.solve(Block(steps=[st], control_type=control_type), "EPICP", EPICP, 8, T_init=290.,
                               **_formulation(c)) for c in (3, "total")}
    assert len(out["total"]) == len(out[3]) == 50                 # no step cut
    assert np.abs(out["total"]["Kirchhoff"][1:, -1]).max() < 1e-6
    assert out["total"]["Kirchhoff"][0, -1] == pytest.approx(out[3]["Kirchhoff"][0, -1], rel=1e-6)
    assert out["total"]["Wm"][0, -1] == pytest.approx(out[3]["Wm"][0, -1], rel=1e-5)


def test_umat_reproduces_the_solver():
    """sim.umat under the total formulation: Cauchy stress in and out, the dual stress never
    crosses the boundary."""
    n = 25
    r = _solve("EPICP", EPICP, 8, [_shear(1.0)], n, "total")
    s, sv, wm, _ = _umat6("EPICP", EPICP, np.eye(3), np.eye(3), np.zeros(6), np.zeros(8), np.zeros(4), start=True)
    F0 = np.eye(3)
    for g in np.linspace(0, 1, n + 1)[1:]:
        s, sv, wm, _ = _umat6("EPICP", EPICP, F0, _shear(g), s, sv, wm)
        F0 = _shear(g)
    np.testing.assert_allclose(s, r["Stress"][:, -1], atol=1e-8)
    np.testing.assert_allclose(wm, r["Wm"][:, -1], atol=1e-9)
    np.testing.assert_allclose(sv[1], r["Statev"][1, -1], atol=1e-12)


def test_umat_refuses_a_rate_law():
    with pytest.raises(ValueError, match="total formulation"):
        _umat6("HYPOO", ORTHO, np.eye(3), _shear(0.1), np.zeros(6), np.zeros(1), np.zeros(4))


def test_umat_takes_no_corate_under_the_total_formulation():
    with pytest.raises(ValueError, match="no objective rate"):
        _umat6("EPICP", EPICP, np.eye(3), _shear(0.1), np.zeros(6), np.zeros(8), np.zeros(4), corate=3)


def test_umat_serves_a_law_built_on_F_like_the_solver():
    """NEOHC through sim.umat under the total formulation: the stress and the work of the solver."""
    props, n = [1000., 10000.], 20
    Fa = _rz(0.6) @ np.diag([1.3, 0.85, 1.0]) @ _shear(0.4)
    r = _solve("NEOHC", props, 1, [Fa], n, "total")
    s, sv, wm, _ = _umat6("NEOHC", props, np.eye(3), np.eye(3), np.zeros(6), np.zeros(1), np.zeros(4), start=True)
    for k in range(n):
        F0 = np.eye(3) if k == 0 else r["F"][:, :, k - 1]
        s, sv, wm, _ = _umat6("NEOHC", props, F0, r["F"][:, :, k], s, sv, wm)
    np.testing.assert_allclose(s, r["Stress"][:, -1], atol=1e-9 * np.abs(r["Stress"]).max())
    np.testing.assert_allclose(wm, r["Wm"][:, -1], atol=1e-9 * np.abs(r["Wm"]).max())


def test_umat_needs_the_deformation_gradients():
    col = lambda a: np.asfortranarray(np.asarray(a, float).reshape(-1, 1))
    empty = np.empty((3, 3, 0), order="F")
    with pytest.raises(ValueError, match="F0 and F1"):
        sim.umat("EPICP", col(np.zeros(6)), col(np.zeros(6)), empty, empty, col(np.zeros(6)),
                 np.asfortranarray(np.eye(3).reshape(3, 3, 1)), col(EPICP), col(np.zeros(8)), 0., 1.,
                 col(np.zeros(4)), n_threads=1, formulation="total")


def test_material_tangent_matches_finite_difference_with_plastic_history():
    """dS/dE = P^T (dT/dH) P + T:L at a plastified, sheared and rotated state, against a finite
    difference of S(E) through the kernel. The geometric term carries a non-coaxial T here, and
    the box is the algorithmic one: the Jacobian of the discrete update."""
    s, sv, wm, _ = _umat6("EPICP", EPICP, np.eye(3), np.eye(3), np.zeros(6), np.zeros(8), np.zeros(4), start=True)
    Fa = _rz(0.3) @ np.diag([1.04, 0.99, 0.98]) @ _shear(0.08)
    s, sv, wm, _ = _umat6("EPICP", EPICP, np.eye(3), Fa, s, sv, wm)
    assert sv[1] > 1e-2
    Fb = _rz(0.5) @ np.diag([1.07, 0.97, 0.99]) @ _shear(0.14)
    Rb, Ub = _polar(Fb)

    def S_of(C):
        c, N = np.linalg.eigh(C)
        F = Rb @ N @ np.diag(np.sqrt(c)) @ N.T
        sig, svb, _, Lt = _umat6("EPICP", EPICP, Fa, F, s, sv, wm, tangent_mode=2, tangent_output="material")
        Fi = np.linalg.inv(F)
        return _t2v(np.linalg.det(F) * Fi @ _v2t(sig) @ Fi.T), svb, Lt

    S, svb, DSDE = S_of(Ub @ Ub)
    assert svb[1] > sv[1] + 1e-3                                    # the increment is plastic
    h = 1e-7
    fd = np.zeros((6, 6))
    for j in range(6):
        dE = np.zeros(6)
        dE[j] = 1.
        dE = np.array([[dE[0], dE[3] / 2, dE[4] / 2], [dE[3] / 2, dE[1], dE[5] / 2], [dE[4] / 2, dE[5] / 2, dE[2]]])
        fd[:, j] = (S_of(Ub @ Ub + 2 * h * dE)[0] - S_of(Ub @ Ub - 2 * h * dE)[0]) / (2 * h)
    scale = np.abs(fd).max()
    np.testing.assert_allclose(DSDE, fd, atol=1e-6 * scale)


# ---------------------------------------------------------------------------------------------
# The measures of one formulation read in a run of the other: ln U <-> ln V, T <-> tau
# ---------------------------------------------------------------------------------------------


def _rotate(R, v, strain=False):
    """R m R^T on a Voigt vector."""
    m = R @ (_v2t(v) if not strain else _v2t(v * np.array([1, 1, 1, .5, .5, .5]))) @ R.T
    return _t2v_strain(m) if strain else _t2v(m)


@pytest.mark.parametrize("corate", [2, 3, 4, "total"])
def test_the_measures_of_the_other_formulation_are_derived_from_F(corate):
    """LogStrainU, LogStrainV and HenckyStress are changes of measure, computed from F and the
    Kirchhoff stress: ln V = R (ln U) R^T, and tau = Hencky2Kirchoff(T)."""
    F = _rz(0.4) @ np.diag([1.15, 0.95, 1.0]) @ _shear(0.4)
    r = _solve("EPICP", EPICP, 8, [F], 20, corate)
    assert r["LogStrainV"] is r["LogStrain"]
    R, _ = _polar(F)
    np.testing.assert_allclose(r["LogStrainU"][:, -1], _t2v_strain(_lnU(F)), atol=1e-12)
    np.testing.assert_allclose(_rotate(R, r["LogStrainU"][:, -1], strain=True), r["LogStrainV"][:, -1], atol=1e-12)
    T = r["HenckyStress"]
    assert T.shape == r["Kirchhoff"].shape
    back = sim.stress_convert(np.asfortranarray(T), np.asfortranarray(r["F"]), "Hencky2Kirchoff")
    np.testing.assert_allclose(back, r["Kirchhoff"], atol=1e-9 * np.abs(r["Kirchhoff"]).max())
    np.testing.assert_allclose(_v2t(T[:, -1]), _dual_stress(_v2t(r["Kirchhoff"][:, -1]), F), atol=1e-8)
    if corate == "total":
        np.testing.assert_allclose(r["Strain"], r["LogStrainU"], atol=1e-12)


@pytest.mark.parametrize("corate", [0, 1, 5])
def test_the_material_measures_are_not_offered_when_the_strain_is_not_a_state_function(corate):
    r = _solve("EPICP", EPICP, 8, [_shear(0.3)], 5, corate)
    assert "LogStrainU" not in r and "HenckyStress" not in r
    assert "LogStrainV" in r
    with pytest.raises(KeyError):
        r["HenckyStress"]


def test_one_box_read_in_the_two_formulations():
    """Plastic simple shear under corate 3 and under the total formulation, on the same measures:
    the box output is the same (R^T tau_rate R = T_total), and the two Kirchhoff stresses differ
    by the A_R map, tau_total = Hencky2Kirchoff(R^T tau_rate R). Neither closes with the step
    beyond the quadrature of the log_R strain."""
    F = _shear(1.0)
    R, _ = _polar(F)
    gap = []
    for n in (100, 200):
        rr, rt = _solve("EPICP", EPICP, 8, [F], n, 3), _solve("EPICP", EPICP, 8, [F], n, "total")
        box_rate = _rotate(R.T, rr["Kirchhoff"][:, -1])
        gap.append(np.abs(box_rate - rt["HenckyStress"][:, -1]).max())
        read_as_T = np.ravel(sim.stress_convert(box_rate, F, "Hencky2Kirchoff"))
        np.testing.assert_allclose(read_as_T, rt["Kirchhoff"][:, -1], atol=2 * gap[-1] + 1e-9)
        np.testing.assert_allclose(rr["Statev"][1, -1], rt["Statev"][1, -1], rtol=1e-3)
    assert gap[1] < 0.6 * gap[0] < 0.1                             # first order in the step
    assert np.abs(rr["Kirchhoff"][:, -1] - rt["Kirchhoff"][:, -1]).max() > 5.   # the reading, not the step


def test_which_laws_the_total_formulation_serves():
    """Log-strain boxes and hyperelastic laws built on F; not the rate law, the modular
    material, nor a small-strain-only law."""
    for name in ("ELISO", "ELORT", "EPICP", "EPCHA", "EPJCK", "PYEXT", "NEOHC", "SNTVE", "HOLZA", "OGDEN"):
        assert sim.total_formulation_serves(name), name
    for name in ("HYPOO", "MODUL", "SMADI", "ZENER", "UMEXT", "no such law"):
        assert not sim.total_formulation_serves(name), name
