"""Finite-difference the box tangent of every finite-strain kernel.

Why this file exists: **a wrong tangent does not move converged solver values.** Newton
finds the same root with a poor Jacobian, so `regression_baseline.py` -- 100 combos of
umat x control type x corate -- is structurally blind to this entire class of error. It
caught neither of the two real tangent defects found in this area:

* ``saint_venant`` fed its stress to a helper that rebuilds ``tau = det(F)*sigma``, i.e.
  expects **Cauchy**, after the kernel had been made Kirchhoff-native -- squaring the J;
* ``HYPOO`` handed over ``dsigma/dD`` where the consumer reads ``d(tau_hat)/dDe`` (it now
  integrates the Kirchhoff rate, for which ``L`` is the box itself).

Both are invisible at J = 1, so every state here carries J well away from 1 and the
docstrings state the relative size a J error would have.
"""

import numpy as np
import pytest
from scipy.linalg import expm, polar

import simcoon as sim

# (name, props, nstatev) -- the finite kernels that build their own tangent
KERNELS = [
    ("SNTVE", [100000.0, 0.3, 0.0], 1),
    ("NEOHC", [1000.0, 10000.0], 1),
    ("NEOHI", [2000.0, 0.40, 0.0], 1),
    ("MOORI", [0.2588, -0.0449, 10000.0], 1),
    ("YEOHH", [0.30, -0.010, 0.0005, 1000.0], 1),
]

# Coaxial: ln V and the perturbation commute, so the box tangent's normal block IS
# d(tau)/d(eps) and no rate conversion is needed. The trace is deliberately NOT zero.
EPS0 = np.array([0.11, -0.04, -0.03, 0.0, 0.0, 0.0])


def _umat(name, props, F1, nstatev=1, corate=3):
    z6 = lambda: np.zeros((6, 1), order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    stress, sv, wm, Lt = sim.umat(
        name, z6(), z6(), eye,
        np.asarray(F1).reshape(3, 3, 1).copy(order="F"), z6(), eye,
        np.asfortranarray(np.asarray(props, dtype=float).reshape(-1, 1)),
        np.zeros((nstatev, 1), order="F"), 0.0, 1.0,
        np.zeros((4, 1), order="F"), n_threads=1, corate=corate)
    return stress[:, 0], Lt[:, :, 0]


def _F(eps):
    return expm(sim.v2t_strain(np.asarray(eps, dtype=float)))


@pytest.mark.parametrize("corate", [2, 3])
@pytest.mark.parametrize("name,props,nstatev", KERNELS, ids=[k[0] for k in KERNELS])
def test_box_tangent_matches_finite_difference(name, props, nstatev, corate):
    """Lt[:3,:3] must be the central difference of tau at a coaxial, J != 1 state.

    Corates 2 (XBM) and 3 (log_R) share the exact spectral map, so the box is
    d(tau)/d(ln V) for both and the same finite difference applies. J = 1.0408 here, so a
    stray or missing power of J shows up at ~4 % -- six orders above the tolerance.
    """
    J = float(np.exp(EPS0[:3].sum()))
    assert abs(J - 1.0) > 0.03, "J must be away from 1 or a J error would be invisible"

    def tau_of(eps):
        sigma, _ = _umat(name, props, _F(eps), nstatev, corate)
        # sim.umat returns CAUCHY at the boundary; the box tangent is Kirchhoff, no J
        return np.exp(np.sum(eps[:3])) * np.asarray(sigma).ravel()

    _, Lt = _umat(name, props, _F(EPS0), nstatev, corate)
    d = 1e-6
    for col in range(3):
        step = np.zeros(6)
        step[col] = d
        fd = (tau_of(EPS0 + step) - tau_of(EPS0 - step)) / (2.0 * d)
        np.testing.assert_allclose(
            fd[:3], Lt[:3, col], rtol=2e-6,
            atol=2e-6 * max(1.0, np.abs(Lt[:3, col]).max()),
            err_msg=f"{name}: d(tau)/d(eps) column {col} at corate {corate}")


@pytest.mark.parametrize("name,props,nstatev", KERNELS, ids=[k[0] for k in KERNELS])
def test_truesdell_box_is_the_convected_tangent(name, props, nstatev):
    """Corate 4 (Truesdell): the box is d(tau_hat)/d(De) with the Kirchhoff stress transported
    upper-convected, tau_hat = DF tau0 DF^T, and De the Almansi increment
    1/2 (I - (DF DF^T)^-1). Its limit at DF -> I is the Lie tangent, which is what the kernels
    return for corate 4; checked by central differences about the reference state F0."""
    F0 = _F(EPS0)
    sigma0, Lt = _umat(name, props, F0, nstatev, corate=4)
    tau0 = sim.v2t_stress(np.exp(EPS0[:3].sum()) * np.asarray(sigma0).ravel())
    d = 1e-6
    for col in range(3):
        dtau, de = [], []
        for sgn in (1.0, -1.0):
            DF = np.eye(3)
            DF[col, col] += sgn * d
            F1 = DF @ F0
            sigma1, _ = _umat(name, props, F1, nstatev, corate=4)
            tau1 = np.linalg.det(F1) * np.asarray(sigma1).ravel()
            dtau.append(tau1 - np.asarray(sim.t2v_stress(DF @ tau0 @ DF.T)).ravel())
            de.append(0.5 * (1.0 - 1.0 / (1.0 + sgn * d) ** 2))
        fd = (dtau[0] - dtau[1]) / (de[0] - de[1])
        np.testing.assert_allclose(
            fd[:3], Lt[:3, col], rtol=1e-5,
            atol=1e-5 * max(1.0, np.abs(Lt[:3, col]).max()),
            err_msg=f"{name}: Truesdell box column {col}")


def test_corates_2_and_3_return_the_same_box():
    """They resolve to the same exact spectral map, so the box is literally identical.

    This is what makes ``corate=3`` a contract-preserving default for ``sim.umat``: it
    returns precisely the box the finite kernels used to bake unconditionally.
    """
    for name, props, nstatev in KERNELS:
        _, Lt2 = _umat(name, props, _F(EPS0), nstatev, corate=2)
        _, Lt3 = _umat(name, props, _F(EPS0), nstatev, corate=3)
        np.testing.assert_allclose(Lt3, Lt2, rtol=1e-13, atol=1e-13,
                                   err_msg=f"{name}: corate 2 and 3 must coincide")


def test_jaumann_and_green_naghdi_differ_from_the_log_box():
    """Otherwise the corate argument would be doing nothing and the tests above vacuous.

    The threshold is deliberately loose. The Jaumann and Green-Naghdi corrections are
    **stress-proportional**, so for a near-incompressible parameter set they are tiny next
    to the largest tangent entry: MOORI here has ``kappa = 10000`` against ``C10 = 0.26``,
    so the bulk term dominates the matrix and the correction is ~3e-3 out of ~1e4. That is
    still ten orders above round-off, and a kernel that ignored the corate would return
    the log box **exactly**.
    """
    for name, props, nstatev in KERNELS:
        _, log_box = _umat(name, props, _F(EPS0), nstatev, corate=3)
        scale = max(1.0, np.abs(log_box).max())
        for corate in (0, 1):
            _, other = _umat(name, props, _F(EPS0), nstatev, corate=corate)
            assert np.abs(other - log_box).max() > 1e-9 * scale, \
                f"{name}: corate {corate} returned the log box unchanged"


HYPOO_PROPS = [70000., 60000., 50000., 0.3, 0.28, 0.25, 26000., 22000., 20000., 0., 0., 0.]


def _hypoo(De, props=HYPOO_PROPS):
    """One HYPOO increment from the reference state, with F1 consistent as exp(De).

    HYPOO is a RATE kernel: it responds to ``Detot``, not to the total strain, so any
    difference is taken along the increment.
    """
    De = np.asarray(De, dtype=float)
    F1 = expm(sim.v2t_strain(De)).reshape(3, 3, 1).copy(order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    stress, _, _, Lt = sim.umat(
        "HYPOO", np.zeros((6, 1), order="F"), np.asfortranarray(De.reshape(6, 1)),
        eye, F1, np.zeros((6, 1), order="F"), eye,
        np.asfortranarray(np.asarray(props, dtype=float).reshape(-1, 1)),
        np.zeros((1, 1), order="F"), 0.0, 1.0,
        np.zeros((4, 1), order="F"), n_threads=1, corate=3)
    return np.asarray(stress).ravel(), Lt[:, :, 0], float(np.linalg.det(F1[:, :, 0]))


def test_hypoelastic_integrates_the_kirchhoff_rate():
    """HYPOO integrates ``tau_{n+1} = tau_n + L : DEel``, so ``L`` IS its box tangent.

    ``sim.umat`` exposes Cauchy, as for every other law: the returned stress is
    ``tau / J1``. At this state J is 4 % away from 1, so a Cauchy-rate integration (or a
    missing boundary conversion) would show here at that size.
    """
    De0 = np.array([0.09, -0.03, -0.02, 0.0, 0.0, 0.0])
    cauchy, Lt, J = _hypoo(De0)
    assert abs(J - 1.0) > 0.03, "J must be away from 1 or a measure error is invisible"

    L_ortho = np.asarray(sim.L_ortho(HYPOO_PROPS[:9], "EnuG"))
    np.testing.assert_allclose(Lt, L_ortho, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(J * cauchy, L_ortho @ De0, rtol=1e-12, atol=1e-9)

    d = 1e-6
    for col in range(6):
        step = np.zeros(6)
        step[col] = d
        cp, _, Jp = _hypoo(De0 + step)
        cm, _, Jm = _hypoo(De0 - step)
        fd = (Jp * cp - Jm * cm) / (2.0 * d)
        np.testing.assert_allclose(fd, Lt[:, col], rtol=1e-6,
                                   atol=1e-6 * np.abs(Lt).max(),
                                   err_msg=f"HYPOO: d(tau)/d(De) column {col}")


@pytest.mark.parametrize("name,props,nstatev", KERNELS, ids=[k[0] for k in KERNELS])
def test_log_F_box_matches_finite_difference(name, props, nstatev):
    """Corate 5 (log_F): De = A^F:D dt and the stress is carried by sym(DF X DF^-1), so the box
    is c^J : (A^F)^-1 (chain rule). Checked by perturbing F along directions whose corate-5
    increment is a unit Voigt vector; for these isotropic kernels it also equals corate 3's box.
    (It used to return the Lie tangent, 8-13 % off.)"""
    from scipy.linalg import logm
    F0 = _F(EPS0)
    sigma0, Lt5 = _umat(name, props, F0, nstatev, corate=5)
    _, Lt3 = _umat(name, props, F0, nstatev, corate=3)
    tau0 = sim.v2t_stress(np.linalg.det(F0) * np.asarray(sigma0).ravel())
    lnV = lambda F: 0.5 * np.real(logm(F @ F.T))
    sym = lambda X: 0.5 * (X + X.T)
    A_F_inv = np.linalg.inv(np.asarray(sim.A_F(F0)))
    d = 1e-6
    fd = np.zeros((6, 6))
    for c in range(6):
        e = np.zeros(6)
        e[c] = 1.0
        M = sim.v2t_strain(A_F_inv @ e)
        num, den = [], []
        for sgn in (1.0, -1.0):
            DF = expm(sgn * d * M)
            F1 = DF @ F0
            sigma1, _ = _umat(name, props, F1, nstatev, corate=5)
            tau1 = np.linalg.det(F1) * np.asarray(sigma1).ravel()
            tau_hat = np.asarray(sim.t2v_stress(sym(DF @ tau0 @ np.linalg.inv(DF)))).ravel()
            De = np.asarray(sim.t2v_strain(lnV(F1) - sym(DF @ lnV(F0) @ np.linalg.inv(DF)))).ravel()
            num.append(tau1 - tau_hat)
            den.append(De[c])
        fd[:, c] = (num[0] - num[1]) / (den[0] - den[1])
    scale = np.abs(fd).max()
    np.testing.assert_allclose(Lt5, fd, atol=1e-6 * scale)
    np.testing.assert_allclose(Lt5, Lt3, atol=1e-9 * scale)


@pytest.mark.parametrize("corate", [-1, 6])
def test_umat_rejects_an_unknown_corate(corate):
    """sim.umat validates the corate up front, before the parallel region."""
    with pytest.raises(ValueError, match="corate"):
        _umat("NEOHC", [1000., 10000.], _F(EPS0), corate=corate)


def _v2t(v):
    return np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])


_WORK_CASES = [
    ("ELISO", [70000., 0.3, 0.], 1),
    ("ELORT", [150000., 10000., 10000., 0.3, 0.3, 0.45, 5000., 5000., 3500., 0., 0., 0.], 1),
    ("EPICP", [200000., 0.3, 0., 300., 1000., 0.5], 8),
    ("EPJCK", [200000., 0.33, 0., 792., 510., 0.26, 0.014, 1.0, 1.03, 293., 1793.], 9),
    ("HYPOO", [150000., 10000., 10000., 0.3, 0.3, 0.45, 5000., 5000., 3500., 0., 0., 0.], 1),
    ("SNTVE", [70000., 0.3, 0.], 1),
    ("NEOHC", [1000., 10000.], 1),
]


@pytest.mark.parametrize("name, props, nstatev", _WORK_CASES)
@pytest.mark.parametrize("corate", [0, 2, 3, 5])
def test_umat_work_is_the_stress_power_under_the_log_corates(name, props, nstatev, corate):
    """Non-coaxial state, arbitrary strain increment: under the log corates the Wm increment of
    sim.umat is the midpoint stress power 1/2 (tau_n + tau_n+1) : D dt with the lab start stress,
    as in the solver; under the others it stays the kernel's 1/2 (tau_n + tau_n+1) : De."""
    rng = np.random.default_rng(0)
    F0 = np.eye(3) + 0.2 * rng.standard_normal((3, 3))
    F0 = F0 if np.linalg.det(F0) > 0 else -F0
    F1 = (np.eye(3) + 0.05 * rng.standard_normal((3, 3))) @ F0
    R0, R1 = polar(F0)[0], polar(F1)[0]   # the polar increment; the caller has transported sig0 already
    De = 0.03 * rng.standard_normal(6)
    sig0 = 100. * rng.standard_normal(6)
    col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
    cube = lambda m: np.asarray(m, dtype=float).reshape(3, 3, 1).copy(order="F")
    statev = np.zeros((nstatev, 1), order="F")
    statev[0] = 290.
    sig1, _, Wm, _ = sim.umat(name, col(np.zeros(6)), col(De), cube(F0), cube(F1), col(sig0),
                              cube(R1 @ R0.T), col(props), statev, 0.5, 1.,
                              np.zeros((4, 1), order="F"), n_threads=1, corate=corate)
    tau0 = sig0 * np.linalg.det(F0)   # as passed: transported to the end frame by the caller
    tau1 = sig1[:, 0] * np.linalg.det(F1)
    if corate in (2, 3, 5):
        # stress power with the lab start stress DR^T tau0 DR, both stresses in one frame
        DR = R1 @ R0.T
        t0 = DR.T @ _v2t(tau0) @ DR
        tau0_lab = np.array([t0[0, 0], t0[1, 1], t0[2, 2], t0[0, 1], t0[0, 2], t0[1, 2]])
        Ldt = 2. * (F1 - F0) @ np.linalg.inv(F1 + F0)
        ref = 0.5 * np.dot(tau0_lab + tau1, np.asarray(sim.t2v_strain(0.5 * (Ldt + Ldt.T))).ravel())
    else:
        ref = 0.5 * np.dot(tau0 + tau1, De)
    np.testing.assert_allclose(Wm[0, 0], ref, rtol=1e-10, atol=1e-10 * abs(ref))


def test_umat_work_with_identity_F_is_the_kernel_work():
    """Identical F0 and F1 (small-strain use with placeholder deformation gradients) carry no
    increment: the log-corate work correction must not apply, Wm is the kernel's trapezoid."""
    col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    De = np.array([1e-3, -3e-4, -3e-4, 5e-4, 0., 0.])
    sig0 = np.array([50., 10., 10., 5., 0., 0.])
    statev = np.zeros((8, 1), order="F")
    statev[0] = 290.
    sig1, _, Wm, _ = sim.umat("EPICP", col(np.zeros(6)), col(De), eye, eye, col(sig0), eye,
                              col([200000., 0.3, 0., 300., 1000., 0.5]), statev, 0.5, 1.,
                              np.zeros((4, 1), order="F"), n_threads=1, corate=3)
    np.testing.assert_allclose(Wm[0, 0], 0.5 * np.dot(sig0 + sig1[:, 0], De), rtol=1e-12)


def test_umat_work_correction_can_be_turned_off():
    """work_correction=False (a caller whose F is not in the basis of its stresses): Wm is the
    kernel's own trapezoid 1/2 (tau_n + tau_n+1) : De even under a log corate with F0 != F1."""
    rng = np.random.default_rng(0)
    F0 = np.eye(3) + 0.2 * rng.standard_normal((3, 3))
    F0 = F0 if np.linalg.det(F0) > 0 else -F0
    F1 = (np.eye(3) + 0.05 * rng.standard_normal((3, 3))) @ F0
    De = 0.03 * rng.standard_normal(6)
    sig0 = 100. * rng.standard_normal(6)
    col = lambda a: np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))
    cube = lambda m: np.asarray(m, dtype=float).reshape(3, 3, 1).copy(order="F")
    statev = np.zeros((1, 1), order="F")
    statev[0] = 290.
    common = (col(np.zeros(6)), col(De), cube(F0), cube(F1), col(sig0), cube(np.eye(3)),
              col([70000., 0.3, 0.]), statev)
    sig1, _, Wm_off, _ = sim.umat("ELISO", *common, 0.5, 1., np.zeros((4, 1), order="F"),
                                  n_threads=1, corate=3, work_correction=False)
    tau0 = sig0 * np.linalg.det(F0)
    tau1 = sig1[:, 0] * np.linalg.det(F1)
    np.testing.assert_allclose(Wm_off[0, 0], 0.5 * np.dot(tau0 + tau1, De), rtol=1e-12)
    _, _, Wm_on, _ = sim.umat("ELISO", *common, 0.5, 1., np.zeros((4, 1), order="F"),
                              n_threads=1, corate=3)
    assert abs(Wm_on[0, 0] - Wm_off[0, 0]) > 1e-6 * abs(Wm_off[0, 0])   # default: corrected
