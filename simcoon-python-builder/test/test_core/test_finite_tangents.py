"""Finite-difference the box tangent of every finite-strain kernel.

Why this file exists: **a wrong tangent does not move converged solver values.** Newton
finds the same root with a poor Jacobian, so `regression_baseline.py` -- 100 combos of
umat x control type x corate -- is structurally blind to this entire class of error. It
caught neither of the two real tangent defects found in this area:

* ``saint_venant`` fed its stress to a helper that rebuilds ``tau = det(F)*sigma``, i.e.
  expects **Cauchy**, after the kernel had been made Kirchhoff-native -- squaring the J;
* ``HYPOO`` hands over ``dsigma/dD`` where the consumer reads ``d(tau_hat)/dDe``.

Both are invisible at J = 1, so every state here carries J well away from 1 and the
docstrings state the relative size a J error would have.
"""

import numpy as np
import pytest
from scipy.linalg import expm

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


def test_hypoelastic_tangent_carries_the_J_and_the_stress_term():
    """HYPOO integrates a corotational CAUCHY rate, so its L is dsigma/dD, not the box.

    With ``tau = J sigma`` and ``dJ/dt = J tr(D)``, the consistent box is
    ``J (L + sigma (x) I)``. Handing over the bare ``L`` -- which this kernel did -- is wrong
    by both factors: at this state 4.1 % from the J and 7.6 % from the stress term, 9.9 %
    together. Both vanish at J = 1 and zero stress, which is how it survived.

    HYPOO is a RATE kernel: it responds to ``Detot``, not to the total strain, so the
    difference is taken along the increment and F1 is kept consistent as exp(De).
    """
    def _hypoo(De):
        De = np.asarray(De, dtype=float)
        F1 = expm(sim.v2t_strain(De)).reshape(3, 3, 1).copy(order="F")
        eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
        stress, sv, wm, Lt = sim.umat(
            "HYPOO", np.zeros((6, 1), order="F"), np.asfortranarray(De.reshape(6, 1)),
            eye, F1, np.zeros((6, 1), order="F"), eye,
            np.asfortranarray(np.asarray(HYPOO_PROPS, dtype=float).reshape(-1, 1)),
            np.zeros((1, 1), order="F"), 0.0, 1.0,
            np.zeros((4, 1), order="F"), n_threads=1, corate=3)
        return stress[:, 0], Lt[:, :, 0], float(np.linalg.det(F1[:, :, 0]))

    De0 = np.array([0.09, -0.03, -0.02, 0.0, 0.0, 0.0])
    _, Lt, J = _hypoo(De0)
    assert abs(J - 1.0) > 0.03, "J must be away from 1 or the J half is invisible"

    def tau_of(De):
        stress, _, Jl = _hypoo(De)
        return Jl * np.asarray(stress).ravel()

    d = 1e-6
    for col in range(3):
        step = np.zeros(6)
        step[col] = d
        fd = (tau_of(De0 + step) - tau_of(De0 - step)) / (2.0 * d)
        np.testing.assert_allclose(
            fd[:3], Lt[:3, col], rtol=2e-6,
            atol=2e-6 * max(1.0, np.abs(Lt[:3, col]).max()),
            err_msg=f"HYPOO: d(tau)/d(De) column {col}")

    # and the bare dsigma/dD really is a different matrix, so the check above has teeth
    L_ortho = np.asarray(sim.L_ortho(HYPOO_PROPS[:9], "EnuG"))
    assert np.linalg.norm(Lt - L_ortho) > 0.05 * np.linalg.norm(Lt)


def test_log_F_box_is_the_spatial_tangent_itself():
    """corate 5 is the convected/Oldroyd-Lie box: an identity, not a spectral map."""
    name, props, nstatev = KERNELS[1]          # NEOHC
    _, lie = _umat(name, props, _F(EPS0), nstatev, corate=5)
    _, log_box = _umat(name, props, _F(EPS0), nstatev, corate=3)
    assert np.all(np.isfinite(lie))
    assert np.abs(lie - log_box).max() > 1e-6 * np.abs(log_box).max()
