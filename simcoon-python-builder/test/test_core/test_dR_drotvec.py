"""Tests for ``dR_drotvec``: derivatives of R(omega) w.r.t. rotation-vector components.

Layout convention (matches simcoon's project-wide cube convention):
    single : shape (3, 3, 3), ``result[:, :, k]`` is dR/d(omega_k)
    batch  : shape (3, 3, 3, N), ``result[:, :, k, n]`` is dR_n/d(omega_k)
"""

import numpy as np
import pytest

import simcoon as sim


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def rotvec_general():
    """Generic non-axis-aligned rotvec (moderate angle ~0.37 rad)."""
    return np.array([0.1, -0.2, 0.3])


@pytest.fixture
def rotvecs_batch():
    """A batch of diverse rotvecs, including identity."""
    return np.array([
        [0.1, -0.2, 0.3],
        [0.0, 0.0, 0.5],
        [0.7, 0.0, 0.0],
        [0.0, 0.0, 0.0],                  # identity
        [np.pi / 2 - 0.05, 0.0, 0.0],     # near pi/2
    ])


# Reference skew matrices [e_k]x for k=0,1,2
E_SKEW = np.array([
    [[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]],   # [e_0]x
    [[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],   # [e_1]x
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],   # [e_2]x
])


# ---------------------------------------------------------------------------
# Shape + API symmetry
# ---------------------------------------------------------------------------

class TestShape:
    def test_single_shape(self, rotvec_general):
        assert sim.dR_drotvec(rotvec_general).shape == (3, 3, 3)

    def test_member_shape(self, rotvec_general):
        r = sim.Rotation.from_rotvec(rotvec_general)
        assert r.dR_drotvec().shape == (3, 3, 3)

    def test_batch_shape(self, rotvecs_batch):
        r = sim.Rotation.from_rotvec(rotvecs_batch)
        assert r.dR_drotvec().shape == (3, 3, 3, len(rotvecs_batch))

    def test_free_eq_member(self, rotvec_general):
        # The member path goes through Rotation._ensure_cls -> from_quat, which
        # re-normalises the quaternion and introduces ~1 ULP precision loss
        # vs. the input rotvec. Compare with a tight tolerance instead of
        # bit-exact equality.
        r = sim.Rotation.from_rotvec(rotvec_general)
        np.testing.assert_allclose(
            sim.dR_drotvec(rotvec_general), r.dR_drotvec(), atol=1e-12
        )


# ---------------------------------------------------------------------------
# Small-angle / identity branch
# ---------------------------------------------------------------------------

class TestSmallAngle:
    def test_zero_rotvec_equals_basis_skews(self):
        """At omega=0: dR/d(omega_k) == [e_k]x for each k."""
        result = sim.dR_drotvec(np.zeros(3))
        for k in range(3):
            np.testing.assert_allclose(result[:, :, k], E_SKEW[k], atol=0.0)

    def test_identity_rotation_matches_zero_rotvec(self):
        result = sim.Rotation.identity().dR_drotvec()
        expected = sim.dR_drotvec(np.zeros(3))
        np.testing.assert_allclose(result, expected, atol=0.0)

    def test_below_iota_threshold_uses_small_angle_branch(self):
        """Very small rotvec should give near-basis-skew slices."""
        tiny = np.array([1e-14, 0.0, 0.0])
        result = sim.dR_drotvec(tiny)
        for k in range(3):
            np.testing.assert_allclose(result[:, :, k], E_SKEW[k], atol=1e-12)


# ---------------------------------------------------------------------------
# Finite-difference agreement
# ---------------------------------------------------------------------------

def _R(omega):
    return sim.Rotation.from_rotvec(omega).as_matrix()


def _fd_slice(omega, k, eps=1e-6):
    pert = omega.copy()
    pert[k] += eps
    return (_R(pert) - _R(omega - np.eye(3)[k] * eps)) / (2.0 * eps)


class TestFiniteDifference:
    @pytest.mark.parametrize("rotvec", [
        np.array([0.1, -0.2, 0.3]),
        np.array([0.0, 0.0, 0.9]),
        np.array([1.2, -0.4, 0.3]),           # |omega| ~ 1.3 rad
        np.array([np.pi / 2 - 0.01, 0.0, 0.0]),
    ])
    def test_fd_matches_analytic(self, rotvec):
        analytic = sim.dR_drotvec(rotvec)
        for k in range(3):
            fd = _fd_slice(rotvec, k)
            np.testing.assert_allclose(analytic[:, :, k], fd, atol=1e-7, rtol=1e-7)


# ---------------------------------------------------------------------------
# Batch consistency
# ---------------------------------------------------------------------------

class TestBatch:
    def test_batch_slice_matches_single_call(self, rotvecs_batch):
        r = sim.Rotation.from_rotvec(rotvecs_batch)
        batch = r.dR_drotvec()
        for n, rv in enumerate(rotvecs_batch):
            single = sim.dR_drotvec(rv)
            np.testing.assert_allclose(batch[:, :, :, n], single, atol=0.0)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestValidation:
    def test_rejects_wrong_size(self):
        with pytest.raises(Exception):
            sim.dR_drotvec(np.zeros(2))

    def test_rejects_2d_input(self):
        with pytest.raises(Exception):
            sim.dR_drotvec(np.zeros((3, 1)))
