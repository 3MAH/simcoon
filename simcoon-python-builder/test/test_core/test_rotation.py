"""Tests for the Rotation class and scipy interoperability."""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation as ScipyRotation

import simcoon as sim


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def angle():
    return np.pi / 4


@pytest.fixture
def sim_rot(angle):
    """A simcoon.Rotation (scipy subclass)."""
    return sim.Rotation.from_axis_angle(angle, 3)


@pytest.fixture
def scipy_rot(angle):
    """A plain scipy Rotation representing the same rotation."""
    return ScipyRotation.from_rotvec([0, 0, angle])


@pytest.fixture
def L():
    """An isotropic stiffness matrix."""
    return sim.L_iso([70000, 0.3], "Enu")


@pytest.fixture
def M():
    """An isotropic compliance matrix."""
    return sim.M_iso([70000, 0.3], "Enu")


@pytest.fixture
def sigma():
    return np.array([100.0, 50.0, 0.0, 25.0, 0.0, 0.0])


@pytest.fixture
def epsilon():
    return np.array([0.01, -0.005, -0.005, 0.002, 0.001, 0.0])


# ===================================================================
# 1. Type identity
# ===================================================================

class TestTypeIdentity:
    """simcoon.Rotation IS a scipy Rotation."""

    def test_isinstance_scipy(self, sim_rot):
        assert isinstance(sim_rot, ScipyRotation)

    def test_isinstance_simcoon(self, sim_rot):
        assert isinstance(sim_rot, sim.Rotation)

    def test_from_euler_returns_simcoon(self):
        r = sim.Rotation.from_euler("ZXZ", [0.5, 0.3, 0.7])
        assert isinstance(r, sim.Rotation)

    def test_from_quat_returns_simcoon(self):
        r = sim.Rotation.from_quat([0, 0, 0.38268343, 0.92387953])
        assert isinstance(r, sim.Rotation)

    def test_from_matrix_returns_simcoon(self):
        R = np.eye(3)
        r = sim.Rotation.from_matrix(R)
        assert isinstance(r, sim.Rotation)

    def test_from_rotvec_returns_simcoon(self):
        r = sim.Rotation.from_rotvec([0, 0, 0.5])
        assert isinstance(r, sim.Rotation)

    def test_identity_returns_simcoon(self):
        r = sim.Rotation.identity()
        assert isinstance(r, sim.Rotation)

    def test_random_returns_simcoon(self):
        r = sim.Rotation.random()
        assert isinstance(r, sim.Rotation)

    def test_composition_preserves_type(self, sim_rot):
        r2 = sim.Rotation.from_axis_angle(np.pi / 6, 1)
        assert isinstance(sim_rot * r2, sim.Rotation)

    def test_inv_preserves_type(self, sim_rot):
        assert isinstance(sim_rot.inv(), sim.Rotation)


# ===================================================================
# 2. from_scipy conversion
# ===================================================================

class TestFromScipy:
    """Upgrade a plain scipy Rotation to simcoon.Rotation."""

    def test_from_scipy_type(self, scipy_rot):
        r = sim.Rotation.from_scipy(scipy_rot)
        assert isinstance(r, sim.Rotation)

    def test_from_scipy_preserves_quaternion(self, scipy_rot):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(r.as_quat(), scipy_rot.as_quat(), atol=1e-15)

    def test_from_scipy_preserves_matrix(self, scipy_rot):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(r.as_matrix(), scipy_rot.as_matrix(), atol=1e-14)

    def test_from_scipy_mechanics_work(self, scipy_rot, L):
        r = sim.Rotation.from_scipy(scipy_rot)
        L_rot = r.apply_stiffness(L)
        assert L_rot.shape == (6, 6)


# ===================================================================
# 3. Mechanics methods
# ===================================================================

class TestMechanics:
    """Mechanics methods on simcoon.Rotation delegate to C++ correctly."""

    def test_apply_stress_shape(self, sim_rot, sigma):
        result = sim_rot.apply_stress(sigma)
        assert result.shape == (6,)

    def test_apply_strain_shape(self, sim_rot, epsilon):
        result = sim_rot.apply_strain(epsilon)
        assert result.shape == (6,)

    def test_apply_stiffness_shape(self, sim_rot, L):
        result = sim_rot.apply_stiffness(L)
        assert result.shape == (6, 6)

    def test_apply_compliance_shape(self, sim_rot, M):
        result = sim_rot.apply_compliance(M)
        assert result.shape == (6, 6)

    def test_apply_tensor_shape(self, sim_rot):
        T = np.eye(3)
        result = sim_rot.apply_tensor(T)
        assert result.shape == (3, 3)

    def test_apply_strain_concentration_shape(self, sim_rot):
        A = np.eye(6)
        result = sim_rot.apply_strain_concentration(A)
        assert result.shape == (6, 6)

    def test_apply_stress_concentration_shape(self, sim_rot):
        B = np.eye(6)
        result = sim_rot.apply_stress_concentration(B)
        assert result.shape == (6, 6)

    def test_voigt_stress_rotation_shape(self, sim_rot):
        QS = sim_rot.as_voigt_stress_rotation()
        assert QS.shape == (6, 6)

    def test_voigt_strain_rotation_shape(self, sim_rot):
        QE = sim_rot.as_voigt_strain_rotation()
        assert QE.shape == (6, 6)

    def test_identity_stiffness_unchanged(self, L):
        r = sim.Rotation.identity()
        np.testing.assert_allclose(r.apply_stiffness(L), L, atol=1e-10)

    def test_roundtrip_stress(self, sim_rot, sigma):
        sigma_rot = sim_rot.apply_stress(sigma)
        sigma_back = sim_rot.apply_stress(sigma_rot, active=False)
        np.testing.assert_allclose(sigma_back, sigma, atol=1e-10)

    def test_roundtrip_stiffness(self, sim_rot, L):
        L_rot = sim_rot.apply_stiffness(L)
        L_back = sim_rot.apply_stiffness(L_rot, active=False)
        np.testing.assert_allclose(L_back, L, atol=1e-10)

    def test_stiffness_via_voigt(self, sim_rot, L):
        """apply_stiffness(L) == QS @ L @ QS.T"""
        QS = sim_rot.as_voigt_stress_rotation()
        L_rot_direct = sim_rot.apply_stiffness(L)
        L_rot_manual = QS @ L @ QS.T
        np.testing.assert_allclose(L_rot_direct, L_rot_manual, atol=1e-10)


# ===================================================================
# 4. Mechanics from scipy-created rotation (via from_scipy)
# ===================================================================

class TestMechanicsFromScipy:
    """Verify that a scipy Rotation upgraded via from_scipy produces
    identical results to a natively-created simcoon Rotation."""

    def test_stiffness_matches(self, sim_rot, scipy_rot, L):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.apply_stiffness(L), sim_rot.apply_stiffness(L), atol=1e-10
        )

    def test_compliance_matches(self, sim_rot, scipy_rot, M):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.apply_compliance(M), sim_rot.apply_compliance(M), atol=1e-10
        )

    def test_stress_matches(self, sim_rot, scipy_rot, sigma):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.apply_stress(sigma), sim_rot.apply_stress(sigma), atol=1e-10
        )

    def test_strain_matches(self, sim_rot, scipy_rot, epsilon):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.apply_strain(epsilon), sim_rot.apply_strain(epsilon), atol=1e-10
        )

    def test_tensor_matches(self, sim_rot, scipy_rot):
        T = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.apply_tensor(T), sim_rot.apply_tensor(T), atol=1e-10
        )

    def test_voigt_stress_rotation_matches(self, sim_rot, scipy_rot):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.as_voigt_stress_rotation(),
            sim_rot.as_voigt_stress_rotation(),
            atol=1e-10,
        )

    def test_voigt_strain_rotation_matches(self, sim_rot, scipy_rot):
        r = sim.Rotation.from_scipy(scipy_rot)
        np.testing.assert_allclose(
            r.as_voigt_strain_rotation(),
            sim_rot.as_voigt_strain_rotation(),
            atol=1e-10,
        )


# ===================================================================
# 5. Utility methods
# ===================================================================

class TestUtilities:

    def test_equals_same(self, sim_rot):
        r2 = sim.Rotation.from_axis_angle(np.pi / 4, 3)
        assert sim_rot.equals(r2)

    def test_equals_different(self, sim_rot):
        r2 = sim.Rotation.from_axis_angle(np.pi / 3, 3)
        assert not sim_rot.equals(r2)

    def test_equals_with_scipy(self, sim_rot, scipy_rot):
        assert sim_rot.equals(scipy_rot)

    def test_is_identity_true(self):
        r = sim.Rotation.identity()
        assert r.is_identity()

    def test_is_identity_false(self, sim_rot):
        assert not sim_rot.is_identity()

    def test_slerp_to_endpoints(self, sim_rot):
        r_id = sim.Rotation.identity()
        r0 = r_id.slerp_to(sim_rot, 0.0)
        r1 = r_id.slerp_to(sim_rot, 1.0)
        np.testing.assert_allclose(r0.as_quat(), r_id.as_quat(), atol=1e-12)
        np.testing.assert_allclose(
            np.abs(r1.as_quat() @ sim_rot.as_quat()), 1.0, atol=1e-12
        )

    def test_slerp_to_with_scipy_target(self, scipy_rot):
        r_id = sim.Rotation.identity()
        r_mid = r_id.slerp_to(scipy_rot, 0.5)
        assert isinstance(r_mid, sim.Rotation)

    def test_from_axis_angle_degrees(self):
        r_rad = sim.Rotation.from_axis_angle(np.pi / 4, 3)
        r_deg = sim.Rotation.from_axis_angle(45, 3, degrees=True)
        np.testing.assert_allclose(r_rad.as_quat(), r_deg.as_quat(), atol=1e-12)

    def test_from_axis_angle_invalid_axis(self):
        with pytest.raises(ValueError, match="axis must be 1, 2, or 3"):
            sim.Rotation.from_axis_angle(0.5, 4)
