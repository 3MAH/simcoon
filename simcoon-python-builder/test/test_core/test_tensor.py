"""Tests for the Tensor2 and Tensor4 classes."""

import numpy as np
import pytest

import simcoon as sim


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sigma_mat():
    return np.array([[100, 30, 20], [30, 200, 40], [20, 40, 300]], dtype=float)


@pytest.fixture
def eps_mat():
    return np.array([[0.01, 0.005, 0.003], [0.005, 0.02, 0.004], [0.003, 0.004, 0.03]], dtype=float)


@pytest.fixture
def sigma_voigt():
    return np.array([100, 200, 300, 30, 20, 40], dtype=float)


@pytest.fixture
def eps_voigt():
    return np.array([0.01, 0.02, 0.03, 0.01, 0.006, 0.008], dtype=float)


@pytest.fixture
def L():
    return sim.L_iso([70000, 0.3], "Enu")


@pytest.fixture
def M():
    return sim.M_iso([70000, 0.3], "Enu")


@pytest.fixture
def F():
    return np.array([[1.1, 0.1, 0.05], [0.02, 0.95, 0.03], [0.01, 0.04, 1.05]])


@pytest.fixture
def rot():
    return sim.Rotation.from_euler("ZXZ", [0.3, 0.5, 0.7])


# ===================================================================
# Tensor2 tests
# ===================================================================

class TestTensor2Construction:

    def test_stress_from_mat(self, sigma_mat):
        t = sim.Tensor2.stress(sigma_mat)
        assert t.vtype == sim.VoigtType.stress
        np.testing.assert_allclose(t.mat, sigma_mat, atol=1e-14)

    def test_strain_from_mat(self, eps_mat):
        t = sim.Tensor2.strain(eps_mat)
        assert t.vtype == sim.VoigtType.strain
        np.testing.assert_allclose(t.mat, eps_mat, atol=1e-14)

    def test_stress_from_voigt(self, sigma_voigt):
        t = sim.Tensor2.stress(sigma_voigt)
        assert t.vtype == sim.VoigtType.stress
        np.testing.assert_allclose(t.voigt, sigma_voigt, atol=1e-12)

    def test_strain_from_voigt(self, eps_voigt):
        t = sim.Tensor2.strain(eps_voigt)
        assert t.vtype == sim.VoigtType.strain
        np.testing.assert_allclose(t.voigt, eps_voigt, atol=1e-12)

    def test_zeros(self):
        t = sim.Tensor2.zeros(sim.VoigtType.strain)
        np.testing.assert_allclose(t.mat, np.zeros((3, 3)), atol=1e-15)

    def test_identity(self):
        t = sim.Tensor2.identity()
        np.testing.assert_allclose(t.mat, np.eye(3), atol=1e-15)

    def test_bad_shape_raises(self):
        with pytest.raises(ValueError):
            sim.Tensor2.stress(np.zeros(5))


class TestTensor2VoigtConventions:

    def test_stress_voigt_matches_t2v_stress(self, sigma_mat):
        t = sim.Tensor2.stress(sigma_mat)
        ref = sim.t2v_stress(sigma_mat).ravel()
        np.testing.assert_allclose(t.voigt, ref, atol=1e-12)

    def test_strain_voigt_matches_t2v_strain(self, eps_mat):
        t = sim.Tensor2.strain(eps_mat)
        ref = sim.t2v_strain(eps_mat).ravel()
        np.testing.assert_allclose(t.voigt, ref, atol=1e-12)

    def test_stress_voigt_roundtrip(self, sigma_voigt):
        t = sim.Tensor2.stress(sigma_voigt)
        np.testing.assert_allclose(t.voigt, sigma_voigt, atol=1e-12)

    def test_strain_voigt_roundtrip(self, eps_voigt):
        t = sim.Tensor2.strain(eps_voigt)
        np.testing.assert_allclose(t.voigt, eps_voigt, atol=1e-12)


class TestTensor2Arithmetic:

    def test_addition(self, sigma_mat):
        t1 = sim.Tensor2.stress(sigma_mat)
        t2 = sim.Tensor2.stress(sigma_mat)
        t3 = t1 + t2
        np.testing.assert_allclose(t3.mat, 2.0 * sigma_mat, atol=1e-12)

    def test_scalar_multiply(self, sigma_mat):
        t = sim.Tensor2.stress(sigma_mat)
        t2 = t * 3.0
        t3 = 3.0 * t
        np.testing.assert_allclose(t2.mat, 3.0 * sigma_mat, atol=1e-12)
        np.testing.assert_allclose(t3.mat, 3.0 * sigma_mat, atol=1e-12)

    def test_equality(self, sigma_mat):
        t1 = sim.Tensor2.stress(sigma_mat)
        t2 = sim.Tensor2.stress(sigma_mat)
        t3 = sim.Tensor2.strain(sigma_mat)
        assert t1 == t2
        assert not (t1 == t3)


class TestTensor2FreeFunctions:

    def test_trace(self, sigma_mat):
        t = sim.Tensor2.stress(sigma_mat)
        assert abs(t.trace() - 600.0) < 1e-10

    def test_mises_uniaxial(self):
        Y = 250.0
        sigma = np.array([Y, 0, 0, 0, 0, 0], dtype=float)
        t = sim.Tensor2.stress(sigma)
        assert abs(t.mises() - Y) < 1e-10

    def test_mises_hydrostatic_zero(self):
        sigma = np.array([100, 100, 100, 0, 0, 0], dtype=float)
        t = sim.Tensor2.stress(sigma)
        assert t.mises() < 1e-10


class TestTensor2Rotation:

    def test_stress_rotation_roundtrip(self, sigma_mat, rot):
        t = sim.Tensor2.stress(sigma_mat)
        t_rot = t.rotate(rot, active=True)
        inv_rot = rot.inv()
        t_back = t_rot.rotate(inv_rot, active=True)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-9)

    def test_strain_rotation_roundtrip(self, eps_mat, rot):
        t = sim.Tensor2.strain(eps_mat)
        t_rot = t.rotate(rot, active=True)
        inv_rot = rot.inv()
        t_back = t_rot.rotate(inv_rot, active=True)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-9)


class TestTensor2PushPull:

    def test_stress_push_pull_roundtrip(self, sigma_mat, F):
        t = sim.Tensor2.stress(sigma_mat)
        t_push = t.push_forward(F)
        t_back = t_push.pull_back(F)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-9)

    def test_strain_push_pull_roundtrip(self, eps_mat, F):
        t = sim.Tensor2.strain(eps_mat)
        t_push = t.push_forward(F)
        t_back = t_push.pull_back(F)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-9)


# ===================================================================
# Tensor4 tests
# ===================================================================

class TestTensor4Construction:

    def test_stiffness(self, L):
        t = sim.Tensor4.stiffness(L)
        assert t.type == sim.Tensor4Type.stiffness
        np.testing.assert_allclose(t.mat, L, atol=1e-12)

    def test_compliance(self, M):
        t = sim.Tensor4.compliance(M)
        assert t.type == sim.Tensor4Type.compliance
        np.testing.assert_allclose(t.mat, M, atol=1e-12)

    def test_identity_matches_Ireal(self):
        t = sim.Tensor4.identity()
        ref = sim.Ireal()
        np.testing.assert_allclose(t.mat, ref, atol=1e-12)

    def test_volumetric_matches_Ivol(self):
        t = sim.Tensor4.volumetric()
        np.testing.assert_allclose(t.mat, sim.Ivol(), atol=1e-12)

    def test_deviatoric_matches_Idev(self):
        t = sim.Tensor4.deviatoric()
        np.testing.assert_allclose(t.mat, sim.Idev(), atol=1e-12)


class TestTensor4Contraction:

    def test_stiffness_contract_gives_stress(self, L):
        stiff = sim.Tensor4.stiffness(L)
        eps_v = np.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0])
        eps = sim.Tensor2.strain(eps_v)
        sigma = stiff.contract(eps)
        assert sigma.vtype == sim.VoigtType.stress
        np.testing.assert_allclose(sigma.voigt, L @ eps_v, atol=1e-10)

    def test_compliance_contract_gives_strain(self, M):
        comp = sim.Tensor4.compliance(M)
        sig_v = np.array([100, 0, 0, 0, 0, 0], dtype=float)
        sig = sim.Tensor2.stress(sig_v)
        eps = comp.contract(sig)
        assert eps.vtype == sim.VoigtType.strain
        np.testing.assert_allclose(eps.voigt, M @ sig_v, atol=1e-10)

    def test_matmul_operator(self, L):
        stiff = sim.Tensor4.stiffness(L)
        eps_v = np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001])
        eps = sim.Tensor2.strain(eps_v)
        sigma = stiff @ eps
        np.testing.assert_allclose(sigma.voigt, L @ eps_v, atol=1e-10)


class TestTensor4Arithmetic:

    def test_isotropic_from_projectors(self, L):
        E = 70000.0
        nu = 0.3
        K = E / (3.0 * (1.0 - 2.0 * nu))
        mu = E / (2.0 * (1.0 + nu))
        L_proj = 3.0 * K * sim.Tensor4.volumetric() + 2.0 * mu * sim.Tensor4.deviatoric()
        np.testing.assert_allclose(L_proj.mat, L, atol=1e-8)

    def test_scalar_multiply(self, L):
        t = sim.Tensor4.stiffness(L)
        t2 = t * 2.0
        t3 = 2.0 * t
        np.testing.assert_allclose(t2.mat, 2.0 * L, atol=1e-12)
        np.testing.assert_allclose(t3.mat, 2.0 * L, atol=1e-12)


class TestTensor4Rotation:

    def test_isotropic_stiffness_invariant(self, L, rot):
        """Isotropic stiffness is invariant under any rotation."""
        t = sim.Tensor4.stiffness(L)
        t_rot = t.rotate(rot, active=True)
        np.testing.assert_allclose(t.mat, t_rot.mat, atol=1e-8)

    def test_isotropic_compliance_invariant(self, M, rot):
        t = sim.Tensor4.compliance(M)
        t_rot = t.rotate(rot, active=True)
        np.testing.assert_allclose(t.mat, t_rot.mat, atol=1e-8)


class TestTensor4PushPull:

    def test_push_pull_roundtrip(self, L, F):
        t = sim.Tensor4.stiffness(L)
        t_push = t.push_forward(F)
        t_back = t_push.pull_back(F)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)


# ===================================================================
# Integration tests
# ===================================================================

class TestIntegration:

    def test_full_cycle_stress_rotate_contract(self, L, rot):
        """L:eps in material frame == (rotate L) : (rotate eps) in global frame."""
        stiff = sim.Tensor4.stiffness(L)
        eps_v = np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001])
        eps = sim.Tensor2.strain(eps_v)

        # Material frame
        sigma = stiff @ eps

        # Global frame
        stiff_rot = stiff.rotate(rot, active=True)
        eps_rot = eps.rotate(rot, active=True)
        sigma_rot_direct = stiff_rot @ eps_rot

        # Rotate material-frame stress to global frame
        sigma_rot_indirect = sigma.rotate(rot, active=True)

        np.testing.assert_allclose(
            sigma_rot_direct.mat, sigma_rot_indirect.mat, atol=1e-8
        )

    def test_repr(self):
        t2 = sim.Tensor2.stress(np.eye(3))
        assert "stress" in repr(t2)

        t4 = sim.Tensor4.stiffness(np.eye(6))
        assert "stiffness" in repr(t4)
