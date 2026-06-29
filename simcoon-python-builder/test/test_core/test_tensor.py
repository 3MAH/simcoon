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
        assert t.vtype == "stress"
        np.testing.assert_allclose(t.mat, sigma_mat, atol=1e-14)

    def test_strain_from_mat(self, eps_mat):
        t = sim.Tensor2.strain(eps_mat)
        assert t.vtype == "strain"
        np.testing.assert_allclose(t.mat, eps_mat, atol=1e-14)

    def test_stress_from_voigt(self, sigma_voigt):
        t = sim.Tensor2.stress(sigma_voigt)
        assert t.vtype == "stress"
        np.testing.assert_allclose(t.voigt, sigma_voigt, atol=1e-12)

    def test_strain_from_voigt(self, eps_voigt):
        t = sim.Tensor2.strain(eps_voigt)
        assert t.vtype == "strain"
        np.testing.assert_allclose(t.voigt, eps_voigt, atol=1e-12)

    def test_zeros(self):
        t = sim.Tensor2.zeros("strain")
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
        assert t.type == "stiffness"
        np.testing.assert_allclose(t.mat, L, atol=1e-12)

    def test_compliance(self, M):
        t = sim.Tensor4.compliance(M)
        assert t.type == "compliance"
        np.testing.assert_allclose(t.mat, M, atol=1e-12)

    def test_strain_concentration(self):
        t = sim.Tensor4.strain_concentration(np.eye(6))
        assert t.type == "strain_concentration"

    def test_stress_concentration(self):
        t = sim.Tensor4.stress_concentration(np.eye(6))
        assert t.type == "stress_concentration"

    # Projector factories are type-driven: the type tag picks the Voigt convention, so a single
    # identity/deviatoric covers every type (the legacy identity2/deviatoric2 are gone). Compared
    # against explicit reference matrices rather than the deprecated free Ireal()/Idev() helpers.
    _IREAL = np.diag([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])   # stiffness identity (Ireal)
    _IREAL2 = np.diag([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])  # compliance identity (Ireal2)

    @staticmethod
    def _ivol():
        v = np.zeros((6, 6)); v[:3, :3] = 1.0 / 3.0
        return v

    def test_identity_stiffness_is_Ireal(self):
        np.testing.assert_allclose(sim.Tensor4.identity("stiffness").mat, self._IREAL, atol=1e-12)

    def test_identity_compliance_is_Ireal2(self):
        # The strain-convention identity (former identity2/Ireal2), now reached via the type tag.
        np.testing.assert_allclose(sim.Tensor4.identity("compliance").mat, self._IREAL2, atol=1e-12)

    def test_identity_concentration_is_eye(self):
        for t in ("strain_concentration", "stress_concentration"):
            np.testing.assert_allclose(sim.Tensor4.identity(t).mat, np.eye(6), atol=1e-12)

    def test_volumetric_is_Ivol(self):
        np.testing.assert_allclose(sim.Tensor4.volumetric("stiffness").mat, self._ivol(), atol=1e-12)

    def test_deviatoric_stiffness_is_Idev(self):
        np.testing.assert_allclose(sim.Tensor4.deviatoric("stiffness").mat,
                                   self._IREAL - self._ivol(), atol=1e-12)

    def test_deviatoric_compliance_is_Idev2(self):
        np.testing.assert_allclose(sim.Tensor4.deviatoric("compliance").mat,
                                   self._IREAL2 - self._ivol(), atol=1e-12)


class TestTensor4Contraction:

    def test_stiffness_contract_gives_stress(self, L):
        stiff = sim.Tensor4.stiffness(L)
        eps_v = np.array([0.01, 0.0, 0.0, 0.0, 0.0, 0.0])
        eps = sim.Tensor2.strain(eps_v)
        sigma = stiff.contract(eps)
        assert sigma.vtype == "stress"
        np.testing.assert_allclose(sigma.voigt, L @ eps_v, atol=1e-10)

    def test_compliance_contract_gives_strain(self, M):
        comp = sim.Tensor4.compliance(M)
        sig_v = np.array([100, 0, 0, 0, 0, 0], dtype=float)
        sig = sim.Tensor2.stress(sig_v)
        eps = comp.contract(sig)
        assert eps.vtype == "strain"
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

    def test_strain_concentration_roundtrip(self, L, rot):
        """Strain concentration: rotate then rotate back gives original."""
        t = sim.Tensor4.strain_concentration(L)
        t_rot = t.rotate(rot, active=True)
        inv_rot = rot.inv()
        t_back = t_rot.rotate(inv_rot, active=True)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)

    def test_stress_concentration_roundtrip(self, L, rot):
        """Stress concentration: rotate then rotate back gives original."""
        t = sim.Tensor4.stress_concentration(L)
        t_rot = t.rotate(rot, active=True)
        inv_rot = rot.inv()
        t_back = t_rot.rotate(inv_rot, active=True)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)

    def test_anisotropic_compliance_roundtrip(self, rot):
        """Anisotropic compliance: rotate and back gives original."""
        L_cubic = sim.L_cubic([200000, 0.3, 80000], "EnuG")
        M = np.linalg.inv(L_cubic)
        comp = sim.Tensor4.compliance(M)
        comp_rot = comp.rotate(rot, active=True)
        inv_rot = rot.inv()
        comp_back = comp_rot.rotate(inv_rot, active=True)
        np.testing.assert_allclose(comp.mat, comp_back.mat, atol=1e-8)


class TestTensor4PushPull:

    def test_push_pull_roundtrip(self, L, F):
        t = sim.Tensor4.stiffness(L)
        t_push = t.push_forward(F)
        t_back = t_push.pull_back(F)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)

    def test_compliance_push_pull_roundtrip(self, M, F):
        t = sim.Tensor4.compliance(M)
        t_push = t.push_forward(F)
        t_back = t_push.pull_back(F)
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)

    def test_compliance_push_forward_consistency(self):
        """Compliance push-forward with simple shear: verify specific value.

        A pure (12,12) shear compliance pushed forward by F^{-T} under x-y simple
        shear maps into the 22-normal component, not the 11. F^{-T} row 1 = [1,0,0]
        leaves [0,0] at 0 (C'_1111 = C_1111 = 0); F^{-T} row 2 = [-1,1,0] gives
        C'_2222 = 4 -> Voigt [1,1] = 4.0 (normal-normal factor 1).
        """
        M = np.zeros((6, 6))
        M[3, 3] = 4.0
        comp = sim.Tensor4.compliance(M)
        F_shear = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        comp_push = comp.push_forward(F_shear)
        assert abs(comp_push.mat[1, 1] - 4.0) < 1e-10
        assert abs(comp_push.mat[0, 0]) < 1e-10


class TestTensor4Inverse:

    def test_stiffness_to_compliance(self, L, M):
        stiff = sim.Tensor4.stiffness(L)
        comp = stiff.inverse()
        assert comp.type == "compliance"
        np.testing.assert_allclose(comp.mat, M, atol=1e-8)

    def test_compliance_to_stiffness(self, L, M):
        comp = sim.Tensor4.compliance(M)
        stiff = comp.inverse()
        assert stiff.type == "stiffness"
        np.testing.assert_allclose(stiff.mat, L, atol=1e-8)

    def test_inverse_roundtrip(self, L):
        t = sim.Tensor4.stiffness(L)
        t_back = t.inverse().inverse()
        assert t_back.type == "stiffness"
        np.testing.assert_allclose(t.mat, t_back.mat, atol=1e-8)


class TestTensor4Dyadic:

    def test_dyadic(self):
        a = sim.Tensor2.stress(np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float))
        b = sim.Tensor2.stress(np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]], dtype=float))
        d = sim.dyadic(a, b)
        assert d.type == "stiffness"
        assert abs(d.mat[0, 1] - 1.0) < 1e-14
        assert abs(d.mat[0, 0]) < 1e-14

    def test_auto_dyadic(self):
        a = sim.Tensor2.stress(np.eye(3))
        ad = sim.auto_dyadic(a)
        assert abs(ad.mat[0, 0] - 1.0) < 1e-14
        assert abs(ad.mat[0, 1] - 1.0) < 1e-14
        assert abs(ad.mat[3, 3]) < 1e-14


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


# ---------------------------------------------------------------------------
# Concentration tensors (strain/stress) — the engineering Voigt convention
# ---------------------------------------------------------------------------

# Voigt index <-> tensor index-pair map (no factor-2 on shears here; the
# concentration factor is applied explicitly per the documented convention).
_VOIGT_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]


def _rand_minor_sym_t4(seed):
    """Random 4th-order tensor with minor symmetries A_ijkl = A_jikl = A_ijlk."""
    rng = np.random.default_rng(seed)
    B = rng.standard_normal((3, 3, 3, 3))
    return 0.25 * (B + B.transpose(1, 0, 2, 3)
                   + B.transpose(0, 1, 3, 2) + B.transpose(1, 0, 3, 2))


def _t4_to_conc_voigt(A, scale_rows):
    """Full-index A_ijkl -> concentration Voigt 6x6. scale_rows=True doubles the shear rows
    (strain_concentration); False doubles the shear cols (stress_concentration)."""
    V = np.empty((6, 6))
    for I, (i, j) in enumerate(_VOIGT_PAIRS):
        for J, (k, l) in enumerate(_VOIGT_PAIRS):
            V[I, J] = A[i, j, k, l] * (2.0 if (I if scale_rows else J) >= 3 else 1.0)
    return V


# Symmetric 3x3 -> Voigt via the shipped conversions (engineering shears for strain, plain
# for stress); ravel() because t2v_* return (6,1). Reused rather than re-implemented inline.
def _sym_to_strain_voigt(e):
    return np.asarray(sim.t2v_strain(e)).ravel()


def _sym_to_stress_voigt(s):
    return np.asarray(sim.t2v_stress(s)).ravel()


class TestTensor4Concentration:
    """contract/inverse for concentration tensors are the sensitive operations —
    these lock the engineering Voigt convention from first principles."""

    def test_strain_concentration_identity_roundtrip(self):
        """identity('strain_concentration') must return engineering strain unchanged."""
        eps_v = np.array([0.02, -0.01, -0.005, 0.012, 0.008, 0.004])
        A = sim.Tensor4.identity("strain_concentration")
        np.testing.assert_allclose(np.diag(A.mat), np.ones(6), atol=1e-14)
        out = A.contract(sim.Tensor2.strain(eps_v))
        assert out.vtype == "strain"
        np.testing.assert_allclose(out.voigt, eps_v, atol=1e-12)

    def test_stress_concentration_identity_roundtrip(self):
        """identity('stress_concentration') must return stress unchanged."""
        sig_v = np.array([100.0, 50, -30, 20, 10, 5])
        B = sim.Tensor4.identity("stress_concentration")
        np.testing.assert_allclose(np.diag(B.mat), np.ones(6), atol=1e-14)
        out = B.contract(sim.Tensor2.stress(sig_v))
        assert out.vtype == "stress"
        np.testing.assert_allclose(out.voigt, sig_v, atol=1e-12)

    def test_eye_strain_concentration_roundtrip(self):
        """A bare eye(6) is the identity in the engineering convention."""
        eps_v = np.array([0.01, -0.02, 0.005, 0.03, -0.01, 0.02])
        A = sim.Tensor4.strain_concentration(np.eye(6))
        np.testing.assert_allclose(A.contract(sim.Tensor2.strain(eps_v)).voigt,
                                   eps_v, atol=1e-12)

    def test_strain_concentration_full_index_ground_truth(self):
        """Typed Voigt contract == full-index A_ijkl : eps_kl, from first principles."""
        A = _rand_minor_sym_t4(seed=1)
        rng = np.random.default_rng(2)
        e = rng.standard_normal((3, 3))
        e = 0.5 * (e + e.T)                          # symmetric global strain
        eps_local = np.einsum("ijkl,kl->ij", A, e)   # full-index ground truth
        V = _t4_to_conc_voigt(A, scale_rows=True)
        typed = sim.Tensor4.strain_concentration(V).contract(
            sim.Tensor2.strain(_sym_to_strain_voigt(e)))
        np.testing.assert_allclose(typed.voigt, _sym_to_strain_voigt(eps_local),
                                   atol=1e-12)

    def test_stress_concentration_full_index_ground_truth(self):
        """Typed Voigt contract == full-index B_ijkl : sig_kl."""
        Bt = _rand_minor_sym_t4(seed=3)
        rng = np.random.default_rng(4)
        s = rng.standard_normal((3, 3))
        s = 0.5 * (s + s.T)
        sig_local = np.einsum("ijkl,kl->ij", Bt, s)
        V = _t4_to_conc_voigt(Bt, scale_rows=False)
        typed = sim.Tensor4.stress_concentration(V).contract(
            sim.Tensor2.stress(_sym_to_stress_voigt(s)))
        np.testing.assert_allclose(typed.voigt, _sym_to_stress_voigt(sig_local),
                                   atol=1e-12)

    def test_strain_concentration_A_R_producer(self, F):
        """A_R (De = A^R:D) wraps and contracts identically to raw matmul."""
        A_R = sim.A_R(F)
        D = np.array([[0.03, 0.01, 0.0], [0.01, -0.02, 0.0], [0.0, 0.0, -0.01]])
        D_eng = np.asarray(sim.t2v_strain(D)).ravel()
        typed = sim.Tensor4.strain_concentration(A_R).contract(
            sim.Tensor2.strain(D_eng))
        np.testing.assert_allclose(typed.voigt, A_R @ D_eng, atol=1e-12)

    def test_strain_concentration_inverse(self, F):
        """inverse keeps the type, equals numpy inverse, and round-trips contract."""
        A_R = sim.A_R(F)
        Ac = sim.Tensor4.strain_concentration(A_R)
        Ainv = Ac.inverse()
        assert Ainv.type == "strain_concentration"
        np.testing.assert_allclose(Ainv.mat, np.linalg.inv(A_R), atol=1e-10)
        x = sim.Tensor2.strain(np.array([0.02, -0.01, -0.005, 0.012, 0.008, 0.004]))
        back = Ainv.contract(Ac.contract(x))
        np.testing.assert_allclose(back.voigt, x.voigt, atol=1e-12)

    def test_stress_concentration_inverse(self):
        """stress concentration inverse keeps type and round-trips."""
        V = _t4_to_conc_voigt(_rand_minor_sym_t4(seed=5), scale_rows=False) + 3.0 * np.eye(6)
        Bc = sim.Tensor4.stress_concentration(V)
        Binv = Bc.inverse()
        assert Binv.type == "stress_concentration"
        np.testing.assert_allclose(Binv.mat, np.linalg.inv(V), atol=1e-10)
        x = sim.Tensor2.stress(np.array([100.0, 50, -30, 20, 10, 5]))
        back = Binv.contract(Bc.contract(x))
        np.testing.assert_allclose(back.voigt, x.voigt, atol=1e-10)

    def test_concentration_push_pull_throws(self, F):
        """push/pull are undefined for concentration tensors (mixed indices)."""
        A = sim.Tensor4.strain_concentration(np.eye(6))
        B = sim.Tensor4.stress_concentration(np.eye(6))
        for t in (A, B):
            with pytest.raises(RuntimeError):
                t.push_forward(F)
            with pytest.raises(RuntimeError):
                t.pull_back(F)

    def test_rotate_accepts_scipy_rotation(self, F):
        """Single Tensor4.rotate accepts a scipy Rotation (not only simcoon.Rotation)."""
        from scipy.spatial.transform import Rotation as ScipyRotation
        A = sim.Tensor4.strain_concentration(sim.A_R(F))
        R = ScipyRotation.from_euler("z", 25, degrees=True)
        rotated = A.rotate(R)
        # objectivity: rotating A_R(F) equals A_R(R F R^T)
        Rm = R.as_matrix()
        np.testing.assert_allclose(rotated.mat, sim.A_R(Rm @ F @ Rm.T), atol=1e-7)

    def test_rotate_rejects_bad_type(self):
        """An unrecognized rotation argument raises TypeError (not an opaque C++ error)."""
        A = sim.Tensor4.strain_concentration(np.eye(6))
        with pytest.raises(TypeError):
            A.rotate(np.eye(3))


# ---------------------------------------------------------------------------
# Anisotropic guards on the Mandel-internal storage (the engineering API is
# unchanged: these exercise the general non-isotropic stiffness/compliance path
# that the isotropic L_iso fixtures cannot discriminate).
# ---------------------------------------------------------------------------

def _aniso_spd_stiffness(seed=7):
    """A generic anisotropic symmetric positive-definite 6x6 test stiffness (builder-free)."""
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((6, 6))
    return A @ A.T + 6.0 * np.eye(6)


class TestTensor4MandelInvariants:

    def test_anisotropic_mat_roundtrip(self):
        """Constructor (eng->Mandel) then .mat (Mandel->eng) is the identity for every type."""
        L = _aniso_spd_stiffness()
        for factory in (sim.Tensor4.stiffness, sim.Tensor4.compliance,
                        sim.Tensor4.strain_concentration, sim.Tensor4.stress_concentration):
            np.testing.assert_allclose(factory(L).mat, L, rtol=1e-11, atol=1e-11)

    def test_anisotropic_inverse_roundtrip(self):
        L = _aniso_spd_stiffness()
        K = sim.Tensor4.stiffness(L)
        np.testing.assert_allclose(K.inverse().inverse().mat, L, rtol=1e-9, atol=1e-9)

    def test_anisotropic_stiffness_compliance_contract(self):
        """M:(L:eps) == eps for an anisotropic stiffness (contract + inverse via Mandel)."""
        L = _aniso_spd_stiffness()
        K = sim.Tensor4.stiffness(L)
        Minv = K.inverse()
        assert Minv.type == "compliance"
        eps_v = np.array([0.01, -0.004, -0.003, 0.006, 0.002, 0.0015])
        eps_back = Minv @ (K @ sim.Tensor2.strain(eps_v))
        np.testing.assert_allclose(eps_back.voigt, eps_v, atol=1e-10)
