"""
Tests for unified Tensor2/Tensor4 batch operations.

Cross-validates every batch operation against single-tensor loops using the
existing C++ single-tensor methods as reference.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation as ScipyRotation

import simcoon as smc
from simcoon.tensor import (
    Tensor2,
    Tensor4,
    double_contract,
    _mat_to_voigt,
    _voigt_to_mat,
)


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

N = 10
TOL = 1e-12
ROTATION_TOL = 1e-10


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def random_sym_matrices(rng):
    raw = rng.standard_normal((N, 3, 3))
    return 0.5 * (raw + raw.transpose(0, 2, 1))


@pytest.fixture
def random_stress_voigt(rng):
    return rng.standard_normal((N, 6))


@pytest.fixture
def random_strain_voigt(rng):
    return rng.standard_normal((N, 6))


@pytest.fixture
def iso_stiffness():
    return smc.L_iso([70000.0, 0.3], "Enu")


@pytest.fixture
def iso_compliance(iso_stiffness):
    return np.linalg.inv(iso_stiffness)


@pytest.fixture
def random_rotations(rng):
    return ScipyRotation.random(N, random_state=rng)


@pytest.fixture
def single_rotation(rng):
    return ScipyRotation.random(random_state=rng)


# ==================================================================
# Voigt conversion helpers
# ==================================================================

class TestVoigtConversion:
    def test_stress_mat_voigt_roundtrip(self, random_sym_matrices):
        v = _mat_to_voigt(random_sym_matrices, "stress")
        m2 = _voigt_to_mat(v, "stress")
        assert_allclose(m2, random_sym_matrices, atol=TOL)

    def test_strain_mat_voigt_roundtrip(self, random_sym_matrices):
        v = _mat_to_voigt(random_sym_matrices, "strain")
        m2 = _voigt_to_mat(v, "strain")
        assert_allclose(m2, random_sym_matrices, atol=TOL)

    def test_stress_voigt_matches_single(self, random_sym_matrices):
        batch_v = _mat_to_voigt(random_sym_matrices, "stress")
        for idx in range(N):
            single = Tensor2.stress(random_sym_matrices[idx].copy())
            assert_allclose(batch_v[idx], single.voigt, atol=TOL)

    def test_strain_voigt_matches_single(self, random_sym_matrices):
        batch_v = _mat_to_voigt(random_sym_matrices, "strain")
        for idx in range(N):
            single = Tensor2.strain(random_sym_matrices[idx].copy())
            assert_allclose(batch_v[idx], single.voigt, atol=TOL)


# ==================================================================
# Tensor2 batch: Sequence-like protocol
# ==================================================================

class TestTensor2BatchSequence:
    def test_batch_has_len(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        assert len(batch) == N

    def test_single_flag(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        assert not batch.single
        single = Tensor2.stress(random_stress_voigt[0])
        assert single.single

    def test_iter(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        tensors = list(batch)
        assert len(tensors) == N
        for i, t in enumerate(tensors):
            assert isinstance(t, Tensor2)
            assert t.single
            assert_allclose(t.voigt, random_stress_voigt[i], atol=TOL)

    def test_reversed(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        rev = list(reversed(batch))
        assert len(rev) == N
        assert_allclose(rev[0].voigt, random_stress_voigt[-1], atol=TOL)
        assert_allclose(rev[-1].voigt, random_stress_voigt[0], atol=TOL)

    def test_contains(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        t0 = Tensor2.from_voigt(random_stress_voigt[0].copy(), "stress")
        assert t0 in batch

    def test_construct_from_list(self):
        tensors = [Tensor2.stress(np.array([float(i)] * 6)) for i in range(5)]
        batch = Tensor2(tensors)
        assert len(batch) == 5
        for i in range(5):
            assert_allclose(batch[i].voigt, tensors[i].voigt, atol=TOL)

    def test_construct_from_list_mixed_type_raises(self):
        t1 = Tensor2.stress(np.zeros(6))
        t2 = Tensor2.strain(np.zeros(6))
        with pytest.raises(ValueError, match="Mixed type"):
            Tensor2([t1, t2])

    def test_negative_index(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        assert_allclose(batch[-1].voigt, batch[N - 1].voigt, atol=TOL)

    def test_index_out_of_range(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        with pytest.raises(IndexError):
            batch[N]

    def test_count(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        t0 = Tensor2.from_voigt(random_stress_voigt[0].copy(), "stress")
        assert batch.count(t0) >= 1


# ==================================================================
# Tensor2 batch: numpy array protocol
# ==================================================================

class TestTensor2BatchNumpy:
    def test_asarray(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        arr = np.asarray(batch)
        assert arr.shape == (N, 6)
        assert_allclose(arr, random_stress_voigt, atol=TOL)

    def test_asarray_dtype(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        arr = np.asarray(batch, dtype=np.float32)
        assert arr.dtype == np.float32
        assert arr.shape == (N, 6)

    def test_np_add(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt * 2)
        c = np.add(a, b)
        assert isinstance(c, Tensor2)
        assert_allclose(c.voigt, random_stress_voigt * 3, atol=TOL)

    def test_np_multiply_scalar(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        result = np.multiply(batch, 3.0)
        assert isinstance(result, Tensor2)
        assert_allclose(result.voigt, random_stress_voigt * 3, atol=TOL)

    def test_np_subtract(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt)
        c = np.subtract(a, b)
        assert isinstance(c, Tensor2)
        assert_allclose(c.voigt, 0.0, atol=TOL)

    def test_np_negative(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        neg = np.negative(batch)
        assert isinstance(neg, Tensor2)
        assert_allclose(neg.voigt, -random_stress_voigt, atol=TOL)

    def test_unsupported_ufunc(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        with pytest.raises(TypeError):
            np.sin(batch)


# ==================================================================
# Tensor2 batch: concatenation
# ==================================================================

class TestTensor2BatchConcatenate:
    def test_concatenate_two(self, random_stress_voigt):
        b1 = Tensor2.stress(random_stress_voigt[:5])
        b2 = Tensor2.stress(random_stress_voigt[5:])
        combined = Tensor2.concatenate([b1, b2])
        assert len(combined) == N
        assert_allclose(combined.voigt, random_stress_voigt, atol=TOL)

    def test_concatenate_three(self, random_stress_voigt):
        parts = [Tensor2.stress(random_stress_voigt[i:i+4])
                 for i in range(0, 8, 4)]
        parts.append(Tensor2.stress(random_stress_voigt[8:]))
        combined = Tensor2.concatenate(parts)
        assert len(combined) == N

    def test_concatenate_mixed_type_raises(self):
        b1 = Tensor2.stress(np.zeros((3, 6)))
        b2 = Tensor2.strain(np.zeros((3, 6)))
        with pytest.raises(ValueError, match="Mixed type"):
            Tensor2.concatenate([b1, b2])


# ==================================================================
# Tensor2 batch: construction from arrays
# ==================================================================

class TestTensor2BatchConstruction:
    def test_from_voigt(self, random_stress_voigt):
        batch = Tensor2.from_voigt(random_stress_voigt, "stress")
        assert_allclose(batch.voigt, random_stress_voigt, atol=TOL)
        assert batch.vtype == "stress"
        assert len(batch) == N

    def test_from_mat(self, random_sym_matrices):
        batch = Tensor2.from_mat(random_sym_matrices, "stress")
        assert len(batch) == N
        assert_allclose(batch.mat, random_sym_matrices, atol=TOL)

    def test_stress_factory(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        assert batch.vtype == "stress"
        assert_allclose(batch.voigt, random_stress_voigt, atol=TOL)

    def test_strain_factory(self, random_strain_voigt):
        batch = Tensor2.strain(random_strain_voigt)
        assert batch.vtype == "strain"

    def test_stress_from_matrices(self, random_sym_matrices):
        batch = Tensor2.stress(random_sym_matrices)
        assert batch.vtype == "stress"
        assert len(batch) == N

    def test_from_tensor_broadcast(self):
        t = Tensor2.stress(np.array([100., 200., 300., 10., 20., 30.]))
        batch = Tensor2.from_tensor(t, 5)
        assert len(batch) == 5
        for i in range(5):
            assert_allclose(batch[i].voigt, t.voigt, atol=TOL)

    def test_from_columns_roundtrip(self, random_stress_voigt):
        arr_6N = random_stress_voigt.T
        batch = Tensor2.from_columns(arr_6N, "stress")
        assert_allclose(batch.voigt_T, arr_6N, atol=TOL)

    def test_from_list_alias(self):
        tensors = [Tensor2.stress(np.array([float(i)] * 6)) for i in range(3)]
        batch = Tensor2.from_list(tensors)
        assert len(batch) == 3


# ==================================================================
# Tensor2 batch: properties and indexing
# ==================================================================

class TestTensor2BatchProperties:
    def test_mat_property(self, random_sym_matrices):
        batch = Tensor2.from_mat(random_sym_matrices, "stress")
        assert batch.mat.shape == (N, 3, 3)
        assert_allclose(batch.mat, random_sym_matrices, atol=TOL)

    def test_voigt_T(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        assert batch.voigt_T.shape == (6, N)

    def test_getitem_int(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        t = batch[3]
        assert isinstance(t, Tensor2)
        assert_allclose(t.voigt, random_stress_voigt[3], atol=TOL)

    def test_getitem_slice(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        sub = batch[2:5]
        assert isinstance(sub, Tensor2)
        assert not sub.single
        assert len(sub) == 3
        assert_allclose(sub.voigt, random_stress_voigt[2:5], atol=TOL)

    def test_getitem_mask(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        mask = np.array([True, False] * 5)
        sub = batch[mask]
        assert isinstance(sub, Tensor2)
        assert len(sub) == 5


# ==================================================================
# Tensor2 batch: scalar invariants (C++ batch loops)
# ==================================================================

class TestTensor2BatchInvariants:
    def test_trace(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        traces = batch.trace()
        for i in range(N):
            assert_allclose(traces[i], batch[i].trace(), atol=TOL)

    def test_mises_stress(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        vm = batch.mises()
        for i in range(N):
            assert_allclose(vm[i], batch[i].mises(), atol=TOL)

    def test_mises_strain(self, random_strain_voigt):
        batch = Tensor2.strain(random_strain_voigt)
        vm = batch.mises()
        for i in range(N):
            assert_allclose(vm[i], batch[i].mises(), atol=TOL)

    def test_dev(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        dev = batch.dev()
        assert_allclose(dev.trace(), 0.0, atol=1e-10)

    def test_norm(self, random_sym_matrices):
        batch = Tensor2.from_mat(random_sym_matrices, "stress")
        norms = batch.norm()
        for i in range(N):
            expected = np.sqrt(np.sum(random_sym_matrices[i] ** 2))
            assert_allclose(norms[i], expected, atol=TOL)


# ==================================================================
# Tensor2 batch: arithmetic
# ==================================================================

class TestTensor2BatchArithmetic:
    def test_add_batch(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt * 2)
        c = a + b
        assert_allclose(c.voigt, random_stress_voigt * 3, atol=TOL)

    def test_add_single(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        single = Tensor2.stress(np.ones(6))
        result = batch + single
        assert_allclose(result.voigt, random_stress_voigt + 1.0, atol=TOL)

    def test_radd_single(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        single = Tensor2.stress(np.ones(6))
        result = single + batch
        assert_allclose(result.voigt, random_stress_voigt + 1.0, atol=TOL)

    def test_sub_batch(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt)
        c = a - b
        assert_allclose(c.voigt, 0.0, atol=TOL)

    def test_neg(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        neg = -batch
        assert_allclose(neg.voigt, -random_stress_voigt, atol=TOL)

    def test_mul_scalar(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        result = batch * 3.0
        assert_allclose(result.voigt, random_stress_voigt * 3.0, atol=TOL)

    def test_rmul_scalar(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        result = 3.0 * batch
        assert_allclose(result.voigt, random_stress_voigt * 3.0, atol=TOL)

    def test_mul_per_point(self, rng, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        factors = rng.standard_normal(N)
        result = batch * factors
        for i in range(N):
            assert_allclose(result.voigt[i], random_stress_voigt[i] * factors[i],
                            atol=TOL)

    def test_div_scalar(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        result = batch / 2.0
        assert_allclose(result.voigt, random_stress_voigt / 2.0, atol=TOL)

    def test_double_contraction_mod(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt * 0.5)
        result = a % b
        for i in range(N):
            ma = a[i].mat
            mb = b[i].mat
            expected = np.sum(ma * mb)
            assert_allclose(result[i], expected, atol=TOL)


# ==================================================================
# Tensor2 batch: rotation (C++ batch loop)
# ==================================================================

class TestTensor2BatchRotation:
    def test_rotate_stress_matches_single(self, random_stress_voigt, random_rotations):
        batch = Tensor2.stress(random_stress_voigt)
        rotated = batch.rotate(random_rotations)
        for i in range(N):
            R_single = random_rotations[i]
            single_rot = batch[i].rotate(smc.Rotation.from_scipy(R_single))
            assert_allclose(rotated[i].voigt, single_rot.voigt, atol=ROTATION_TOL)

    def test_rotate_strain_matches_single(self, random_strain_voigt, random_rotations):
        batch = Tensor2.strain(random_strain_voigt)
        rotated = batch.rotate(random_rotations)
        for i in range(N):
            R_single = random_rotations[i]
            single_rot = batch[i].rotate(smc.Rotation.from_scipy(R_single))
            assert_allclose(rotated[i].voigt, single_rot.voigt, atol=ROTATION_TOL)

    def test_rotate_single_broadcast(self, random_stress_voigt, single_rotation):
        batch = Tensor2.stress(random_stress_voigt)
        rotated = batch.rotate(single_rotation)
        R = smc.Rotation.from_scipy(single_rotation)
        for i in range(N):
            single_rot = batch[i].rotate(R)
            assert_allclose(rotated[i].voigt, single_rot.voigt, atol=ROTATION_TOL)

    def test_rotate_roundtrip(self, random_stress_voigt, random_rotations):
        batch = Tensor2.stress(random_stress_voigt)
        rotated = batch.rotate(random_rotations)
        restored = rotated.rotate(random_rotations.inv())
        assert_allclose(restored.voigt, batch.voigt, atol=ROTATION_TOL)


# ==================================================================
# Tensor2 batch: push-forward / pull-back (C++ batch loop)
# ==================================================================

class TestTensor2BatchPushPull:
    def test_push_forward_pull_back_stress(self, rng, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        F = np.eye(3) + 0.1 * rng.standard_normal((N, 3, 3))
        pushed = batch.push_forward(F)
        pulled = pushed.pull_back(F)
        assert_allclose(pulled.voigt, batch.voigt, atol=1e-8)

    def test_push_forward_pull_back_strain(self, rng, random_strain_voigt):
        batch = Tensor2.strain(random_strain_voigt)
        F = np.eye(3) + 0.1 * rng.standard_normal((N, 3, 3))
        pushed = batch.push_forward(F)
        pulled = pushed.pull_back(F)
        assert_allclose(pulled.voigt, batch.voigt, atol=1e-8)

    def test_push_forward_single_F_broadcast(self, rng, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        F = np.eye(3) + 0.05 * rng.standard_normal((3, 3))
        pushed = batch.push_forward(F)
        for i in range(N):
            single_pushed = batch[i].push_forward(F)
            assert_allclose(pushed[i].voigt, single_pushed.voigt, atol=1e-10)


# ==================================================================
# Tensor4 batch: Sequence-like protocol
# ==================================================================

class TestTensor4BatchSequence:
    def test_batch_has_len(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        assert len(batch) == 5

    def test_iter(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        tensors = list(batch)
        assert len(tensors) == 5
        for t in tensors:
            assert isinstance(t, Tensor4)
            assert t.single
            assert_allclose(t.mat, iso_stiffness, atol=TOL)

    def test_construct_from_list(self, iso_stiffness):
        L = Tensor4.stiffness(iso_stiffness)
        batch = Tensor4([L, L, L])
        assert len(batch) == 3
        for t in batch:
            assert_allclose(t.mat, iso_stiffness, atol=TOL)

    def test_negative_index(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        assert_allclose(batch[-1].mat, batch[4].mat, atol=TOL)

    def test_index_out_of_range(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        with pytest.raises(IndexError):
            batch[5]


# ==================================================================
# Tensor4 batch: numpy array protocol
# ==================================================================

class TestTensor4BatchNumpy:
    def test_asarray(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        arr = np.asarray(batch)
        assert arr.shape == (5, 6, 6)
        assert_allclose(arr[0], iso_stiffness, atol=TOL)

    def test_np_add(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        result = np.add(L, L)
        assert isinstance(result, Tensor4)
        assert_allclose(result.voigt, L.voigt * 2, atol=TOL)

    def test_np_multiply_scalar(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), 5)
        result = np.multiply(L, 2.0)
        assert isinstance(result, Tensor4)
        assert_allclose(result.voigt, L.voigt * 2, atol=TOL)


# ==================================================================
# Tensor4 batch: concatenation
# ==================================================================

class TestTensor4BatchConcatenate:
    def test_concatenate(self, iso_stiffness):
        L = Tensor4.stiffness(iso_stiffness)
        b1 = Tensor4.from_tensor(L, 3)
        b2 = Tensor4.from_tensor(L, 4)
        combined = Tensor4.concatenate([b1, b2])
        assert len(combined) == 7

    def test_concatenate_mixed_type_raises(self, iso_stiffness, iso_compliance):
        b1 = Tensor4.stiffness(np.broadcast_to(iso_stiffness[np.newaxis], (3, 6, 6)).copy())
        b2 = Tensor4.compliance(np.broadcast_to(iso_compliance[np.newaxis], (3, 6, 6)).copy())
        with pytest.raises(ValueError, match="Mixed type"):
            Tensor4.concatenate([b1, b2])


# ==================================================================
# Tensor4 batch: construction from arrays
# ==================================================================

class TestTensor4BatchConstruction:
    def test_stiffness(self, iso_stiffness):
        data = np.broadcast_to(iso_stiffness[np.newaxis], (N, 6, 6)).copy()
        batch = Tensor4.stiffness(data)
        assert len(batch) == N
        assert batch.type == "stiffness"

    def test_compliance(self, iso_compliance):
        data = np.broadcast_to(iso_compliance[np.newaxis], (N, 6, 6)).copy()
        batch = Tensor4.compliance(data)
        assert batch.type == "compliance"

    def test_from_tensor(self, iso_stiffness):
        single = Tensor4.stiffness(iso_stiffness)
        batch = Tensor4.from_tensor(single, N)
        assert len(batch) == N
        for i in range(N):
            assert_allclose(batch[i].mat, iso_stiffness, atol=TOL)

    def test_from_list_alias(self, iso_stiffness):
        tensors = [Tensor4.stiffness(iso_stiffness * (1 + 0.01 * i))
                    for i in range(5)]
        batch = Tensor4.from_list(tensors)
        assert len(batch) == 5
        for i in range(5):
            assert_allclose(batch[i].mat, tensors[i].mat, atol=TOL)

    def test_getitem_int(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        t = batch[0]
        assert isinstance(t, Tensor4)
        assert_allclose(t.mat, iso_stiffness, atol=TOL)

    def test_getitem_slice(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        sub = batch[2:5]
        assert isinstance(sub, Tensor4)
        assert not sub.single
        assert len(sub) == 3


# ==================================================================
# Tensor4 batch: contraction (C++ batch loop)
# ==================================================================

class TestTensor4BatchContraction:
    def test_contraction_matches_single(self, iso_stiffness, random_strain_voigt):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        eps = Tensor2.strain(random_strain_voigt)
        sigma = L @ eps
        assert isinstance(sigma, Tensor2)
        assert not sigma.single
        assert sigma.vtype == "stress"
        for i in range(N):
            single_sigma = L[i] @ eps[i]
            assert_allclose(sigma[i].voigt, single_sigma.voigt, atol=TOL)

    def test_contraction_broadcast_single_tensor2(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        eps = Tensor2.strain(np.array([0.01, -0.003, -0.003, 0.005, 0.002, 0.001]))
        sigma = L @ eps
        assert len(sigma) == N
        for i in range(1, N):
            assert_allclose(sigma[i].voigt, sigma[0].voigt, atol=TOL)

    def test_compliance_contraction_gives_strain(self, iso_stiffness):
        M_mat = np.linalg.inv(iso_stiffness)
        M = Tensor4.from_tensor(Tensor4.compliance(M_mat), N)
        sigma_v = np.zeros((N, 6))
        sigma_v[:, 0] = 100.0
        sigma = Tensor2.stress(sigma_v)
        eps = M @ sigma
        assert eps.vtype == "strain"


# ==================================================================
# Tensor4 batch: rotation (C++ batch loop)
# ==================================================================

class TestTensor4BatchRotation:
    def test_rotate_stiffness_matches_single(self, iso_stiffness, random_rotations):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        rotated = L.rotate(random_rotations)
        for i in range(N):
            R_single = smc.Rotation.from_scipy(random_rotations[i])
            single_rot = L[i].rotate(R_single)
            assert_allclose(rotated[i].mat, single_rot.mat, atol=ROTATION_TOL)

    def test_rotate_compliance_matches_single(self, iso_compliance, random_rotations):
        M = Tensor4.from_tensor(Tensor4.compliance(iso_compliance), N)
        rotated = M.rotate(random_rotations)
        for i in range(N):
            R_single = smc.Rotation.from_scipy(random_rotations[i])
            single_rot = M[i].rotate(R_single)
            assert_allclose(rotated[i].mat, single_rot.mat, atol=ROTATION_TOL)

    def test_rotate_roundtrip(self, iso_stiffness, random_rotations):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        rotated = L.rotate(random_rotations)
        restored = rotated.rotate(random_rotations.inv())
        for i in range(N):
            assert_allclose(restored[i].mat, iso_stiffness, atol=1e-8)

    def test_material_frame_equivalence(self, rng, random_rotations):
        """L @ eps in material frame == (rotate L) @ (rotate eps) in global frame."""
        L_mat = smc.L_iso([70000.0, 0.3], "Enu")
        L_batch = Tensor4.from_tensor(Tensor4.stiffness(L_mat), N)
        eps_voigt = rng.standard_normal((N, 6)) * 0.01
        eps_batch = Tensor2.strain(eps_voigt)

        sigma_mat = L_batch @ eps_batch
        L_global = L_batch.rotate(random_rotations)
        eps_global = eps_batch.rotate(random_rotations)
        sigma_global = L_global @ eps_global
        sigma_mat_to_global = sigma_mat.rotate(random_rotations)

        assert_allclose(sigma_global.voigt, sigma_mat_to_global.voigt, atol=1e-6)


# ==================================================================
# Tensor4 batch: inverse (C++ batch loop)
# ==================================================================

class TestTensor4BatchInverse:
    def test_inverse_type_flip(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        M = L.inverse()
        assert M.type == "compliance"
        L2 = M.inverse()
        assert L2.type == "stiffness"
        assert_allclose(L2.voigt, L.voigt, atol=1e-8)

    def test_inverse_matches_single(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        M = L.inverse()
        M_single = Tensor4.stiffness(iso_stiffness).inverse()
        for i in range(N):
            assert_allclose(M[i].mat, M_single.mat, atol=1e-8)


# ==================================================================
# Tensor4 batch: arithmetic
# ==================================================================

class TestTensor4BatchArithmetic:
    def test_add(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        result = L + L
        assert_allclose(result.voigt, L.voigt * 2, atol=TOL)

    def test_add_single(self, iso_stiffness):
        L_batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        L_single = Tensor4.stiffness(iso_stiffness)
        result = L_batch + L_single
        assert_allclose(result.voigt, L_batch.voigt * 2, atol=TOL)

    def test_sub(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        result = L - L
        assert_allclose(result.voigt, 0.0, atol=TOL)

    def test_neg(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        result = -L
        assert_allclose(result.voigt, -L.voigt, atol=TOL)

    def test_mul_scalar(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        result = L * 2.0
        assert_allclose(result.voigt, L.voigt * 2, atol=TOL)

    def test_mul_per_point(self, rng, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        factors = rng.standard_normal(N)
        result = L * factors
        for i in range(N):
            assert_allclose(result[i].mat, iso_stiffness * factors[i], atol=TOL)

    def test_div_scalar(self, iso_stiffness):
        L = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        result = L / 2.0
        assert_allclose(result.voigt, L.voigt / 2, atol=TOL)


# ==================================================================
# Free function: double_contract
# ==================================================================

class TestDoubleContract:
    def test_matches_mod_operator(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        b = Tensor2.stress(random_stress_voigt * 0.5)
        result1 = double_contract(a, b)
        result2 = a % b
        assert_allclose(result1, result2, atol=TOL)

    def test_self_contraction_positive(self, random_stress_voigt):
        a = Tensor2.stress(random_stress_voigt)
        result = double_contract(a, a)
        assert np.all(result >= 0)


# ==================================================================
# Interop
# ==================================================================

class TestInterop:
    def test_columns_roundtrip(self, rng):
        arr_6N = rng.standard_normal((6, N))
        batch = Tensor2.from_columns(arr_6N, "stress")
        assert_allclose(batch.voigt_T, arr_6N, atol=TOL)

    def test_single_to_batch_and_back(self):
        t = Tensor2.stress(np.array([100., -50., -50., 30., 20., 10.]))
        batch = Tensor2.from_tensor(t, 5)
        for i in range(5):
            assert_allclose(batch[i].voigt, t.voigt, atol=TOL)
            assert batch[i].vtype == t.vtype


# ==================================================================
# Repr
# ==================================================================

class TestRepr:
    def test_tensor2_batch_repr(self, random_stress_voigt):
        batch = Tensor2.stress(random_stress_voigt)
        r = repr(batch)
        assert "Tensor2" in r
        assert str(N) in r

    def test_tensor4_batch_repr(self, iso_stiffness):
        batch = Tensor4.from_tensor(Tensor4.stiffness(iso_stiffness), N)
        r = repr(batch)
        assert "Tensor4" in r
        assert str(N) in r


# ==================================================================
# Transparency tests (unified API)
# ==================================================================

class TestUnifiedTransparency:
    def test_single_flag(self):
        t = Tensor2.stress(np.zeros(6))
        assert t.single is True

    def test_batch_flag(self):
        t = Tensor2.stress(np.zeros((5, 6)))
        assert t.single is False

    def test_single_len_raises(self):
        t = Tensor2.stress(np.zeros(6))
        with pytest.raises(TypeError):
            len(t)

    def test_single_no_getitem(self):
        t = Tensor2.stress(np.zeros(6))
        with pytest.raises(TypeError):
            t[0]

    def test_batch_indexing_returns_single(self):
        t = Tensor2.stress(np.zeros((5, 6)))
        assert t[0].single is True

    def test_batch_iteration(self):
        t = Tensor2.stress(np.zeros((3, 6)))
        singles = list(t)
        assert len(singles) == 3
        assert all(s.single for s in singles)

    def test_single_matmul_returns_single(self):
        L = Tensor4.stiffness(smc.L_iso([70000, 0.3], "Enu"))
        eps = Tensor2.strain(np.array([0.01, 0, 0, 0, 0, 0]))
        sigma = L @ eps
        assert sigma.single is True

    def test_batch_matmul_returns_batch(self):
        L = Tensor4.stiffness(smc.L_iso([70000, 0.3], "Enu"))
        eps = Tensor2.strain(np.zeros((5, 6)))
        sigma = L @ eps
        assert sigma.single is False
        assert len(sigma) == 5

