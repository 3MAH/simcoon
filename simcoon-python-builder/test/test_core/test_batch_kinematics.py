"""Batch (3,3,N) forms of RU_decomposition, VR_decomposition, isochoric_pstretch and stress_convert:
the same values as the single-matrix calls, point by point (they run in parallel past 100 points)."""

import numpy as np
import pytest

import simcoon as sim

N = 150  # past the parallel cutoff (100)


def _F(n=N, seed=3):
    rng = np.random.default_rng(seed)
    return np.asfortranarray(np.eye(3)[:, :, None] + 0.2 * rng.standard_normal((3, 3, n)))


def _slice(a, i):
    return np.ascontiguousarray(a[:, :, i])


@pytest.mark.parametrize("name", ["RU_decomposition", "VR_decomposition"])
def test_polar_batch_is_the_single_call(name):
    F = _F()
    A, B = getattr(sim, name)(F)
    assert A.shape == B.shape == (3, 3, N)
    for i in range(N):
        a, b = getattr(sim, name)(_slice(F, i))
        np.testing.assert_array_equal(A[:, :, i], a)
        np.testing.assert_array_equal(B[:, :, i], b)


@pytest.mark.parametrize("tensor", ["V", "b"])
def test_isochoric_pstretch_batch(tensor):
    F = _F()
    V, _ = sim.VR_decomposition(F)
    X = V if tensor == "V" else np.asfortranarray(np.einsum("ijn,kjn->ikn", F, F))
    lam = sim.isochoric_pstretch(X, tensor)
    assert lam.shape == (3, N)
    for i in range(N):
        np.testing.assert_array_equal(lam[:, i], np.ravel(sim.isochoric_pstretch(_slice(X, i), tensor)))
    # Voigt columns (6, N), as the single (6,) vector
    Xv = np.asfortranarray(np.stack([np.ravel(sim.t2v_strain(_slice(X, i))) for i in range(N)], axis=1))
    np.testing.assert_array_equal(sim.isochoric_pstretch(Xv, tensor), lam)


def test_isochoric_invariants_batch():
    F = _F()
    b = np.asfortranarray(np.einsum("ijn,kjn->ikn", F, F))
    inv = sim.isochoric_invariants(b)
    assert inv.shape == (3, N)
    for i in range(N):
        np.testing.assert_array_equal(inv[:, i], np.ravel(sim.isochoric_invariants(_slice(b, i))))
    bv = np.asfortranarray(np.stack([np.ravel(sim.t2v_strain(_slice(b, i))) for i in range(N)], axis=1))
    np.testing.assert_array_equal(sim.isochoric_invariants(bv), inv)
    # the same J for one point and for a batch
    np.testing.assert_array_equal(sim.isochoric_invariants(b, 1.3)[:, 0],
                                  np.ravel(sim.isochoric_invariants(_slice(b, 0), 1.3)))


def test_isochoric_pstretch_single_point_reads_J():
    """A single point honours J as the batch does (it used to recompute it from the tensor)."""
    b = np.asfortranarray(np.einsum("ijn,kjn->ikn", _F(n=2), _F(n=2)))
    for X in (_slice(b, 0), np.ravel(sim.t2v_strain(_slice(b, 0)))):
        single = np.ravel(sim.isochoric_pstretch(X, "b", 1.3))
        np.testing.assert_array_equal(single, sim.isochoric_pstretch(b, "b", 1.3)[:, 0])
        assert not np.array_equal(single, np.ravel(sim.isochoric_pstretch(X, "b")))
    with pytest.raises(ValueError, match="input string"):
        sim.isochoric_pstretch(_slice(b, 0), "C")


@pytest.mark.parametrize("key", ["Cauchy2PKI", "Cauchy2PKII", "PKII2Cauchy", "Kirchoff2LogStressU"])
def test_stress_convert_tensor_batch(key):
    F = _F()
    rng = np.random.default_rng(5)
    s = 100.0 * rng.standard_normal((3, 3, N))
    s = np.asfortranarray(0.5 * (s + s.transpose(1, 0, 2)))
    out = sim.stress_convert(s, F, key)
    assert out.shape == (3, 3, N)
    for i in range(N):
        np.testing.assert_array_equal(out[:, :, i], sim.stress_convert(_slice(s, i), _slice(F, i), key))
    # F given as (N, 3, 3)
    np.testing.assert_array_equal(sim.stress_convert(s, np.ascontiguousarray(F.transpose(2, 0, 1)), key), out)


def test_batch_shape_errors():
    with pytest.raises(Exception):
        sim.VR_decomposition(np.zeros((3, 2, 4)))
    with pytest.raises(Exception):
        sim.stress_convert(np.zeros((3, 3, 4)), _F(5), "Cauchy2PKI")
