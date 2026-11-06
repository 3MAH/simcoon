"""rotate_voigt: batch rotation of Voigt tensors by (3,3,N) rotation matrices.

It must match Rotation.from_matrix(DR).apply_<kind> to round-off, for one point and for
batches past the parallel cutoff, and give the same bits whatever n_threads."""

import numpy as np
import pytest

import simcoon as sim

_VECTOR_KINDS = ["stress", "strain"]
_TENSOR_KINDS = ["stiffness", "compliance", "strain_concentration", "stress_concentration"]


def _rotations(n, seed=0):
    """(3,3,n) F-ordered rotation matrices, as objective_rate returns DR."""
    mats = sim.Rotation.random(n, random_state=seed).as_matrix().reshape(n, 3, 3)
    return np.asfortranarray(mats.transpose(1, 2, 0))


def _values(kind, n, seed=1):
    rng = np.random.default_rng(seed)
    shape = (6, n) if kind in _VECTOR_KINDS else (6, 6, n)
    return np.asfortranarray(100.0 * rng.standard_normal(shape))


def _reference(values, DR, kind, active):
    """The scipy path fedoo used: one Rotation object for the batch, then apply_<kind>."""
    rot = sim.Rotation.from_matrix(DR.transpose(2, 0, 1) if DR.ndim == 3 else DR)
    return getattr(rot, "apply_" + kind)(values, active=active)


def _assert_roundoff(result, reference):
    scale = max(1.0, float(np.abs(reference).max()))
    np.testing.assert_allclose(result, reference, rtol=0, atol=1e-14 * scale)


@pytest.mark.parametrize("active", [True, False])
@pytest.mark.parametrize("kind", _VECTOR_KINDS + _TENSOR_KINDS)
@pytest.mark.parametrize("n", [1, 7, 1000])
def test_matches_rotation_apply(kind, n, active):
    DR = _rotations(n)
    values = _values(kind, n)
    result = sim.rotate_voigt(values, DR, kind, active)
    assert result.shape == values.shape
    _assert_roundoff(result, _reference(values, DR, kind, active))


@pytest.mark.parametrize("kind", _VECTOR_KINDS + _TENSOR_KINDS)
def test_single_tensor_single_rotation(kind):
    DR = _rotations(1)[:, :, 0]
    values = _values(kind, 1)[..., 0]
    result = sim.rotate_voigt(values, DR, kind)
    assert result.shape == values.shape
    _assert_roundoff(result, _reference(values, DR, kind, True))


@pytest.mark.parametrize("kind", ["stress", "stiffness"])
@pytest.mark.parametrize("one_rotation", [lambda DR: DR[:, :, 0], lambda DR: DR[:, :, :1]])
def test_one_rotation_for_every_tensor(kind, one_rotation):
    n = 250
    DR = _rotations(n)
    values = _values(kind, n)
    repeated = np.asfortranarray(np.repeat(DR[:, :, :1], n, axis=2))
    np.testing.assert_array_equal(sim.rotate_voigt(values, one_rotation(DR), kind),
                                  sim.rotate_voigt(values, repeated, kind))


@pytest.mark.parametrize("kind", ["strain", "stiffness"])
def test_same_bits_for_any_thread_count(kind):
    n = 5000
    DR = _rotations(n)
    values = _values(kind, n)
    serial = sim.rotate_voigt(values, DR, kind, n_threads=1)
    for n_threads in (0, 2, 4):
        np.testing.assert_array_equal(sim.rotate_voigt(values, DR, kind, n_threads=n_threads), serial)


@pytest.mark.parametrize("kind", ["stress", "compliance"])
def test_any_memory_layout(kind):
    n = 300
    DR = _rotations(n)
    values = _values(kind, n)
    expected = sim.rotate_voigt(values, DR, kind)
    # C-ordered copies, and DR as the (3,3,N) view of an (N,3,3) stack
    DR_view = np.ascontiguousarray(DR.transpose(2, 0, 1)).transpose(1, 2, 0)
    np.testing.assert_array_equal(
        sim.rotate_voigt(np.ascontiguousarray(values), DR_view, kind), expected)


def test_inputs_are_not_modified():
    DR = _rotations(200)
    values = _values("stiffness", 200)
    DR_copy, values_copy = DR.copy(), values.copy()
    sim.rotate_voigt(values, DR, "stiffness")
    np.testing.assert_array_equal(DR, DR_copy)
    np.testing.assert_array_equal(values, values_copy)


def test_objective_rate_increment():
    """The fedoo use: transport start strain and stress with the DR of objective_rate."""
    n = 400
    rng = np.random.default_rng(5)
    F0 = np.asfortranarray(np.eye(3)[:, :, None] + 0.05 * rng.standard_normal((3, 3, n)))
    F1 = np.asfortranarray(F0 + 0.02 * rng.standard_normal((3, 3, n)))
    _, DR, _ = sim.objective_rate("log_R", F0, F1, 1.0)
    for kind in ("strain", "stress"):
        values = _values(kind, n)
        _assert_roundoff(sim.rotate_voigt(values, DR, kind, n_threads=0),
                         _reference(values, DR, kind, True))


@pytest.mark.parametrize("values, DR, kind, match", [
    (np.zeros((6, 4)), np.zeros((3, 3, 4)), "tangent", "unknown kind"),
    (np.zeros((3, 4)), np.zeros((3, 3, 4)), "stress", r"shape \(6,\) or \(6, N\)"),
    (np.zeros((6, 6, 4)), np.zeros((3, 3, 4)), "strain", r"shape \(6,\) or \(6, N\)"),
    (np.zeros((6, 4)), np.zeros((3, 3, 4)), "stiffness", r"shape \(6, 6\) or \(6, 6, N\)"),
    (np.zeros((6, 4)), np.zeros((3, 4)), "stress", r"DR must have shape"),
    (np.zeros((6, 4)), np.zeros((3, 3, 3)), "stress", "expected 1 or 4"),
    (np.zeros(6), np.zeros((3, 3, 2)), "stress", "expected 1 or 1"),
])
def test_rejects_bad_shapes(values, DR, kind, match):
    with pytest.raises(ValueError, match=match):
        sim.rotate_voigt(values, DR, kind)
