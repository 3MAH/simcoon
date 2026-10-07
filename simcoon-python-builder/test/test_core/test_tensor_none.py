"""Tensor2 type "none": a non-symmetric 3x3 stored as 9 components.

F, R, L or PK1 do not fit the 6 Voigt components of the symmetric types; the
type "none" (no Voigt convention, the C++ Tensor2Type::none) holds them as lab-lab
components, with
no Voigt vector, no variance and no transport.
"""

import numpy as np
import pytest

import simcoon as sim

F = np.array([[1.20, 0.15, -0.05], [0.10, 0.90, 0.20], [0.00, -0.10, 1.10]])
R = sim.Rotation.from_euler("zxz", [25.0, 40.0, -15.0], degrees=True)


def test_none_holds_a_non_symmetric_matrix_exactly():
    t = sim.Tensor2.from_mat(F, "none")
    assert t.type == "none" and np.asarray(t).shape == (9,)
    assert np.array_equal(t.mat, F)
    with pytest.raises(ValueError, match="symmetric tensors only"):
        sim.Tensor2.from_mat(F, "stress")
    for bad in (lambda: t.voigt, lambda: t.mandel, lambda: sim.Tensor2.from_voigt(np.zeros(6), "none")):
        with pytest.raises(ValueError, match="no Voigt"):
            bad()


def test_none_invariants_and_arithmetic():
    t = sim.Tensor2.from_mat(F, "none")
    assert np.isclose(t.trace(), np.trace(F))
    assert np.isclose(t.det(), np.linalg.det(F))
    assert np.isclose(t.norm(), np.linalg.norm(F))
    np.testing.assert_allclose(t.dev().mat, F - np.trace(F) / 3 * np.eye(3))
    np.testing.assert_allclose(np.sort_complex(np.atleast_1d(t.eigvals())),
                               np.sort_complex(np.linalg.eigvals(F)))     # complex pair here
    sym = sim.Tensor2.from_mat(F @ F.T, "none")
    np.testing.assert_allclose(sym.eigvals(), np.linalg.eigvalsh(F @ F.T))  # real spectrum: real, sorted
    np.testing.assert_allclose((2 * t - t).mat, F)
    assert np.isclose(t % t, np.sum(F * F))
    assert not t.is_symmetric()
    np.testing.assert_array_equal(sim.Tensor2.identity("none").mat, np.eye(3))
    with pytest.raises(ValueError, match="Mises"):
        t.mises()


def test_none_rotation_and_orthonormal_basis():
    t = sim.Tensor2.from_mat(F, "none")
    Q = R.as_matrix()
    np.testing.assert_allclose(t.rotate(R).mat, Q @ F @ Q.T, atol=1e-13)
    passive = t.rotate(R, active=False)
    np.testing.assert_allclose(passive.mat, Q.T @ F @ Q, atol=1e-13)
    assert passive.basis.orthonormal
    np.testing.assert_allclose(passive.to_basis(None).mat, F, atol=1e-13)
    n = 5
    Fb = np.eye(3) + 0.1 * np.random.default_rng(0).standard_normal((n, 3, 3))
    tb = sim.Tensor2.from_mat(Fb, "none")
    rots = sim.Rotation.random(n, random_state=3)
    for i in range(n):
        np.testing.assert_allclose(tb.rotate(rots)[i].mat, tb[i].rotate(rots[i]).mat, atol=1e-13)
    np.testing.assert_allclose(tb.rotate(R).mat, Q @ Fb @ Q.T, atol=1e-13)


def test_none_has_no_variance_and_no_transport():
    t = sim.Tensor2.from_mat(F, "none")
    for bad in (lambda: t.push_forward(F), lambda: t.pull_back(F),
                lambda: t.with_basis(sim.Basis.from_F(F)), lambda: t.to_basis(sim.Basis.from_F(F))):
        with pytest.raises(ValueError, match="no variance"):
            bad()
    L = sim.Tensor4.stiffness(sim.L_iso([70000.0, 0.3], "Enu"))
    with pytest.raises(ValueError, match="symmetric Tensor2"):
        L @ t
    with pytest.raises(ValueError, match="no Voigt"):
        sim.dyadic(t, t)


def test_none_batches():
    n = 4
    Fb = np.eye(3) + 0.1 * np.random.default_rng(1).standard_normal((n, 3, 3))
    tb = sim.Tensor2.from_mat(Fb, "none")
    assert len(tb) == 4 and np.asarray(tb).shape == (n, 9)
    np.testing.assert_array_equal(tb[2].mat, Fb[2])
    np.testing.assert_array_equal(sim.Tensor2.from_list([tb[i] for i in range(n)]).mat, Fb)
    np.testing.assert_array_equal(sim.Tensor2.concatenate([tb[:2], tb[2:]]).mat, Fb)
    cols = sim.Tensor2.from_columns(np.asarray(tb).T, "none")
    np.testing.assert_array_equal(cols.mat, Fb)
    np.testing.assert_allclose(tb.det(), np.linalg.det(Fb))
    with pytest.raises(ValueError, match="Mixed type"):
        sim.Tensor2.from_list([tb[0], sim.Tensor2.stress(np.zeros(6))])
