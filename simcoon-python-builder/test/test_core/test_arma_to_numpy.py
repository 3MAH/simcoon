"""Armadillo -> NumPy output conversion (arma_to_numpy.hpp) is an index identity.

The converter changes no mechanics, so it can only break a result by transposing a matrix,
permuting cube axes or handing out a dangling buffer. Symmetric inputs hide the first two,
hence the non-symmetric probes (simple shear) and slice-distinct batches below. Bitwise
parity with the carma build it replaced (Lt_convert, A_R/A_F, dR_drotvec, sim.umat EPICP)
was checked once, on one platform, before the switch.
"""

import gc

import numpy as np
import pytest

import simcoon as sim

GAMMA = 0.5
F_SHEAR = np.array([[1.0, GAMMA, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
EPICP = [200000.0, 0.3, 0.0, 300.0, 1000.0, 0.5]


def _owned_by_capsule(a):
    assert type(a.base).__name__ == "PyCapsule"
    assert not a.flags.owndata
    assert a.flags.writeable
    assert a.flags.f_contiguous


# --- layout -------------------------------------------------------------------------------

def test_mat_col_cube_layout():
    m = sim.L_iso([70000.0, 0.3], "Enu")
    assert m.shape == (6, 6)
    _owned_by_capsule(m)
    v = sim.t2v_stress(np.eye(3))
    assert v.shape == (6, 1)  # Col stays (n,1): callers .ravel() it
    _owned_by_capsule(v)
    c = sim.dR_drotvec(np.array([0.1, 0.2, 0.3]))
    assert c.shape == (3, 3, 3)
    _owned_by_capsule(c)


def test_small_buffers_survive_heap_churn():
    # n_elem <= 16 lives inside the Armadillo object (mem_local): the array must point
    # into the heap object the capsule owns, never into a dead stack frame.
    keep = [sim.G_UdX(F_SHEAR), sim.t2v_strain(F_SHEAR), sim.dR_drotvec(np.array([0.1, 0.2, 0.3]))]
    want = [k.copy() for k in keep]
    for _ in range(10_000):
        sim.t2v_stress(np.eye(3))
    gc.collect()
    _ = [np.ones((4, 4)) for _ in range(100_000)]
    for k, w in zip(keep, want):
        np.testing.assert_array_equal(k, w)


# --- transposition probes (exact, simple shear) -------------------------------------------

def test_3x3_outputs_are_not_transposed():
    g = sim.G_UdX(F_SHEAR)
    assert g[0, 1] == GAMMA and g[1, 0] == 0.0
    assert sim.R_Cauchy_Green(F_SHEAR)[0, 0] == pytest.approx(1.0)
    assert sim.L_Cauchy_Green(F_SHEAR)[0, 0] == pytest.approx(1.0 + GAMMA**2)
    L = sim.finite_L(np.eye(3), F_SHEAR, 1.0)
    assert L[0, 1] == pytest.approx(GAMMA) and L[1, 0] == 0.0
    W = sim.finite_W(np.eye(3), F_SHEAR, 1.0)
    assert W[0, 1] == pytest.approx(GAMMA / 2) and W[1, 0] == pytest.approx(-GAMMA / 2)
    R, U = sim.RU_decomposition(F_SHEAR)
    assert R[0, 1] == pytest.approx(GAMMA / np.sqrt(4 + GAMMA**2))
    np.testing.assert_allclose(R @ U, F_SHEAR, atol=1e-14)
    np.testing.assert_allclose(U, U.T, atol=1e-14)


def test_6x6_output_is_not_transposed():
    nu = 0.3
    S = sim.Eshelby_cylinder(nu)
    assert S[1, 0] == pytest.approx(nu / (2 * (1 - nu)))
    assert S[0, 1] == 0.0


def test_voigt_slots_and_factors():
    T = np.array([[1.0, 4.0, 5.0], [4.0, 2.0, 6.0], [5.0, 6.0, 3.0]])
    np.testing.assert_array_equal(sim.t2v_stress(T).ravel(), [1, 2, 3, 4, 5, 6])
    np.testing.assert_array_equal(sim.t2v_strain(T).ravel(), [1, 2, 3, 8, 10, 12])


def test_cube_slices_are_not_permuted():
    n = 5
    F = np.repeat(np.eye(3)[:, :, None], n, axis=2)
    F[0, 1, :] = 0.1 * (np.arange(n) + 1)
    F = np.asfortranarray(F)
    batch = sim.Log_strain(F)
    for k in range(n):
        np.testing.assert_allclose(batch[:, :, k], sim.Log_strain(F[:, :, k]), atol=1e-15)
    I = np.asfortranarray(np.repeat(np.eye(3)[:, :, None], n, axis=2))
    D, DR, Omega = sim.objective_rate("green_naghdi", I, F, 1.0)
    for k in range(n):
        d, dr, om = sim.objective_rate("green_naghdi", np.eye(3), F[:, :, k], 1.0)
        np.testing.assert_array_equal(D[:, :, k], d)
        np.testing.assert_array_equal(DR[:, :, k], dr)
        np.testing.assert_array_equal(Omega[:, :, k], om)


# --- fedoo pattern ------------------------------------------------------------------------

def _umat_batch(n):
    col = lambda a: np.asfortranarray(np.tile(np.asarray(a, float).reshape(-1, 1), (1, n)))
    de = np.asfortranarray(1e-3 * np.vstack([np.linspace(1, 4, n), -np.linspace(0.3, 1, n),
                                             -np.linspace(0.3, 1, n), np.linspace(0.5, 2, n),
                                             np.linspace(0, 1, n), np.linspace(1, 0, n)]))
    I = np.asfortranarray(np.repeat(np.eye(3)[:, :, None], n, axis=2))
    args = [col(np.zeros(6)), de, I, I.copy(order="F"), col(np.zeros(6)), I.copy(order="F"),
            col(EPICP), col([290.0] + [0.0] * 7)]
    return args, col(np.zeros(4))


def test_umat_outputs_own_their_memory():
    args, wm = _umat_batch(1000)
    out = sim.umat("EPICP", *args, 0.5, 1.0, wm, n_threads=4)
    for o in out:
        _owned_by_capsule(o)
        for a in args + [wm]:
            assert not np.shares_memory(o, a)
    want = [o.copy() for o in out]
    del args, wm
    gc.collect()
    for o, w in zip(out, want):
        np.testing.assert_array_equal(o, w)


def test_chained_calls_leave_inputs_untouched():
    # fedoo feeds the outputs of one increment back as the inputs of the next.
    args, wm = _umat_batch(8)
    sigma, statev, wm, Lt = sim.umat("EPICP", *args, 0.5, 1.0, wm, n_threads=1)
    snapshot = [sigma.copy(), statev.copy(), wm.copy()]
    args[4], args[7] = sigma, statev
    sim.umat("EPICP", *args, 1.0, 1.0, wm, n_threads=1)
    for a, s in zip([sigma, statev, wm], snapshot):
        np.testing.assert_array_equal(a, s)
