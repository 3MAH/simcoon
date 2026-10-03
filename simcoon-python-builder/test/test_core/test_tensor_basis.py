"""The basis a Tensor2 / Tensor4 is written in (simcoon.Basis).

``basis is None`` is the lab orthonormal basis. An orthonormal Basis carries a
rotation (identity metric, every formula as in the lab); a natural Basis carries
any three vectors -- typically the convected basis of F -- and its metric enters
the invariants. The type tag gives the variance: stress contravariant, strain
covariant. Closed forms are checked on simple shear F = I + gamma e1 (x) e2.
"""

import copy
import pickle

import numpy as np
import pytest

import simcoon as sim
from simcoon.tensor import _voigt_operators

GAMMA = 0.5
F_SHEAR = np.array([[1.0, GAMMA, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
F_GEN = np.array([[1.20, 0.15, -0.05], [0.10, 0.90, 0.20], [0.00, -0.10, 1.10]])
R1 = sim.Rotation.from_euler("zxz", [25.0, 40.0, -15.0], degrees=True)
R2 = sim.Rotation.from_euler("zxz", [-60.0, 10.0, 35.0], degrees=True)

SIG = sim.Tensor2.stress(np.array([100.0, -40.0, 25.0, 30.0, -12.0, 8.0]))
EPS = sim.Tensor2.strain(np.array([0.010, -0.004, 0.002, 0.006, -0.003, 0.001]))
L_ISO = sim.Tensor4.stiffness(sim.L_iso([70000.0, 0.3], "Enu"))
L_ANI = sim.Tensor4.stiffness(sim.L_cubic([180000.0, 110000.0, 75000.0], "Cii")).rotate(R2)


def close(a, b, tol=1e-10):
    np.testing.assert_allclose(np.asarray(a), np.asarray(b), rtol=tol, atol=tol * 1e3)


# --- lab default and orthonormal frames ----------------------------------------------------

def test_lab_is_the_default_and_free():
    assert SIG.basis is None and L_ISO.basis is None
    assert (SIG + SIG).basis is None and (L_ISO @ EPS).basis is None
    assert SIG.rotate(R1).basis is None                 # active: another tensor, still lab
    close(SIG.rotate(R1).voigt, R1.apply_stress(SIG.voigt).ravel())
    assert SIG.to_basis(None) is SIG


def test_passive_rotation_ships_its_frame():
    sig_r = SIG.rotate(R1, active=False)
    assert sig_r.basis.orthonormal and sig_r.basis.rotation.equals(R1)
    close(sig_r.voigt, R1.apply_stress(SIG.voigt, active=False).ravel())   # numbers as before
    close(sig_r.to_basis(None).voigt, SIG.voigt)                           # same tensor
    twice = sig_r.rotate(R2, active=False)                                 # R2 read in the R1 frame
    assert twice.basis.rotation.equals(R1 * R2)
    close(twice.to_basis(None).voigt, SIG.voigt)
    L_r = L_ANI.rotate(R1, active=False)
    close(L_r.to_basis(None).mat, L_ANI.mat)


def test_to_basis_orthonormal_matches_passive_rotation():
    b = sim.Basis(rotation=R1, name="material")
    close(SIG.to_basis(b).voigt, SIG.rotate(R1, active=False).voigt)
    close(L_ANI.to_basis(b).mat, L_ANI.rotate(R1, active=False).mat)
    b2 = sim.Basis(rotation=R2)
    close(SIG.to_basis(b).to_basis(b2).voigt, SIG.to_basis(b2).voigt)
    assert "material" in repr(SIG.to_basis(b))


def test_active_rotation_of_a_framed_tensor_turns_its_basis():
    b = sim.Basis(rotation=R1)
    sig_b = SIG.to_basis(b)
    turned = sig_b.rotate(R2)
    assert np.array_equal(turned.voigt, sig_b.voigt)            # components kept
    assert turned.basis.rotation.equals(R2 * R1)
    close(turned.to_basis(None).voigt, SIG.rotate(R2).voigt)    # the lab rotation of the tensor


def test_mixed_basis_raises():
    b = sim.Basis(rotation=R1)
    sig_b, eps_b, L_b = SIG.to_basis(b), EPS.to_basis(b), L_ISO.to_basis(b)
    for op in (lambda: SIG + sig_b, lambda: sig_b - SIG, lambda: SIG % sig_b,
               lambda: L_ISO @ eps_b, lambda: L_b @ EPS, lambda: np.add(SIG, sig_b),
               lambda: sim.dyadic(SIG, sig_b), lambda: sim.double_contract(SIG, sig_b),
               lambda: sim.Tensor2.concatenate([SIG, sig_b])):
        with pytest.raises(ValueError, match="Mixed basis"):
            op()
    assert (SIG == sig_b) is False
    same = SIG.to_basis(sim.Basis(rotation=R1))                 # another object, equal basis
    close((sig_b + same).voigt, 2 * sig_b.voigt)


def test_operations_propagate_the_basis_object():
    b = sim.Basis(rotation=R1)
    s, e, L = SIG.to_basis(b), EPS.to_basis(b), L_ANI.to_basis(b)
    results = [s + s, s - s, -s, 2.0 * s, s / 3.0, np.add(s, s), s.dev(),
               L @ e, L * e, L.inverse(), L + L, 0.5 * L,
               sim.dyadic(s, s), sim.auto_dyadic(s), sim.sym_dyadic(s, s), sim.auto_sym_dyadic(s),
               sim.Tensor2.from_tensor(s, 4), sim.Tensor2.from_tensor(s, 4)[1],
               sim.Tensor2.concatenate([s, s]), sim.Tensor2.from_list([s, s]),
               sim.Tensor4.from_list([L, L])]
    for r in results:
        assert r.basis is b
    close((L @ e).to_basis(None).voigt, (L_ANI @ EPS).voigt)
    close(L.inverse().to_basis(None).mat, L_ANI.inverse().mat)
    for value, lab in [(s.trace(), SIG.trace()), (s.mises(), SIG.mises()), (s.norm(), SIG.norm()),
                       (s.det(), SIG.det()), (s % e, SIG % EPS), (e.mises(), EPS.mises())]:
        close(value, lab)
    close(s.eigvals(), np.linalg.eigvalsh(SIG.mat))


def test_one_basis_for_a_batch_and_per_tensor_bases():
    n = 6
    rng = np.random.default_rng(3)
    batch = sim.Tensor2.stress(rng.standard_normal((n, 6)))
    shared = sim.Basis(rotation=R1)
    in_shared = batch.to_basis(shared)
    assert in_shared.basis is shared and in_shared[2].basis is shared       # held once
    for i in range(n):
        close(in_shared[i].voigt, batch[i].to_basis(shared).voigt)
    rots = sim.Rotation.random(n, random_state=5)
    per_point = batch.to_basis(sim.Basis(rotation=rots))
    assert len(per_point.basis) == n
    for i in range(n):
        assert per_point[i].basis.single and per_point[i].basis.rotation.equals(rots[i])
        close(per_point[i].voigt, batch[i].rotate(rots[i], active=False).voigt)
    close(per_point.to_basis(None).voigt, batch.voigt)
    with pytest.raises(ValueError, match="batch size"):
        batch.with_basis(sim.Basis(rotation=rots[:3]))
    with pytest.raises(ValueError, match="single"):
        SIG.with_basis(sim.Basis(rotation=rots))


def test_stacking_merges_bases():
    b1, b2 = sim.Basis(rotation=R1), sim.Basis(rotation=R2)
    stacked = sim.Tensor2.from_list([SIG.to_basis(b1), SIG.to_basis(b2)])
    assert len(stacked.basis) == 2 and stacked.basis.orthonormal
    close(stacked.to_basis(None).voigt, np.vstack([SIG.voigt, SIG.voigt]))
    with pytest.raises(ValueError, match="Mixed basis"):
        sim.Tensor2.from_list([SIG, SIG.to_basis(b1)])


def test_pickle_and_copy_keep_the_basis():
    s = SIG.to_basis(sim.Basis(rotation=R1, name="material"))
    for clone in (pickle.loads(pickle.dumps(s)), copy.deepcopy(s)):
        assert clone.basis.name == "material" and clone.basis.equals(s.basis)
        assert clone == s


# --- natural basis: the metric ------------------------------------------------------------------

def test_voigt_operators_are_dual():
    P_sharp, P_flat = _voigt_operators(F_GEN)
    close(P_sharp.T @ P_flat, np.eye(6), 1e-13)
    close(P_sharp @ SIG.voigt, sim.Tensor2.stress(F_GEN @ SIG.mat @ F_GEN.T).voigt, 1e-13)
    Fi = np.linalg.inv(F_GEN)
    close(P_flat @ EPS.voigt, sim.Tensor2.strain(Fi.T @ EPS.mat @ Fi).voigt, 1e-13)


def test_convected_stress_on_simple_shear():
    s = 80.0
    S = sim.Tensor2.stress(np.array([0.0, s, 0.0, 0.0, 0.0, 0.0]))        # PK2 = s e2 (x) e2
    conv = sim.Basis.from_F(F_SHEAR)
    tau = S.with_basis(conv)                                               # tau^ij = S^IJ
    assert np.array_equal(tau.voigt, S.voigt) and not tau.basis.orthonormal
    close(tau.basis.metric, F_SHEAR.T @ F_SHEAR)                           # g = C
    close(tau.trace(), s * (1 + GAMMA**2))                                 # tr = g_ij tau^ij
    close(tau.mises(), s * (1 + GAMMA**2))                                 # uniaxial along g_2
    close(tau.norm(), s * (1 + GAMMA**2))
    tau_lab = tau.to_basis(None)
    close(tau_lab.voigt, S.push_forward(F_SHEAR, metric=False).voigt)
    close(tau.det(), tau_lab.det())
    close(tau.eigvals(), np.linalg.eigvalsh(tau_lab.mat))
    close(tau.dev().to_basis(None).voigt, tau_lab.dev().voigt)
    close(tau.dev().trace(), 0.0)


def test_convected_strain_on_simple_shear():
    C = F_SHEAR.T @ F_SHEAR
    E = sim.Tensor2.strain(0.5 * (C - np.eye(3)))                          # Green-Lagrange
    e = E.with_basis(sim.Basis.from_F(F_SHEAR))                            # e_ij = E_IJ
    almansi = e.to_basis(None).mat
    close(almansi, 0.5 * (np.eye(3) - np.linalg.inv(F_SHEAR @ F_SHEAR.T)))
    close(almansi[0, 1], GAMMA / 2)
    close(almansi[1, 1], -GAMMA**2 / 2)
    close(e.to_basis(None).voigt, E.push_forward(F_SHEAR).voigt)
    close(e.trace(), np.trace(almansi))                                    # tr = g^ij e_ij
    close(e.mises(), sim.Tensor2.strain(almansi).mises())


def test_dual_contraction_needs_no_metric_same_variance_does():
    conv = sim.Basis.from_F(F_GEN)
    tau, e = SIG.with_basis(conv), EPS.with_basis(conv)
    close(tau % e, SIG % EPS)                                              # sigma^ij eps_ij
    close(tau % e, tau.to_basis(None) % e.to_basis(None))
    close(tau % tau, tau.to_basis(None) % tau.to_basis(None))              # metric on both indices
    close(e % e, e.to_basis(None) % e.to_basis(None))
    close(sim.double_contract(tau, tau), [tau.norm() ** 2])


def test_tensor4_in_a_natural_basis():
    conv = sim.Basis.from_F(F_GEN)
    # convected components of the spatial tensor = lab components of its pull-back
    close(L_ANI.to_basis(conv).mat, L_ANI.pull_back(F_GEN, metric=False).mat)
    M_ANI = L_ANI.inverse()
    close(M_ANI.to_basis(conv).mat, M_ANI.pull_back(F_GEN, metric=False).mat)
    L_c, eps_c = L_ANI.to_basis(conv), EPS.to_basis(conv)
    close((L_c @ eps_c).to_basis(None).voigt, (L_ANI @ EPS).voigt)
    close(L_c.inverse().to_basis(None).mat, M_ANI.mat)
    A = sim.Tensor4.strain_concentration(sim.A_R(F_GEN))
    close((A.to_basis(conv) @ eps_c).to_basis(None).voigt, (A @ EPS).voigt)
    B = sim.Tensor4.stress_concentration(np.linalg.inv(sim.A_R(F_GEN)).T)
    close((B.to_basis(conv) @ SIG.to_basis(conv)).to_basis(None).voigt, (B @ SIG).voigt)
    with pytest.raises(ValueError, match="natural basis"):
        L_c @ SIG.to_basis(conv)                                           # not dual


def test_metric_identity_and_projectors():
    conv = sim.Basis.from_F(F_GEN)
    g, ginv = conv.metric, conv.inverse_metric
    close(sim.Tensor2.identity("stress", basis=conv).mat, ginv)            # g^ij
    close(sim.Tensor2.identity("strain", basis=conv).mat, g)               # g_ij
    I4 = sim.Tensor4.identity("stiffness", basis=conv)
    # 1/2 (g^ik g^jl + g^il g^jk): spot values in Voigt, then its action (raises both indices)
    close(I4.mat[0, 0], ginv[0, 0] ** 2)
    close(I4.mat[0, 3], ginv[0, 0] * ginv[0, 1])
    close(I4.mat[3, 3], 0.5 * (ginv[0, 0] * ginv[1, 1] + ginv[0, 1] ** 2))
    eps_c = EPS.to_basis(conv)
    close((I4 @ eps_c).mat, ginv @ eps_c.mat @ ginv)
    close(sim.Tensor4.identity("strain_concentration", basis=conv).mat, np.eye(6))   # mixed
    P_vol = sim.Tensor4.volumetric("stiffness", basis=conv)
    P_dev = sim.Tensor4.deviatoric("stiffness", basis=conv)
    close((P_vol + P_dev).mat, I4.mat)
    K, mu = 70000.0 / (3 * (1 - 0.6)), 70000.0 / 2.6
    close(L_ISO.to_basis(conv).mat, (3 * K * P_vol + 2 * mu * P_dev).mat)


def test_reciprocal_stretch_and_polar():
    conv = sim.Basis.from_F(F_GEN)
    close(np.swapaxes(conv.matrix, -1, -2) @ conv.reciprocal, np.eye(3), 1e-13)   # g^i . g_j = delta
    close(conv.stretch @ conv.stretch, conv.metric, 1e-12)                           # U^2 = g
    R_F, U_F = sim.RU_decomposition(F_GEN)
    close(conv.stretch, U_F, 1e-12)
    close(conv.polar.matrix, R_F, 1e-12)
    close(conv.polar.matrix @ conv.stretch, F_GEN, 1e-12)                            # A = R U
    assert conv.polar.orthonormal and conv.polar is conv.polar                       # cached
    # the polar frame is the rotated orthonormal frame a convected tensor is often read in
    tau = SIG.with_basis(conv)
    close(tau.to_basis(conv.polar).voigt, tau.to_basis(None).rotate(sim.Rotation.from_matrix(R_F), active=False).voigt)
    ortho = sim.Basis(rotation=R1)
    assert ortho.polar is ortho
    close(ortho.stretch, np.eye(3))
    close(ortho.reciprocal, ortho.matrix)                                            # orthonormal: g^i = g_i
    n = 4
    batch = sim.Basis.from_F(np.eye(3) + 0.2 * np.random.default_rng(2).standard_normal((n, 3, 3)))
    close(batch.polar.matrix @ batch.stretch, batch.matrix, 1e-12)
    assert len(batch.polar) == n


def test_a_rotation_as_natural_basis_reduces_to_the_orthonormal_case():
    ortho, natural = sim.Basis(rotation=R1), sim.Basis(vectors=R1.as_matrix())
    assert ortho.orthonormal and not natural.orthonormal
    for t in (SIG, EPS):
        a, b = t.to_basis(ortho), t.to_basis(natural)
        close(a.voigt, b.voigt)
        for name in ("trace", "mises", "norm", "det"):
            close(getattr(a, name)(), getattr(b, name)())
    for t in (L_ANI, L_ANI.inverse(), sim.Tensor4.strain_concentration(sim.A_R(F_GEN))):
        close(t.to_basis(ortho).mat, t.to_basis(natural).mat)


def test_superposed_rotation_is_objective():
    tau = SIG.with_basis(sim.Basis.from_F(F_GEN))
    turned = tau.rotate(R1)                                    # g_i -> Q g_i, components kept
    assert np.array_equal(turned.voigt, tau.voigt)
    close(turned.basis.matrix, R1.as_matrix() @ F_GEN)
    close(turned.basis.metric, tau.basis.metric)
    for name in ("trace", "mises", "norm", "det"):
        close(getattr(turned, name)(), getattr(tau, name)())
    close(turned.to_basis(None).voigt, tau.to_basis(None).rotate(R1).voigt)
    passive = tau.rotate(R1, active=False)                     # same tensor, basis A.Q
    close(passive.to_basis(None).voigt, tau.to_basis(None).voigt)


# --- push-forward / pull-back as a change of basis ---------------------------------------------

@pytest.mark.parametrize("tensor", [SIG, EPS, L_ANI, L_ANI.inverse()],
                         ids=["stress", "strain", "stiffness", "compliance"])
@pytest.mark.parametrize("metric", [False, True])
def test_push_forward_is_a_change_of_basis(tensor, metric):
    eager = tensor.push_forward(F_GEN, metric=metric)          # lab tensor: lab components
    assert eager.basis is None
    if not metric:
        close(np.asarray(tensor.with_basis(sim.Basis.from_F(F_GEN)).to_basis(None)),
              np.asarray(eager))
    reference = tensor.with_basis(sim.Basis(rotation=sim.Rotation.identity()))
    lazy = reference.push_forward(F_GEN, metric=metric)        # own basis: convected, lazy
    close(lazy.basis.matrix, F_GEN)
    if not metric:
        assert np.array_equal(np.asarray(lazy), np.asarray(tensor))    # components untouched
    close(np.asarray(lazy.to_basis(None)), np.asarray(eager))
    back = lazy.pull_back(F_GEN, metric=metric)
    close(back.basis.matrix, np.eye(3))
    close(np.asarray(back), np.asarray(tensor))


def test_batched_convected_basis_aliases_F():
    n = 5
    rng = np.random.default_rng(11)
    F = np.eye(3) + 0.15 * rng.standard_normal((n, 3, 3))
    conv = sim.Basis.from_F(F)
    assert conv._A is F and len(conv) == n                     # no copy of the caller's array
    S = sim.Tensor2.stress(rng.standard_normal((n, 6)) * 50.0)
    tau = S.with_basis(conv)
    close(tau.to_basis(None).voigt, S.push_forward(F, metric=False).voigt)
    lab = tau.to_basis(None)
    close(tau.trace(), lab.trace())
    close(tau.mises(), lab.mises())
    close(tau.eigvals(), np.linalg.eigvalsh(lab.mat))
    close(tau[3].mises(), lab[3].mises())
    L = sim.Tensor4.from_tensor(L_ANI, n)
    close(L.with_basis(conv).to_basis(None).mat, L.push_forward(F, metric=False).mat)


def test_types_without_variance_are_refused_in_a_natural_basis():
    generic = sim.Tensor2.from_voigt(SIG.voigt, "generic")
    close(generic.to_basis(sim.Basis(rotation=R1)).to_basis(None).voigt, generic.voigt)   # fine
    with pytest.raises(ValueError, match="no variance"):
        generic.with_basis(sim.Basis.from_F(F_GEN))
    with pytest.raises(ValueError, match="no variance"):
        generic.to_basis(sim.Basis.from_F(F_GEN))
    with pytest.raises(TypeError):
        SIG.with_basis("material")                             # a name is not a basis


def test_numpy_view_is_components_in_the_own_basis():
    tau = SIG.with_basis(sim.Basis.from_F(F_GEN))
    assert np.array_equal(np.asarray(tau), SIG.voigt)
