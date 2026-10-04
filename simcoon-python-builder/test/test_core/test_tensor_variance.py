"""The variance tag of Tensor2 / Tensor4 components (index raising and lowering).

The variance is a convention of the representation, defaulted from the type
(stress contravariant, strain covariant, stiffness (contra, contra), compliance
(co, co), concentrations mixed) and changed with to_variance(): a retag in the lab
or an orthonormal basis, a contraction with the metric in a natural basis.
"""

import numpy as np
import pytest

import simcoon as sim

F = np.array([[1.20, 0.15, -0.05], [0.10, 0.90, 0.20], [0.00, -0.10, 1.10]])
CONV = sim.Basis.from_F(F)
G, GINV = CONV.metric, CONV.inverse_metric
SIG = sim.Tensor2.stress(np.array([100.0, -40.0, 25.0, 30.0, -12.0, 8.0]))
EPS = sim.Tensor2.strain(np.array([0.010, -0.004, 0.002, 0.006, -0.003, 0.001]))
L = sim.Tensor4.stiffness(sim.L_cubic([180000.0, 110000.0, 75000.0], "Cii")).rotate(
    sim.Rotation.from_euler("zxz", [-60.0, 10.0, 35.0], degrees=True))
CO, CONTRA = "covariant", "contravariant"


def close(a, b, tol=1e-10):
    np.testing.assert_allclose(np.asarray(a), np.asarray(b), rtol=tol, atol=tol * 1e3)


def test_defaults_follow_the_type():
    assert SIG.variance == CONTRA and EPS.variance == CO
    assert sim.Tensor2.from_voigt(SIG.voigt, "symmetric").variance is None
    assert sim.Tensor2.from_mat(F, "none").variance is None
    assert L.variance == (CONTRA, CONTRA) and L.inverse().variance == (CO, CO)
    assert sim.Tensor4.strain_concentration(np.eye(6)).variance == (CO, CONTRA)
    assert sim.Tensor4.stress_concentration(np.eye(6)).variance == (CONTRA, CO)
    assert "variance" not in repr(SIG) and "covariant" in repr(SIG.to_variance(CO))


def test_in_the_lab_to_variance_only_retags():
    low = SIG.to_variance(CO)
    assert np.array_equal(low.voigt, SIG.voigt) and low.variance == CO and low.type == "stress"
    assert low.to_variance(CONTRA) == SIG and (low == SIG) is False
    assert SIG.to_variance(CONTRA) is SIG
    with pytest.raises(ValueError, match="variance must be"):
        SIG.to_variance("mixed")
    with pytest.raises(ValueError, match="no variance"):
        sim.Tensor2.from_mat(F, "none").to_variance(CO)


def test_lowering_and_raising_with_the_metric():
    tau = SIG.with_basis(CONV)                                  # tau^ij
    low = tau.to_variance(CO)                                   # tau_ij = g_ik tau^kl g_lj
    close(low.mat, G @ tau.mat @ G)
    close(low.to_variance(CONTRA).mat, tau.mat)
    for name in ("trace", "mises", "norm", "det"):
        close(getattr(low, name)(), getattr(tau, name)())       # the tensor did not change
    close(low.eigvals(), tau.eigvals())
    close(low.to_basis(None).voigt, tau.to_basis(None).voigt)   # lab components agree
    close(low.dev().to_basis(None).voigt, tau.dev().to_basis(None).voigt)
    e = EPS.with_basis(CONV)
    close(e.to_variance(CONTRA).mat, GINV @ e.mat @ GINV)
    # the identity: g^ij lowered is g_ij
    close(sim.Tensor2.identity("stress", basis=CONV).to_variance(CO).mat, G)


def test_contractions_use_the_tags():
    tau, e = SIG.with_basis(CONV), EPS.with_basis(CONV)
    low = tau.to_variance(CO)
    close(tau % e, SIG % EPS)                                   # dual: plain sum
    close(low % tau, tau % tau)                                 # dual now: equals the metric form
    close(low % e.to_variance(CONTRA), SIG % EPS)               # both tags flipped
    close(low % low, tau % tau)                                 # same variance: metric on both
    with pytest.raises(ValueError, match="Mixed variance"):
        tau + low
    mixed_lab = SIG + SIG.to_variance(CO)                       # harmless in the lab: left tag wins
    assert mixed_lab.variance == CONTRA


def test_tensor4_pairs():
    Lc = L.to_basis(CONV)
    e = EPS.to_basis(CONV)
    s_ref = (L @ EPS)
    low = Lc.to_variance((CO, CO))
    assert low.variance == (CO, CO) and low.type == "stiffness"
    # lowered stiffness contracts a contravariant strain and yields a covariant stress
    with pytest.raises(ValueError, match="to_variance"):
        low @ e
    s_low = low @ e.to_variance(CONTRA)
    assert s_low.variance == CO
    close(s_low.to_basis(None).voigt, s_ref.voigt)
    close(low.to_basis(None).mat, L.mat)                        # lab components of the same tensor
    close(low.to_variance((CONTRA, CONTRA)).mat, Lc.mat)
    half = Lc.to_variance((CO, CONTRA))                         # one pair only
    close((half @ e).to_basis(None).voigt, s_ref.voigt)
    assert low.inverse().variance == (CONTRA, CONTRA) and half.inverse().variance == (CO, CONTRA)
    close(low.inverse().to_basis(None).mat, L.inverse().mat)
    A = sim.Tensor4.strain_concentration(sim.A_R(F)).to_basis(CONV)
    A_flipped = A.to_variance((CONTRA, CO))
    close((A_flipped @ e.to_variance(CONTRA)).to_basis(None).voigt,
          (sim.Tensor4.strain_concentration(sim.A_R(F)) @ EPS).voigt)
    with pytest.raises(ValueError, match="pair"):
        L.to_variance(CO)


def test_symmetric_type_declares_its_variance():
    sym = sim.Tensor2.from_voigt(SIG.voigt, "symmetric")
    with pytest.raises(ValueError, match="to_variance"):
        sym.with_basis(CONV)
    declared = sym.to_variance(CONTRA)
    assert declared.variance == CONTRA and np.array_equal(declared.voigt, sym.voigt)
    close(declared.with_basis(CONV).trace(), SIG.with_basis(CONV).trace())
    close(declared.push_forward(F, metric=False).voigt, SIG.push_forward(F, metric=False).voigt)
    assert sim.dyadic(SIG.to_variance(CO), EPS).variance == (CO, CO)
    # a lowered stress or stiffness in the lab transports as covariant components
    low = SIG.to_variance(CO).push_forward(F, metric=False)
    Fi = np.linalg.inv(F)
    close(low.mat, Fi.T @ SIG.mat @ Fi)
    close(SIG.to_variance(CO).push_forward(F).mat, Fi.T @ SIG.mat @ Fi / np.linalg.det(F))
    close(low.pull_back(F, metric=False).voigt, SIG.voigt)
    # a covariantly tagged stiffness pushes as covariant components (F^-T on the four
    # indices): the same numbers a compliance kernel would produce. Lowering in the lab
    # and pushing is not the push of the contravariant tensor: the variance is part of
    # what "push-forward" means.
    L_low = L.to_variance((CO, CO))
    as_compliance = sim.Tensor4.from_mandel(L.mandel, "compliance")   # same tensor, covariant kernel
    close(L_low.push_forward(F, metric=False).mandel,
          as_compliance.push_forward(F, metric=False).mandel)
    close(L_low.push_forward(F, metric=False).pull_back(F, metric=False).mat, L.mat)
    assert sim.Tensor2.from_list([SIG, SIG]).variance == CONTRA
    with pytest.raises(ValueError, match="Mixed variance"):
        sim.Tensor2.from_list([SIG, SIG.to_variance(CO)])
