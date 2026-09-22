"""The HOLZA potential: Gasser-Ogden-Holzapfel anisotropic hyperelasticity.

The fibre term is checked three ways that do not share machinery with it:
against a numerically differentiated energy, against NEOHC in the two limits
where it degenerates (no fibre stiffness, fully dispersed fibres), and against
the frame indifference it must satisfy.
"""

import numpy as np
import pytest
from scipy.linalg import expm

import simcoon as sim
from simcoon.modular import (
    HolzapfelElasticity,
    Damage,
    LinearIsotropicHardening,
    ModularMaterial,
    NeoHookeanElasticity,
    Plasticity,
    Viscoelasticity,
    VonMisesYield,
)

# a representative artery-like parameter set (MPa)
C10, K1, K2, KAPPA = 0.0354, 0.0107, 7.48, 1000.0

# two families at +/- 40 deg about e3, both in the e1-e2 plane
FIBRES = sim.Rotation.from_euler("zxz", [[0.0, 0.0, 40.0], [0.0, 0.0, -40.0]], degrees=True)


def _props(c10=C10, k1=K1, k2=K2, kappa_d=0.0, fibres=FIBRES, kappa=KAPPA):
    """The standalone HOLZA props, straight from the Python dataclass."""
    return HolzapfelElasticity(C10=c10, k1=k1, k2=k2, kappa_d=kappa_d,
                               fibres=fibres, kappa=kappa).potential_params()


def _umat(name, props, F1):
    n, nstatev = 1, 1
    z6 = lambda: np.zeros((6, n), order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    stress, sv, wm, Lt = sim.umat(
        name, z6(), z6(), eye,
        np.asarray(F1).reshape(3, 3, 1).copy(order="F"),
        z6(), eye,
        np.asfortranarray(np.asarray(props, dtype=float).reshape(-1, 1)),
        np.zeros((nstatev, n), order="F"), 0.0, 1.0,
        np.zeros((4, n), order="F"), n_threads=1)
    return stress[:, 0], Lt[:, :, 0]


def _energy(F, a0, kappa_d=0.0, c10=C10, k1=K1, k2=K2, kappa=KAPPA):
    """W(F) written out independently of the kernel, for the finite difference."""
    J = np.linalg.det(F)
    C_bar = J ** (-2.0 / 3.0) * (F.T @ F)
    I1_bar = np.trace(C_bar)
    W = c10 * (I1_bar - 3.0)
    for i in range(a0.shape[1]):
        a = a0[:, i]
        I4_star = kappa_d * I1_bar + (1.0 - 3.0 * kappa_d) * (a @ C_bar @ a)
        E = I4_star - 1.0
        if E > 0.0:
            W += k1 / (2.0 * k2) * (np.exp(k2 * E * E) - 1.0)
    return W + kappa * (J * np.log(J) - J + 1.0)   # the default "log" U(J)


def _F(eps):
    return expm(sim.v2t_strain(np.asarray(eps, dtype=float)))


GENERIC = [0.18, -0.06, -0.05, 0.05, -0.02, 0.03]


# ---------------------------------------------------------------- degeneracies

def test_no_fibre_stiffness_is_neo_hookean():
    """k1 = 0 removes the fibre term: HOLZA must BE NEOHC with mu = 2 C10."""
    F = _F(GENERIC)
    sig_h, Lt_h = _umat("HOLZA", _props(k1=0.0), F)
    sig_n, Lt_n = _umat("NEOHC", [2.0 * C10, KAPPA], F)
    np.testing.assert_allclose(sig_h, sig_n, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(Lt_h, Lt_n, rtol=1e-12, atol=1e-12)


def test_full_dispersion_is_isotropic():
    """kappa_d = 1/3 makes A = b_bar/3: the response cannot depend on the fibres."""
    F = _F(GENERIC)
    ref = _umat("HOLZA", _props(kappa_d=1.0 / 3.0), F)
    for angles in ([[0.0, 0.0, 0.0], [0.0, 90.0, 0.0]], [[0.0, 35.0, 12.0], [0.0, 70.0, 5.0]]):
        other = sim.Rotation.from_euler("zxz", angles, degrees=True)
        sig, Lt = _umat("HOLZA", _props(kappa_d=1.0 / 3.0, fibres=other), F)
        np.testing.assert_allclose(sig, ref[0], rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(Lt, ref[1], rtol=1e-11, atol=1e-11)


def test_fibres_carry_no_compression():
    """Where every I4* < 1 the fibre term is inactive and only the matrix answers."""
    # equibiaxial compression in the fibre plane: both families are shortened
    F = _F([-0.10, -0.10, 0.22, 0.0, 0.0, 0.0])
    a0 = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=0.0,
                             fibres=FIBRES, kappa=KAPPA).directions
    J = np.linalg.det(F)
    C_bar = J ** (-2.0 / 3.0) * (F.T @ F)
    assert all(a0[:, i] @ C_bar @ a0[:, i] < 1.0 for i in range(a0.shape[1]))

    sig_h, Lt_h = _umat("HOLZA", _props(), F)
    sig_n, Lt_n = _umat("NEOHC", [2.0 * C10, KAPPA], F)
    np.testing.assert_allclose(sig_h, sig_n, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(Lt_h, Lt_n, rtol=1e-12, atol=1e-12)


def test_fibres_stiffen_their_own_direction():
    """A stretch along a single fibre is stiffer than the same stretch across it."""
    along = sim.Rotation.identity()                              # a0 = e1
    across = sim.Rotation.from_euler("zxz", [0.0, 0.0, 90.0], degrees=True)   # a0 = e2
    F = _F([0.20, -0.09, -0.09, 0.0, 0.0, 0.0])                  # stretch along e1
    sig_along, _ = _umat("HOLZA", _props(fibres=along), F)
    sig_across, _ = _umat("HOLZA", _props(fibres=across), F)
    assert sig_along[0] > sig_across[0]


# ------------------------------------------------------------------- the model

@pytest.mark.parametrize("kappa_d", [0.0, 0.1, 1.0 / 3.0])
@pytest.mark.parametrize("state", [GENERIC, [0.25, 0.03, -0.10, 0.0, 0.04, 0.0]])
def test_stress_matches_finite_difference_of_the_energy(kappa_d, state):
    """sigma = J^-1 (dW/dF) F^T, with dW/dF taken numerically. No shared machinery."""
    F = _F(state)
    a0 = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=kappa_d,
                             fibres=FIBRES, kappa=KAPPA).directions

    d = 1e-6
    P = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Fp, Fm = F.copy(), F.copy()
            Fp[i, j] += d
            Fm[i, j] -= d
            P[i, j] = (_energy(Fp, kappa_d=kappa_d, a0=a0)
                       - _energy(Fm, kappa_d=kappa_d, a0=a0)) / (2.0 * d)
    sigma_fd = np.asarray(sim.t2v_stress((P @ F.T) / np.linalg.det(F))).ravel()

    sigma, _ = _umat("HOLZA", _props(kappa_d=kappa_d), F)
    np.testing.assert_allclose(sigma, sigma_fd, rtol=2e-6,
                               atol=2e-6 * max(1.0, np.abs(sigma).max()))


@pytest.mark.parametrize("kappa_d", [0.0, 0.1])
def test_fibre_tangent_matches_finite_difference(kappa_d):
    """The fibre TANGENT, with the fibres active.

    Every other Lt assertion in this file sits where the fibre term is identically zero
    (k1 = 0, compression) or compares HOLZA against HOLZA, so a systematic factor error in
    gamma_linear/gamma_quadratic would pass all of them. This one differences the stress
    the kernel returns against the tangent the kernel returns.

    Measured detection power at this state: a 1 % error in the fibre tangent exceeds the
    tolerance by ~4800x (kappa_d = 0) and ~2800x (kappa_d = 0.1).
    """
    # A COAXIAL state with the fibres on principal axes: ln V and the perturbation commute,
    # so the box tangent's normal block IS d(tau)/d(eps) and the difference needs no rate
    # conversion -- the same recipe as the volumetric check in test_modular_hyper.py.
    # A soft bulk modulus and a stretch well into the exponential put the fibre terms at
    # 56-96 % of the normal block; at the default kappa = 1000 they are 0.03 % of it and
    # the check would pass whatever the fibre tangent returned.
    fibres = sim.Rotation.from_euler("zxz", [[0.0, 0.0, 0.0], [0.0, 0.0, 90.0]], degrees=True)
    eps0 = np.array([0.30, 0.22, -0.52, 0.0, 0.0, 0.0])
    kappa = 5.0
    props = _props(kappa_d=kappa_d, fibres=fibres, kappa=kappa)

    a0 = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=kappa_d,
                             fibres=fibres, kappa=kappa).directions
    F0 = _F(eps0)
    Cb = np.linalg.det(F0) ** (-2.0 / 3.0) * (F0.T @ F0)
    active = [kappa_d * np.trace(Cb) + (1 - 3 * kappa_d) * (a0[:, i] @ Cb @ a0[:, i]) > 1.0
              for i in range(a0.shape[1])]
    assert all(active), f"the fibres must be ACTIVE for this to test anything: {active}"

    def tau_of(eps):
        sigma, _ = _umat("HOLZA", props, _F(eps))
        return np.exp(np.sum(eps[:3])) * np.asarray(sigma).ravel()   # J * sigma

    _, Lt = _umat("HOLZA", props, F0)
    d = 1e-6
    for col in range(3):
        step = np.zeros(6)
        step[col] = d
        fd = (tau_of(eps0 + step) - tau_of(eps0 - step)) / (2.0 * d)
        np.testing.assert_allclose(
            fd[:3], Lt[:3, col], rtol=2e-6,
            atol=2e-6 * max(1.0, np.abs(Lt[:3, col]).max()),
            err_msg=f"d(tau)/d(eps) column {col} (kappa_d={kappa_d})")


def test_rotating_the_fibres_and_the_motion_rotates_the_stress():
    """Frame indifference, which also exercises Rotation -> props end to end."""
    R = sim.Rotation.from_euler("zxz", [25.0, 40.0, -15.0], degrees=True)
    Rm = R.as_matrix()
    F = _F(GENERIC)

    sigma, _ = _umat("HOLZA", _props(), F)
    sigma_rot, _ = _umat("HOLZA", _props(fibres=R * FIBRES), Rm @ F @ Rm.T)

    expected = Rm @ sim.v2t_stress(sigma) @ Rm.T
    np.testing.assert_allclose(sim.v2t_stress(sigma_rot), expected, rtol=1e-9, atol=1e-11)


# ------------------------------------------------------------------- the module

def _modul(mat, eps):
    n = 1
    z6 = lambda: np.zeros((6, n), order="F")
    eye = np.eye(3).reshape(3, 3, 1).copy(order="F")
    sigma, sv, wm, Lt = sim.umat(
        "MODUL", np.asfortranarray(np.asarray(eps, float).reshape(6, 1)), z6(),
        np.array([]), np.array([]), z6(), eye,
        np.asfortranarray(np.asarray(mat.props, dtype=float).reshape(-1, 1)),
        np.zeros((mat.nstatev, n), order="F"), 0.0, 1.0,
        np.zeros((4, n), order="F"), n_threads=1)
    return sigma[:, 0], Lt[:, :, 0], sv[:, 0]


def test_modul_accepts_plasticity_and_the_return_mapping_is_consistent():
    """An anisotropic potential composed with plasticity is ACCEPTED, and the return
    mapping is correct: it gets the anisotropic tangent, so the consistency condition
    sigma_vm = sigma_Y + H p holds exactly. Only the fibre convection is approximate
    (the fibres are not reoriented by the plastic strain) -- a documented modelling
    assumption, not a run-time error."""
    # a soft bulk modulus so the deviatoric part is not swamped by the volumetric one
    block = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=0.0,
                                fibres=FIBRES, kappa=5.0)
    eps = np.array([0.30, -0.14, -0.13, 0.09, -0.04, 0.05])

    def von_mises(s):
        d = s - s[:3].mean() * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        return np.sqrt(1.5 * (d[:3] @ d[:3] + 2.0 * d[3:] @ d[3:]))

    sigma_el, _, _ = _modul(ModularMaterial(elasticity=block), eps)
    H = 1.0
    for factor in (0.5, 0.2):
        sigma_Y = factor * von_mises(sigma_el)
        mat = ModularMaterial(elasticity=block, mechanisms=[
            Plasticity(sigma_Y=sigma_Y, yield_criterion=VonMisesYield(),
                       isotropic_hardening=LinearIsotropicHardening(H=H))])
        sigma, Lt, sv = _modul(mat, eps)
        p = sv[1]                                 # statev = [T_init, p, ...]
        assert p > 0.0, "plasticity should have been activated"
        assert np.all(np.isfinite(sigma)) and np.all(np.isfinite(Lt))
        np.testing.assert_allclose(von_mises(sigma), sigma_Y + H * p, rtol=1e-8)


def test_modul_accepts_viscoelasticity():
    """Also accepted; the fibre-convection caveat and the isotropic Prony reference
    are documented, not enforced."""
    block = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=0.0,
                                fibres=FIBRES, kappa=KAPPA)
    mat = ModularMaterial(elasticity=block,
                          mechanisms=[Viscoelasticity(terms=((0.02, 0.3, 10.0, 10.0),))])
    sigma, Lt, _ = _modul(mat, np.array(GENERIC))
    assert np.all(np.isfinite(sigma)) and np.all(np.isfinite(Lt))


def test_modul_allows_damage():
    """Damage contributes NO inelastic strain, so Eel stays the total strain and the
    fibre push-forward is exact: the classic anisotropic-tissue-with-softening pairing."""
    block = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=0.0,
                                fibres=FIBRES, kappa=KAPPA)
    mat = ModularMaterial(elasticity=block,
                          mechanisms=[Damage(Y_0=1.0e8, Y_c=1.0e9)])   # never triggers
    eps = np.array(GENERIC)
    sigma, Lt, _ = _modul(mat, eps)

    # undamaged, it must still BE the standalone kernel
    sigma_ref, Lt_ref = _umat("HOLZA", _props(), _F(eps))
    tau_ref = np.exp(eps[:3].sum()) * sigma_ref
    scale = max(1.0, np.abs(tau_ref).max())
    assert np.abs(sigma - tau_ref).max() < 1e-10 * scale
    assert np.abs(Lt - Lt_ref).max() < 1e-9 * np.abs(Lt_ref).max()


def test_damage_scales_the_anisotropic_stress_uniformly():
    """(1 - D) is a SCALAR on the stress, so it applies to an anisotropic response
    exactly as to an isotropic one: the damaged stress is the undamaged one times
    (1 - D), component by component, with no preferred direction."""
    block = HolzapfelElasticity(C10=C10, k1=K1, k2=K2, kappa_d=0.0,
                                fibres=FIBRES, kappa=KAPPA)
    eps = np.array(GENERIC)
    sigma_0, Lt_0, _ = _modul(ModularMaterial(elasticity=block), eps)

    # Y = 1/2 sigma : L^-1 : sigma, the mechanism's driving force at this state
    Y = 0.5 * float(sigma_0 @ np.linalg.solve(Lt_0, sigma_0))
    Y_0, Y_c = 0.5 * Y, 4.0 * Y
    D_expected = (Y - Y_0) / (Y_c - Y_0)          # the LINEAR evolution law

    mat = ModularMaterial(elasticity=block, mechanisms=[Damage(Y_0=Y_0, Y_c=Y_c)])
    sigma, _, _ = _modul(mat, eps)

    # component by component, with no preferred direction: the whole point
    np.testing.assert_allclose(sigma, (1.0 - D_expected) * sigma_0, rtol=1e-12, atol=1e-12)


# ------------------------------------------------------------------- the props

def test_props_layout():
    block = HolzapfelElasticity(C10=1.0, k1=2.0, k2=3.0, kappa_d=0.2,
                                fibres=sim.Rotation.identity(), kappa=7.0)
    assert block.potential_params() == [1.0, 2.0, 3.0, 0.2, 1.0, 1.0, 0.0, 0.0, 7.0]


def test_directions_come_from_the_rotation():
    """a0 = R e1, so a 90 deg turn about e3 sends the fibre from e1 to e2."""
    block = HolzapfelElasticity(
        C10=1.0, k1=1.0, k2=1.0, kappa_d=0.0, kappa=1.0,
        fibres=sim.Rotation.from_euler("zxz", [0.0, 0.0, 90.0], degrees=True))
    np.testing.assert_allclose(block.directions[:, 0], [0.0, 1.0, 0.0], atol=1e-15)


def test_directions_may_be_given_as_components():
    block = HolzapfelElasticity(C10=1.0, k1=1.0, k2=1.0, kappa_d=0.0, kappa=1.0,
                                fibres=np.array([[2.0, 0.0, 0.0], [0.0, 0.0, -3.0]]))
    np.testing.assert_allclose(block.directions, [[1.0, 0.0], [0.0, 0.0], [0.0, -1.0]])


@pytest.mark.parametrize("kappa_d", [-1e-9, 0.34, 1.0])
def test_dispersion_is_validated(kappa_d):
    with pytest.raises(ValueError, match="kappa_d"):
        HolzapfelElasticity(C10=1.0, k1=1.0, k2=1.0, kappa_d=kappa_d,
                            fibres=sim.Rotation.identity(), kappa=1.0)


def test_zero_norm_direction_is_rejected():
    with pytest.raises(ValueError, match="zero norm"):
        HolzapfelElasticity(C10=1.0, k1=1.0, k2=1.0, kappa_d=0.0, kappa=1.0,
                            fibres=np.array([[0.0, 0.0, 0.0]]))
