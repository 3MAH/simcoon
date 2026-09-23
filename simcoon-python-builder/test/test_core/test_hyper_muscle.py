"""The MUSCL potential: activated skeletal muscle.

MUSCL is one potential covering a family of published laws, so it is checked
family by family: each fibre law against the reference implementation it
reproduces (ArtiSynth's, transcribed here from the Java source), the whole
potential against a numerically differentiated energy written independently of
the kernel, the degenerate settings against the isotropic potentials they must
collapse onto, and the activation against the two things that make it unusual --
it is a driven input, and it scales a passive stiffness.

The reference formulas are transcribed from the ArtiSynth SOURCE, not from its
Modeling Guide: the guide's eq. (6.32) claims the passive fibre force vanishes
below the optimal length in GenericMuscle, where the gate is in fact hard-coded
off, and eq. (6.33) omits the zero band the code applies outside [0.4, 1.6].
"""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import expm

import simcoon as sim
from simcoon.modular import (
    Damage,
    ModularMaterial,
    MooneyRivlinElasticity,
    MuscleElasticity,
    MuscleFibreLaw,
    Viscoelasticity,
    YeohElasticity,
)

# Nazari et al. (2010) Table 1, in MPa, and the fibre set ArtiSynth ships
C10, C20, KAPPA = 2.5e-3, 1.175e-3, 2.5
SIGMA_MAX, LAMBDA_OPT, LAMBDA_STAR, P1, P2 = 0.3, 1.0, 1.4, 0.05, 6.6

ALONG = sim.Rotation.identity()                                             # a0 = e1
ACROSS = sim.Rotation.from_euler("zxz", [0.0, 0.0, 90.0], degrees=True)     # a0 = e2


def _law(**kw):
    """A BLEMKER muscle with the published parameters, overridable field by field."""
    base = dict(fibre_law=MuscleFibreLaw.BLEMKER, C10=C10, C20=C20, sigma_max=SIGMA_MAX,
                lambda_opt=LAMBDA_OPT, lambda_star=LAMBDA_STAR, P1=P1, P2=P2,
                zero_below_opt=True, fibres=ALONG, kappa=KAPPA)
    if "law" in kw:
        kw["fibre_law"] = kw.pop("law")
    base.update(kw)
    return MuscleElasticity(**base)


def _umat(name, props, F1, n_points=1):
    n = n_points
    z6 = lambda: np.zeros((6, n), order="F")
    eye = np.tile(np.eye(3)[:, :, None], (1, 1, n)).copy(order="F")
    F = np.asarray(F1, dtype=float)
    F = np.tile(F[:, :, None], (1, 1, n)) if F.ndim == 2 else F
    props = np.asarray(props, dtype=float)
    props = props.reshape(-1, 1) if props.ndim == 1 else props
    stress, sv, wm, Lt = sim.umat(
        name, z6(), z6(), eye, np.asfortranarray(F), z6(), eye,
        np.asfortranarray(props), np.zeros((1, n), order="F"), 0.0, 1.0,
        np.zeros((4, n), order="F"), n_threads=1)
    return (stress, Lt) if n_points > 1 else (stress[:, 0], Lt[:, :, 0])


def _F(eps):
    return expm(sim.v2t_strain(np.asarray(eps, dtype=float)))


def _iso(F):
    """(I1_bar, I2_bar, C_bar, J) -- the kinematics every formula below shares."""
    J = np.linalg.det(F)
    C_bar = J ** (-2.0 / 3.0) * (F.T @ F)
    I1 = np.trace(C_bar)
    I2 = 0.5 * (I1 * I1 - np.trace(C_bar @ C_bar))
    return I1, I2, C_bar, J


GENERIC_STATE = [0.18, -0.06, -0.05, 0.05, -0.02, 0.03]


# ------------------------------------------------- the ArtiSynth reference kernels

def _artisynth_fd(law, lam, act):
    """f_d(lambda_bar) exactly as the ArtiSynth muscle materials compute it."""
    if law is MuscleFibreLaw.SIMPLE:
        return act * SIGMA_MAX
    if law is MuscleFibreLaw.GENERIC:
        # GenericMuscle: P1 is a STRESS, there is no optimal length, and the gate
        # below lambda_bar = 1 is hard-coded off.
        if lam <= LAMBDA_STAR:
            fp = P1 * (np.exp(P2 * (lam - 1.0)) - 1.0) / lam
        else:
            E = np.exp(P2 * (LAMBDA_STAR - 1.0))
            P3 = P1 * P2 * E
            P4 = P1 * (E - 1.0) - P3 * LAMBDA_STAR
            fp = (P3 * lam + P4) / lam
        return act * SIGMA_MAX + fp
    # BlemkerMuscle: P1 dimensionless, Hill force-length, zero band, gate ON
    r = lam / LAMBDA_OPT
    if lam <= LAMBDA_OPT:
        fpas = 0.0
    elif lam <= LAMBDA_STAR:
        fpas = P1 * (np.exp(P2 * (r - 1.0)) - 1.0)
    else:
        E = np.exp(P2 * (LAMBDA_STAR / LAMBDA_OPT - 1.0))
        P3 = P1 * P2 * E
        P4 = P1 * (E - 1.0) - P3 * LAMBDA_STAR / LAMBDA_OPT
        fpas = P3 * r + P4
    if r <= 0.6:
        fact = 9.0 * (r - 0.4) ** 2
    elif r < 1.4:
        fact = 1.0 - 4.0 * (1.0 - r) ** 2
    else:
        fact = 9.0 * (r - 1.6) ** 2
    if r < 0.4 or r > 1.6:
        fact = 0.0
    return SIGMA_MAX * (fpas + act * fact) / LAMBDA_OPT


def _artisynth_stress(F, a0, law, act):
    """sigma = 2 W4 I4 / J (a a^T - I/3), the shared tail of every muscle material."""
    a = F @ a0
    mag = np.linalg.norm(a)
    a = a / mag
    J = np.linalg.det(F)
    lam = mag * J ** (-1.0 / 3.0)
    f_d = _artisynth_fd(law, lam, act)
    return (f_d * lam / J) * (np.outer(a, a) - np.eye(3) / 3.0)


def _energy(F, law, act, s_max=1.0, kappa_d=0.0, a0=None, c10=C10, c01=0.0, c20=C20,
            c11=0.0, c02=0.0, kappa=KAPPA):
    """W(F) written out independently of the kernel, for the finite difference.

    The fibre potential is the quadrature of the PUBLISHED f_d, so states must stay
    clear of its corners, where the kernel uses a C1 continuation instead.
    """
    I1, I2, C_bar, J = _iso(F)
    s = 1.0 + (s_max - 1.0) * act
    e1, e2 = I1 - 3.0, I2 - 3.0
    W = s * (c10 * e1 + c01 * e2 + c20 * e1 * e1 + c11 * e1 * e2 + c02 * e2 * e2)
    if law is not MuscleFibreLaw.NONE:
        for i in range(a0.shape[1]):
            a = a0[:, i]
            I4 = kappa_d * I1 + (1.0 - 3.0 * kappa_d) * (a @ C_bar @ a)
            lam = np.sqrt(I4)
            pts = [p for p in (LAMBDA_OPT, LAMBDA_STAR, 0.4, 0.6, 1.4, 1.6)
                   if min(LAMBDA_OPT, lam) < p < max(LAMBDA_OPT, lam)]
            W += quad(lambda u: _artisynth_fd(law, u, act), LAMBDA_OPT, lam,
                      points=sorted(pts), limit=200, epsabs=1e-14)[0]
    return W + s * kappa * (J * np.log(J) - J + 1.0)   # the default "log" U(J)


# ------------------------------------------------------------------ degeneracies

def test_no_fibre_law_is_yeoh():
    """fibre_law = NONE with only C10 and C20 IS the second-order Yeoh potential."""
    F = _F(GENERIC_STATE)
    m = MuscleElasticity(fibre_law="none", C10=C10, C20=C20, kappa=KAPPA)
    y = YeohElasticity(C10=C10, C20=C20, C30=0.0, kappa=KAPPA)
    sig_m, Lt_m = _umat("MUSCL", m.potential_params() + [0.0], F)
    sig_y, Lt_y = _umat("YEOHH", y.potential_params() + [0.0], F)
    np.testing.assert_allclose(sig_m, sig_y, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(Lt_m, Lt_y, rtol=1e-13, atol=1e-13)


def test_no_fibre_law_is_mooney_rivlin():
    """C10 and C01 alone IS Mooney-Rivlin -- the I2 channel carries its term."""
    F = _F(GENERIC_STATE)
    m = MuscleElasticity(fibre_law="none", C10=C10, C01=7.0e-4, kappa=KAPPA)
    mr = MooneyRivlinElasticity(C10=C10, C01=7.0e-4, kappa=KAPPA)
    sig_m, Lt_m = _umat("MUSCL", m.potential_params() + [0.0], F)
    sig_r, Lt_r = _umat("MOORI", mr.potential_params() + [0.0], F)
    np.testing.assert_allclose(sig_m, sig_r, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(Lt_m, Lt_r, rtol=1e-13, atol=1e-13)


def test_inactive_muscle_with_no_passive_fibre_is_the_matrix_alone():
    """act = 0 and P1 = 0: nothing is left of the fibre term, at any stretch."""
    F = _F([0.22, -0.10, -0.09, 0.03, 0.0, 0.0])
    lawful = _law(P1=0.0, activation=0.0)
    matrix = MuscleElasticity(fibre_law="none", C10=C10, C20=C20, kappa=KAPPA)
    sig_f, Lt_f = _umat("MUSCL", lawful.potential_params() + [0.0], F)
    sig_m, Lt_m = _umat("MUSCL", matrix.potential_params() + [0.0], F)
    np.testing.assert_allclose(sig_f, sig_m, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(Lt_f, Lt_m, rtol=1e-13, atol=1e-13)


def test_passive_fibre_is_silent_below_the_optimal_length():
    """With the gate on, a shortened unactivated fibre contributes nothing."""
    F = _F([-0.12, 0.07, 0.06, 0.0, 0.0, 0.0])          # a0 = e1 is shortened
    a0 = _law().directions
    _, _, C_bar, _ = _iso(F)
    assert a0[:, 0] @ C_bar @ a0[:, 0] < 1.0, "the fibre must be SHORTENED here"
    sig_f, Lt_f = _umat("MUSCL", _law(activation=0.0).potential_params() + [0.0], F)
    matrix = MuscleElasticity(fibre_law="none", C10=C10, C20=C20, kappa=KAPPA)
    sig_m, Lt_m = _umat("MUSCL", matrix.potential_params() + [0.0], F)
    np.testing.assert_allclose(sig_f, sig_m, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(Lt_f, Lt_m, rtol=1e-13, atol=1e-13)


def test_full_dispersion_is_isotropic():
    """kappa_d = 1/3 makes A = b_bar/3: the response cannot depend on the fibres."""
    F = _F(GENERIC_STATE)
    ref = _umat("MUSCL", _law(kappa_d=1.0 / 3.0, activation=0.4).potential_params() + [0.0], F)
    for angles in ([0.0, 0.0, 90.0], [0.0, 35.0, 12.0]):
        other = sim.Rotation.from_euler("zxz", angles, degrees=True)
        sig, Lt = _umat("MUSCL", _law(kappa_d=1.0 / 3.0, activation=0.4,
                                      fibres=other).potential_params() + [0.0], F)
        np.testing.assert_allclose(sig, ref[0], rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(Lt, ref[1], rtol=1e-11, atol=1e-11)


# ------------------------------------------------------------------- the activation

def test_isometric_active_stress_is_exactly_act_sigma_max():
    """At the optimal length the active fibre stress is a SET NUMBER.

    The identity is deviatoric -- the fibre term produces
    (f_d lambda/J)(a a^T - I/3) -- so it is sigma_aa - sigma_tt that equals
    a sigma_max, free of the pressure the lateral constraint would add.
    """
    for act in (0.0, 0.35, 1.0):
        sig, _ = _umat("MUSCL", _law(activation=act).potential_params() + [0.0], np.eye(3))
        assert abs((sig[0] - sig[1]) - act * SIGMA_MAX) < 1e-14


def test_activation_scales_the_matrix_by_exactly_s_max():
    """Nazari's law: at full activation every matrix constant, and kappa, are x s_max."""
    F = _F(GENERIC_STATE)
    off, Lt_off = _umat("MUSCL", MuscleElasticity.nazari(activation=0.0).potential_params()
                        + [1.0], F)
    on, Lt_on = _umat("MUSCL", MuscleElasticity.nazari(activation=1.0).potential_params()
                      + [1.0], F)
    np.testing.assert_allclose(on, 10.0 * off, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(Lt_on, 10.0 * Lt_off, rtol=1e-13, atol=1e-13)


def test_activation_is_driven_per_material_point():
    """The contract the FE coupling rests on: one call, one activation per point.

    umat() takes props as (nprops, n_points), and umat_modular re-parses props on
    every call, so a driver supplies a per-element, per-increment activation by
    rewriting one row. Nothing else in the interface has to change.
    """
    F = _F(GENERIC_STATE)
    acts = np.array([0.0, 0.5, 1.0])
    props = np.tile(np.asarray(_law().potential_params() + [0.0])[:, None], (1, len(acts)))
    props[MuscleElasticity.activation_index, :] = acts
    sig, _ = _umat("MUSCL", props, F, n_points=len(acts))
    for k, act in enumerate(acts):
        one, _ = _umat("MUSCL", _law(activation=act).potential_params() + [0.0], F)
        np.testing.assert_allclose(sig[:, k], one, rtol=1e-13, atol=1e-15)
    assert sig[0, 0] != sig[0, 1] != sig[0, 2]


def test_modular_material_reports_where_the_activation_sits():
    mat = ModularMaterial(elasticity=_law(activation=0.42))
    assert mat.props[mat.activation_index] == pytest.approx(0.42)
    with pytest.raises(TypeError, match="no activation"):
        ModularMaterial(elasticity=YeohElasticity(C10=C10, C20=C20, C30=0.0,
                                                  kappa=KAPPA)).activation_index


# ------------------------------------------------------ against ArtiSynth's kernels

@pytest.mark.parametrize("law", [MuscleFibreLaw.SIMPLE, MuscleFibreLaw.GENERIC,
                                 MuscleFibreLaw.BLEMKER])
@pytest.mark.parametrize("lam", [0.5, 0.7, 0.9, 1.05, 1.25, 1.5, 1.8])
def test_matches_the_artisynth_muscle_materials(law, lam):
    """The fibre stress, term for term, against the transcribed Java kernels.

    The matrix constants are zeroed and the motion is isochoric, so U'(1) = 0 and
    what the kernel returns IS the fibre term ArtiSynth would return.
    """
    act = 0.6
    F = np.diag([lam, 1.0 / np.sqrt(lam), 1.0 / np.sqrt(lam)])      # J = 1 exactly
    props = _law(law=law, C10=0.0, C20=0.0, activation=act,
                 zero_below_opt=(law is MuscleFibreLaw.BLEMKER)).potential_params() + [0.0]
    sig, _ = _umat("MUSCL", props, F)
    got = np.asarray(sim.v2t_stress(sig))
    ref = _artisynth_stress(F, np.array([1.0, 0.0, 0.0]), law, act)
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-12)


# --------------------------------------------------------------------- the model

@pytest.mark.parametrize("law", [MuscleFibreLaw.NONE, MuscleFibreLaw.SIMPLE,
                                 MuscleFibreLaw.GENERIC, MuscleFibreLaw.BLEMKER])
@pytest.mark.parametrize("state", [GENERIC_STATE, [0.25, 0.03, -0.10, 0.0, 0.04, 0.0]])
def test_stress_matches_finite_difference_of_the_energy(law, state):
    """sigma = J^-1 (dW/dF) F^T, with dW/dF taken numerically. No shared machinery."""
    act = 0.55
    s_max = 10.0 if law is MuscleFibreLaw.NONE else 1.0
    kw = dict(law=law, activation=act, s_max=s_max,
              zero_below_opt=(law is MuscleFibreLaw.BLEMKER))
    props = _law(**kw).potential_params() + [0.0]
    a0 = _law(**kw).directions

    d = 1e-6
    P = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Fp, Fm = _F(state), _F(state)
            Fp[i, j] += d
            Fm[i, j] -= d
            P[i, j] = (_energy(Fp, law, act, s_max=s_max, a0=a0)
                       - _energy(Fm, law, act, s_max=s_max, a0=a0)) / (2.0 * d)
    F = _F(state)
    sigma_fd = np.asarray(sim.t2v_stress((P @ F.T) / np.linalg.det(F))).ravel()

    sigma, _ = _umat("MUSCL", props, F)
    np.testing.assert_allclose(sigma, sigma_fd, rtol=5e-6,
                               atol=5e-6 * max(1.0, np.abs(sigma).max()))


@pytest.mark.parametrize("law", [MuscleFibreLaw.SIMPLE, MuscleFibreLaw.GENERIC,
                                 MuscleFibreLaw.BLEMKER])
def test_fibre_tangent_matches_finite_difference(law):
    """The fibre TANGENT, with the fibre active and dominating the response.

    Same recipe as the HOLZA fibre-tangent check: a COAXIAL state with the fibres on
    principal axes, so ln V and the perturbation commute and the box tangent's normal
    block IS d(tau)/d(eps), needing no rate conversion. The bulk modulus is softened so
    the fibre terms are a large fraction of that block -- at the stiff default they are
    a rounding error and the check would pass whatever the tangent returned.
    """
    act = 0.6
    kappa = 0.5
    eps0 = np.array([0.26, -0.12, -0.11, 0.0, 0.0, 0.0])
    props = _law(law=law, activation=act, kappa=kappa,
                 zero_below_opt=(law is MuscleFibreLaw.BLEMKER)).potential_params() + [0.0]

    F0 = _F(eps0)
    _, _, C_bar, _ = _iso(F0)
    lam = np.sqrt(C_bar[0, 0])
    assert lam > LAMBDA_OPT, "the fibre must be STRETCHED for the passive term to answer"

    def tau_of(eps):
        sigma, _ = _umat("MUSCL", props, _F(eps))
        return np.exp(np.sum(eps[:3])) * np.asarray(sigma).ravel()   # J * sigma

    _, Lt = _umat("MUSCL", props, F0)
    d = 1e-6
    for col in range(3):
        step = np.zeros(6)
        step[col] = d
        fd = (tau_of(eps0 + step) - tau_of(eps0 - step)) / (2.0 * d)
        np.testing.assert_allclose(
            fd[:3], Lt[:3, col], rtol=5e-6,
            atol=5e-6 * max(1.0, np.abs(Lt[:3, col]).max()),
            err_msg=f"d(tau)/d(eps) column {col} ({law.name})")


@pytest.mark.parametrize("corner", [LAMBDA_OPT, 0.6 * LAMBDA_OPT, 1.4 * LAMBDA_OPT,
                                    LAMBDA_STAR, 0.4 * LAMBDA_OPT, 1.6 * LAMBDA_OPT])
def test_the_force_length_curves_are_c1_at_every_corner(corner):
    """As published, f_p and f_a are only C0 at their junctions.

    A slope that jumps there is a tangent decided by round-off in lambda_bar, so the
    kernel continues the curves across each junction. This differences the fibre stress
    on either side and requires the one-sided slopes to agree -- which they cannot do
    for the published curves at lambda_opt, 0.6 and 1.4.
    """
    act = 0.7
    props = _law(C10=0.0, C20=0.0, activation=act).potential_params() + [0.0]

    def sig11(lam):
        F = np.diag([lam, 1.0 / np.sqrt(lam), 1.0 / np.sqrt(lam)])
        return _umat("MUSCL", props, F)[0][0]

    # C1 is a statement about the limit AT the junction, so the two secants straddle it
    # with a step well inside the continuation's half-width (1e-3). Sampling further out
    # would measure how fast the slope varies, which is not what continuity means.
    h = 1e-5
    left = (sig11(corner) - sig11(corner - h)) / h
    right = (sig11(corner + h) - sig11(corner)) / h
    # Unrepaired, the passive curve's slope jumps by sigma_max P1 P2 / lambda_opt at
    # lambda_opt, i.e. 0.066 in sigma_11 here -- 6x this tolerance.
    assert abs(left - right) < 1e-2 * max(1.0, abs(left)), \
        f"slope jumps {left:.6f} -> {right:.6f} across lambda = {corner}"


def test_rotating_the_fibres_and_the_motion_rotates_the_stress():
    """Frame indifference, which also exercises Rotation -> props end to end."""
    R = sim.Rotation.from_euler("zxz", [25.0, 40.0, -15.0], degrees=True)
    Rm = R.as_matrix()
    F = _F(GENERIC_STATE)
    sig, _ = _umat("MUSCL", _law(activation=0.5).potential_params() + [0.0], F)
    sig_rot, _ = _umat("MUSCL", _law(activation=0.5, fibres=R * ALONG).potential_params()
                       + [0.0], Rm @ F @ Rm.T)
    np.testing.assert_allclose(np.asarray(sim.v2t_stress(sig_rot)),
                               Rm @ np.asarray(sim.v2t_stress(sig)) @ Rm.T,
                               rtol=1e-10, atol=1e-12)


def test_fibres_stiffen_their_own_direction():
    """A stretch along an activated fibre carries more stress than the same across it."""
    F = _F([0.20, -0.09, -0.09, 0.0, 0.0, 0.0])
    along, _ = _umat("MUSCL", _law(fibres=ALONG, activation=0.8).potential_params() + [0.0], F)
    across, _ = _umat("MUSCL", _law(fibres=ACROSS, activation=0.8).potential_params() + [0.0], F)
    assert along[0] > across[0]


# ----------------------------------------------------------- props and validation

def test_props_layout():
    """The flat layout the C++ reads back, slot by slot."""
    law = _law(activation=0.25, kappa_d=0.1)
    p = law.potential_params()
    assert p[:16] == [float(MuscleFibreLaw.BLEMKER), C10, 0.0, C20, 0.0, 0.0, 1.0,
                      0.25, SIGMA_MAX, LAMBDA_OPT, LAMBDA_STAR, P1, P2, 1.0, 0.1, 1.0]
    np.testing.assert_allclose(p[16:19], [1.0, 0.0, 0.0])       # a0 = e1
    assert p[19] == KAPPA
    assert len(p) == 17 + 3 * 1
    # NONE needs no direction at all, and the layout shortens by exactly the triplet
    assert len(MuscleElasticity.nazari().potential_params()) == 17
    assert MuscleElasticity.activation_index == 7


@pytest.mark.parametrize("kw, match", [
    (dict(activation=1.5), "activation must lie"),
    (dict(activation=-0.1), "activation must lie"),
    (dict(s_max=0.5), "s_max must be >= 1"),
    (dict(s_max=10.0), "double-counts"),
    (dict(kappa_d=0.5), "kappa_d must lie"),
    (dict(lambda_opt=-1.0), "lambda_opt must be > 0"),
    (dict(lambda_star=0.5), "lambda_star must exceed"),
    (dict(P2=0.0), "P2 must be > 0"),
    (dict(fibre_law="wobbly"), "fibre_law must be"),
])
def test_rejected_at_construction(kw, match):
    with pytest.raises(ValueError, match=match):
        _law(**kw)


def test_a_fibre_law_needs_a_direction():
    with pytest.raises(ValueError, match="needs a fibre"):
        MuscleElasticity(fibre_law="blemker", C10=C10, kappa=KAPPA)


def test_the_kernel_rejects_the_same_things_the_dataclass_does():
    """The validation is duplicated on purpose: the C++ serves callers that bypass Python."""
    props = _law().potential_params() + [0.0]
    bad = list(props)
    bad[MuscleElasticity.activation_index] = 1.5
    with pytest.raises(ValueError, match="activation"):
        _umat("MUSCL", bad, np.eye(3))
    bad = list(props)
    bad[0] = 9.0
    with pytest.raises(ValueError, match="fibre_law"):
        _umat("MUSCL", bad, np.eye(3))


def test_equality_and_hashing_survive_a_rotation_field():
    a = _law(activation=0.3)
    b = _law(activation=0.3)
    assert a == b and hash(a) == hash(b)
    assert a != _law(activation=0.4)
    assert len({a, b}) == 1


# ------------------------------------------------------------------ compositions

def test_composes_with_damage_and_viscoelasticity():
    """MUSCL is an elasticity block like any other, so the mechanisms ride on it."""
    F = _F(GENERIC_STATE)
    for mechs in ([Damage(Y_0=0.0, Y_c=1.0e3)],
                  [Viscoelasticity(terms=((0.05, 0.3, 7.0e2, 3.0e2),))]):
        mat = ModularMaterial(elasticity=_law(activation=0.5), mechanisms=mechs)
        n = 1
        z6 = lambda: np.zeros((6, n), order="F")
        eye = np.tile(np.eye(3)[:, :, None], (1, 1, n)).copy(order="F")
        stress, sv, wm, Lt = sim.umat(
            "MODUL", z6(), z6(), eye, np.asfortranarray(np.tile(F[:, :, None], (1, 1, n))),
            z6(), eye, np.asfortranarray(mat.props.reshape(-1, 1)),
            np.zeros((mat.nstatev, n), order="F"), 0.0, 1.0,
            np.zeros((4, n), order="F"), n_threads=1)
        assert np.all(np.isfinite(stress)) and np.all(np.isfinite(Lt))
