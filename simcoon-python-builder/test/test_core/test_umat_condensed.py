"""sim.umat / sim.umat_T (ndi=..., nshr=...): condensed stress states.

With ndi < 3 the components outside the (ndi, nshr) state are stress-free: the law runs 3D and
the free strain components are solved around it (Newton on the law's own tangent), so the
answer is the 3D law at the strain that makes them stress-free, whatever the law (plasticity,
closed-form viscoelasticity, SMA, MODUL, hyperelasticity under finite strain). The tangent and
the thermal derivatives leave condensed on the kept components (zero on the others): the block
an element uses as it is. The solved free components are written back into Detot (in place) and
the caller carries them with etot += Detot. At small strain a caller that does not gets the same
answer (the free increment absorbs the whole free strain); under finite strain F0 and F1 take
their free stretches from etot, which must carry them.
"""

import numpy as np
import pytest
from scipy.linalg import expm, logm

import simcoon as sim
from simcoon.modular import IsotropicElasticity, LinearIsotropicHardening, ModularMaterial, Plasticity

E, NU, SY = 70000.0, 0.3, 200.0
G = E / (2.0 * (1.0 + NU))
H = 1000.0
EPICP = [E, NU, 0.0, SY, H, 1.0]   # linear hardening: R(p) = H p
STATES = [(1, 0), (1, 2), (2, 1), (2, 3)]
I3 = np.asfortranarray(np.eye(3)[:, :, None])


def _col(a):
    return np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))


def _kept(ndi, nshr):
    K = [i for i in range(ndi)] + [3 + i for i in range(nshr)]
    return K, [i for i in range(6) if i not in K]


def _drive(name, props, nstatev, de, n_inc, **kwargs):
    etot, sigma = np.zeros(6), _col(np.zeros(6))
    statev, wm = _col(np.zeros(nstatev)), _col(np.zeros(4))
    for i in range(n_inc):
        sigma, statev, wm, Lt = sim.umat(name, _col(etot), _col(de), np.empty(0), np.empty(0), sigma, I3,
                                         _col(props), statev, float(i), 1.0, wm, temp=np.full(1, 293.15),
                                         n_threads=1, start=(i == 0), **kwargs)
        etot = etot + de
    return sigma.ravel(), statev.ravel(), etot


def test_truss_is_youngs_modulus():
    de = np.array([1e-3, 0, 0, 0, 0, 0])
    sigma, _, _ = _drive("EPICP", EPICP, 8, de, 1, ndi=1)
    np.testing.assert_allclose(sigma, [E * 1e-3, 0, 0, 0, 0, 0], rtol=1e-12, atol=1e-12)


def test_beam_fibre_keeps_11_12_13():
    de = np.array([1e-3, 0, 0, 5e-4, -3e-4, 0])
    sigma, _, _ = _drive("EPICP", EPICP, 8, de, 1, ndi=1, nshr=2)
    np.testing.assert_allclose(sigma, [E * 1e-3, 0, 0, G * 5e-4, -G * 3e-4, 0], rtol=1e-12, atol=1e-12)


def test_beam_fibre_tension_torsion_J2():
    de = np.array([1e-3, 0, 0, 1e-3, 0, 0])
    sigma, statev, etot = _drive("EPICP", EPICP, 8, de, 20, ndi=1, nshr=2)
    s, t, p = sigma[0], sigma[3], statev[1]
    ep = statev[2:8]
    np.testing.assert_allclose(sigma[[1, 2, 4, 5]], 0.0, atol=1e-12)
    assert p > 1e-3
    np.testing.assert_allclose(np.sqrt(s * s + 3.0 * t * t), SY + H * p, rtol=1e-9)
    np.testing.assert_allclose([s, t], [E * (etot[0] - ep[0]), G * (etot[3] - ep[3])], rtol=1e-9)


def test_beam_fibre_MODUL_matches_EPICP():
    mat = ModularMaterial(elasticity=IsotropicElasticity(C1=E, C2=NU),
                          mechanisms=[Plasticity(sigma_Y=SY, isotropic_hardening=LinearIsotropicHardening(H=H))])
    de = np.array([1e-3, 0, 0, 1e-3, 0, 0])
    s_ref, _, _ = _drive("EPICP", EPICP, 8, de, 20, ndi=1, nshr=2)
    s_mod, _, _ = _drive("MODUL", mat.props, mat.nstatev, de, 20, ndi=1, nshr=2)
    np.testing.assert_allclose(s_mod, s_ref, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("ndi,nshr", [(1, 4), (2, 0), (2, 2), (1, 3)])
def test_invalid_stress_state(ndi, nshr):
    """Only the states an element passes: (2, 0) would make shear 12 stress-free silently."""
    with pytest.raises(ValueError):
        _drive("EPICP", EPICP, 8, np.zeros(6), 1, ndi=ndi, nshr=nshr)


def test_planar_beam_keeps_11_12():
    de = np.array([1e-3, 0, 0, 5e-4, 0, 0])
    sigma, _, _ = _drive("EPICP", EPICP, 8, de, 1, ndi=1, nshr=1)
    np.testing.assert_allclose(sigma, [E * 1e-3, 0, 0, G * 5e-4, 0, 0], rtol=1e-12, atol=1e-12)


# --- the condensed point is the 3D law at the stress-free strain, with its exact tangent -------

ELORT = [150000., 9000., 8000., 0.3, 0.25, 0.4, 5000., 4500., 3000., 0., 0., 0.]
EPKCP = [E, NU, 0., SY, H, 1.0, 20000.]
EPCHA = [210000.0, 0.3, 0.0, 300.0, 200.0, 20.0, 30000.0, 172.0, 19500.0, 301.0]
EPHIL = [E, NU, 0., SY, H, 1.0, 0.5, 0.6, 0.7, 1.5, 1.5, 1.5]
EPJCK = [200000., 0.33, 0., 792., 510., 0.26, 0.014, 1.0, 1.03, 293., 1793.]
SMADI = [0, 70000.0, 70000.0, 0.3, 0.3, 1e-6, 1e-6, 0.0, 0.05, 0.021, 0.0, 6.0, 5.0,
         293.15, 273.15, 313.15, 333.15, 0.2, 0.2, 0.2, 0.2, 300.0, 1.4, 2.0, 1e-6, 1e-3, 1.0, 1e8]
_MOD = ModularMaterial(elasticity=IsotropicElasticity(C1=E, C2=NU),
                       mechanisms=[Plasticity(sigma_Y=SY, isotropic_hardening=LinearIsotropicHardening(H=H))])
# name: (props, nstatev, temperature, strain scale, tangent_mode); every tangent here is exact in 3D
LAWS = {
    "EPICP": (EPICP, 8, 293.15, 1.0, 2),
    "ZENER": ([3000., 0.35, 0., 1500., 0.35, 3000., 1200.], 8, 293.15, 10.0, 2),
    "PRONK": ([3000., 0.35, 0., 1., 1500., 0.35, 3000., 1200.], 14, 293.15, 10.0, 2),
    "SMADI": (SMADI, 17, 353.15, 6.0, 3),
    "MODUL": (list(_MOD.props), _MOD.nstatev, 293.15, 1.0, 2),
    "ELISO": ([E, NU, 0.0], 1, 293.15, 1.0, 2),
    "ELIST": ([3, 230000., 15000., 0.02, 0.4, 50000., 0., 0.], 1, 293.15, 1.0, 2),
    "ELORT": (ELORT, 1, 293.15, 1.0, 2),
    "ZENNK": ([3000., 0.35, 0., 1., 1500., 0.35, 3000., 1200.], 14, 293.15, 10.0, 2),
    "EPJCK": (EPJCK, 9, 293.15, 3.0, 2),
    # kinematic hardening and Hill: the algorithmic (mode 2) tangent is approximate in 3D already
    "EPKCP": (EPKCP, 14, 293.15, 1.0, 3),
    "EPCHA": (EPCHA, 33, 293.15, 1.0, 3),
    "EPHIL": (EPHIL, 33, 293.15, 1.0, 3),
}
PATH = np.array([3e-4, -1e-4, 0.0, 2e-4, -1.5e-4, 1e-4])


def _umat(name, props, nstatev, e, de, s, sv, wm, T, tm, start, ndi=3, nshr=None):
    out = sim.umat(name, _col(e), _col(de), np.empty(0), np.empty(0), _col(s), I3, _col(props), _col(sv),
                   0.0, 0.05, _col(wm), temp=np.array([T]), ndi=ndi, nshr=nshr, n_threads=1, start=start,
                   tangent_mode=tm)
    return out[0].ravel(), out[1].ravel(), out[2].ravel(), out[3][:, :, 0]


@pytest.mark.parametrize("ndi,nshr", STATES)
@pytest.mark.parametrize("name", list(LAWS))
def test_condensed_point_is_the_3D_law_with_its_exact_tangent(name, ndi, nshr):
    props, nstatev, T, scale, tm = LAWS[name]
    K, C = _kept(ndi, nshr)
    de = np.zeros(6)
    de[K] = scale * PATH[K]
    e = np.zeros(6)                                              # condensed: the caller keeps eps_K only
    cn = [np.zeros(6), np.zeros(nstatev), np.zeros(4)]
    ref = [np.zeros(6), np.zeros(6), np.zeros(nstatev), np.zeros(4)]   # 3D reference, full strain
    dC = np.zeros(len(C))
    for i in range(12):
        start_state, e_last = [x.copy() for x in cn], e.copy()
        s, sv, wm, Lt = _umat(name, props, nstatev, e, de, *cn, T, tm, i == 0, ndi, nshr)
        d3 = de.copy()
        d3[C] = dC
        for _ in range(60):
            s3, sv3, wm3, L3 = _umat(name, props, nstatev, ref[0], d3, *ref[1:], T, tm, i == 0)
            if np.abs(s3[C]).max() <= 1e-11 * np.abs(s3).max():
                break
            d3[C] -= np.linalg.solve(L3[np.ix_(C, C)], s3[C])
        dC = d3[C]
        cn = [s, sv, wm]
        e = e + de
        ref = [ref[0] + d3, s3, sv3, wm3]
    np.testing.assert_allclose(s, s3, rtol=1e-9, atol=1e-9 * np.abs(s3).max())
    np.testing.assert_array_equal(s[C], 0.0)
    np.testing.assert_allclose(sv, sv3, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(wm, wm3, rtol=1e-9, atol=1e-12)
    # the tangent: condensed on K, zero elsewhere, the derivative of the condensed stress
    np.testing.assert_array_equal(np.delete(Lt, K, axis=0), 0.0)
    np.testing.assert_array_equal(np.delete(Lt, K, axis=1), 0.0)
    h = 1e-8
    fd = np.zeros((len(K), len(K)))
    for j, k in enumerate(K):
        dp, dm = de.copy(), de.copy()
        dp[k] += h
        dm[k] -= h
        fd[:, j] = (_umat(name, props, nstatev, e_last, dp, *start_state, T, tm, False, ndi, nshr)[0][K]
                    - _umat(name, props, nstatev, e_last, dm, *start_state, T, tm, False, ndi, nshr)[0][K]) / (2 * h)
    assert np.abs(Lt[np.ix_(K, K)] - fd).max() < 1e-6 * np.abs(fd).max()


@pytest.mark.parametrize("ndi,nshr", [(1, 0), (2, 1)])
def test_free_strain_written_back_into_Detot(ndi, nshr):
    """Detot comes back with the solved free increments: etot + Detot is the total strain of the
    3D law at the stress-free point, and keeping it changes nothing to the answer."""
    K, C = _kept(ndi, nshr)
    de = np.zeros(6)
    de[K] = PATH[K]
    e = np.zeros(6)
    cn = [np.zeros(6), np.zeros(8), np.zeros(4)]
    forgetful = [np.zeros(6), np.zeros(8), np.zeros(4)]   # carries the kept strain only
    for i in range(12):
        Detot = _col(de)
        out = sim.umat("EPICP", _col(e), Detot, np.empty(0), np.empty(0), _col(cn[0]), I3, _col(EPICP),
                       _col(cn[1]), 0.0, 0.05, _col(cn[2]), temp=np.array([293.15]), ndi=ndi, nshr=nshr,
                       n_threads=1, start=(i == 0))
        np.testing.assert_array_equal(Detot.ravel()[K], de[K])
        e = e + Detot.ravel()
        cn = [out[0].ravel(), out[1].ravel(), out[2].ravel()]
        s_f, sv_f, wm_f, _ = _umat("EPICP", EPICP, 8, i * de, de, *forgetful, 293.15, 2, i == 0, ndi, nshr)
        forgetful = [s_f, sv_f, wm_f]
    np.testing.assert_allclose(cn[0], forgetful[0], rtol=1e-9, atol=1e-9 * np.abs(cn[0]).max())
    # the total strain is the 3D law's: elastic + plastic, sigma_C = 0 -> eps_C = (M sigma + EP)_C
    M = np.linalg.inv(sim.L_iso([E, NU], "Enu"))
    np.testing.assert_allclose(e[C], (M @ cn[0] + cn[1][2:8])[C], rtol=1e-8, atol=1e-12)


def test_free_strain_write_back_layouts():
    """C-ordered (6, N) arrays are written back through their strides; a read-only one is refused."""
    de = np.zeros((6, 2))
    de[0] = [1e-3, 2e-3]
    args = dict(temp=np.array([293.15]), ndi=1, n_threads=1, start=True)
    Detot = np.ascontiguousarray(de)
    sim.umat("ELISO", np.zeros((6, 2)), Detot, np.empty(0), np.empty(0), np.zeros((6, 2)),
             np.asfortranarray(np.repeat(np.eye(3)[:, :, None], 2, axis=2)), _col([E, NU, 0.0]),
             np.zeros((1, 2)), 0.0, 1.0, np.zeros((4, 2)), **args)
    np.testing.assert_allclose(Detot[[1, 2]], -NU * de[[0, 0]], rtol=1e-12)
    np.testing.assert_allclose(Detot[3:], 0.0, atol=1e-15)
    Detot = _col(de[:, 0])   # umat_T: the same write-back
    sim.umat_T("ELISO", _col(np.zeros(6)), Detot, _col(np.zeros(6)), I3, _col([RHO, C_P, E, NU, 0.0]), _col([0.0]),
               0.0, 1.0, _col(np.zeros(4)), _col(np.zeros(3)), np.array([293.15]), np.array([0.0]), ndi=1,
               n_threads=1, start=True)
    np.testing.assert_allclose(Detot.ravel()[[1, 2]], -NU * 1e-3, rtol=1e-12)
    Detot = _col(de[:, 0])
    Detot.setflags(write=False)
    with pytest.raises(ValueError, match="writable"):
        sim.umat("ELISO", _col(np.zeros(6)), Detot, np.empty(0), np.empty(0), _col(np.zeros(6)), I3,
                 _col([E, NU, 0.0]), _col([0.0]), 0.0, 1.0, _col(np.zeros(4)), **args)


def test_condensed_refusals():
    F = np.asfortranarray(np.repeat(np.eye(3)[:, :, None], 1, axis=2))
    args = ("ELISO", _col(np.zeros(6)), _col(np.zeros(6)), F, F, _col(np.zeros(6)), I3, _col([E, NU, 0.0]),
            _col([0.0]), 0.0, 1.0, _col(np.zeros(4)))
    with pytest.raises(ValueError, match="box"):
        sim.umat(*args, ndi=2, tangent_output="material", n_threads=1)
    with pytest.raises(ValueError, match="material"):
        sim.umat(*args, ndi=2, condense="exact", n_threads=1)
    with pytest.raises(ValueError, match="'rate'"):
        sim.umat(*args, ndi=2, condense="exact", tangent_output="material", corate="jaumann", n_threads=1)
    with pytest.raises(ValueError, match="approximate"):
        sim.umat(*args, ndi=2, condense="bogus", n_threads=1)


# --- finite strain: F1 follows the solved stretch through ln V ----------------------------------

R_ = expm(np.array([[0, -0.3, 0.1], [0.3, 0, -0.2], [-0.1, 0.2, 0]]))


def _F(lnv):
    return np.asfortranarray((expm(sim.v2t_strain(lnv)) @ R_)[:, :, None])


def _hyper(name, props, lnv, ndi, nshr=None):
    s, _, _, Lt = sim.umat(name, _col(np.zeros(6)), _col(np.zeros(6)), I3, _F(lnv), _col(np.zeros(6)), I3,
                           _col(props), _col([0.0]), 0.0, 1.0, _col(np.zeros(4)), ndi=ndi, nshr=nshr,
                           n_threads=1, corate=3)
    return s.ravel(), Lt[:, :, 0]


HYPER = {"NEOHC": [1000., 10000.], "NEOHI": [2000., 0.40, 0.0], "SNTVE": [100000.0, 0.3, 0.0],
         "MOORI": [0.2588, -0.0449, 10000.0], "YEOHH": [0.30, -0.010, 0.0005, 1000.0]}


@pytest.mark.parametrize("ndi,nshr", STATES)
@pytest.mark.parametrize("name,props", list(HYPER.items()))
def test_hyperelastic_condensed_through_ln_V(name, props, ndi, nshr):
    """An F-based law: the caller's F1 knows the kept stretches only; the point solves the others
    (ln V1 on the free components) and returns the Cauchy stress of that F1, with the box tangent
    d(tau)/d(ln V) condensed on the kept components."""
    K, C = _kept(ndi, nshr)
    lnv = np.zeros(6)
    lnv[K] = np.array([0.25, -0.08, 0.0, 0.12, -0.07, 0.05])[K]   # large, non-coaxial

    def solved(l):   # 3D reference: Newton on ln V1(C) so that tau_C = 0
        l3 = l.copy()
        for _ in range(60):
            s3, L3 = _hyper(name, props, l3, 3)
            tau = np.linalg.det(_F(l3)[:, :, 0]) * s3
            if np.abs(tau[C]).max() <= 1e-13 * np.abs(tau).max():
                break
            l3[C] -= np.linalg.solve(L3[np.ix_(C, C)], tau[C])
        return s3, np.linalg.det(_F(l3)[:, :, 0])

    s, Lt = _hyper(name, props, lnv, ndi, nshr)
    s3, J3 = solved(lnv)
    np.testing.assert_allclose(s, s3, rtol=1e-9, atol=1e-9 * np.abs(s3).max())
    np.testing.assert_array_equal(s[C], 0.0)
    h = 1e-6
    fd = np.zeros((len(K), len(K)))
    for j, k in enumerate(K):
        lp, lm = lnv.copy(), lnv.copy()
        lp[k] += h
        lm[k] -= h
        (sp, Jp), (sm, Jm) = solved(lp), solved(lm)
        fd[:, j] = (Jp * sp[K] - Jm * sm[K]) / (2 * h)
    assert np.abs(Lt[np.ix_(K, K)] - fd).max() < 1e-6 * np.abs(fd).max()


@pytest.mark.parametrize("ndi,nshr", STATES)
def test_MODUL_hyperelastic_condensed_matches_NEOHC(ndi, nshr):
    """MODUL with a neo-Hookean block is a log-strain box (rate form, F1 for J only): condensed,
    it is the NEOHC kernel condensed, and the free increments it writes back into Detot are the
    ln V components of the 3D stress-free point."""
    from simcoon.modular import NeoHookeanElasticity
    mod = ModularMaterial(elasticity=NeoHookeanElasticity(mu=1000.0, kappa=10000.0))
    K, C = _kept(ndi, nshr)
    lnv = np.zeros(6)
    lnv[K] = np.array([0.25, -0.08, 0.0, 0.12, -0.07, 0.05])[K]
    F1 = np.asfortranarray(expm(sim.v2t_strain(lnv))[:, :, None])   # no rotation: box = lab

    def modul(d, ndi_, nshr_):
        Detot = _col(d)
        s, _, _, Lt = sim.umat("MODUL", _col(np.zeros(6)), Detot, I3, F1 if ndi_ < 3 else
                               np.asfortranarray(expm(sim.v2t_strain(d))[:, :, None]), _col(np.zeros(6)), I3,
                               _col(mod.props), _col(np.zeros(mod.nstatev)), 0.0, 1.0, _col(np.zeros(4)),
                               ndi=ndi_, nshr=nshr_, n_threads=1, corate=3)
        return s.ravel(), Lt[:, :, 0], Detot.ravel()

    s, Lt, d_solved = modul(lnv, ndi, nshr)
    np.testing.assert_array_equal(s[C], 0.0)
    np.testing.assert_array_equal(d_solved[K], lnv[K])
    I3_ = np.asfortranarray(np.eye(3)[:, :, None])
    s_neo, Lt_neo = sim.umat("NEOHC", _col(np.zeros(6)), _col(np.zeros(6)), I3_, F1, _col(np.zeros(6)), I3_,
                             _col([1000.0, 10000.0]), _col([0.0]), 0.0, 1.0, _col(np.zeros(4)), ndi=ndi, nshr=nshr,
                             n_threads=1, corate=3)[0::3]
    np.testing.assert_allclose(s, s_neo.ravel(), rtol=1e-8, atol=1e-8 * np.abs(s).max())
    np.testing.assert_allclose(Lt[np.ix_(K, K)], Lt_neo[:, :, 0][np.ix_(K, K)], rtol=1e-7,
                               atol=1e-7 * np.abs(Lt).max())
    # the 3D law at the written-back strain: the same stress, stress-free on C
    s3, _, _ = modul(d_solved, 3, None)
    np.testing.assert_allclose(s3, s, rtol=1e-8, atol=1e-8 * np.abs(s).max())


# --- finite strain over several increments: the free stretch travels in etot -------------------

_NEO_MOD = ModularMaterial(elasticity=__import__("simcoon.modular", fromlist=["x"]).NeoHookeanElasticity(
    mu=1000.0, kappa=10000.0))
FINITE_LAWS = {"NEOHC": ([1000.0, 10000.0], 1), "MODUL": (list(_NEO_MOD.props), _NEO_MOD.nstatev)}


def _finite_path(name, ndi, nshr, n, form="rate", condense="approximate"):
    """Proportional ln V path on the kept components in n increments. The caller's F knows the kept
    stretches only; it accumulates etot += Detot (the solved strain written back)."""
    props, nstatev = FINITE_LAWS[name]
    K, _ = _kept(ndi, nshr)
    target = np.zeros(6)
    target[K] = np.array([0.25, -0.08, 0.0, 0.12, -0.07, 0.05])[K]
    e, s, sv, wm, F0 = np.zeros(6), np.zeros(6), np.zeros(nstatev), np.zeros(4), np.eye(3)
    for i in range(n):
        F1 = expm(sim.v2t_strain(target * (i + 1) / n))
        de = _col(target / n)
        out = sim.umat(name, _col(e), de, np.asfortranarray(F0[:, :, None]), np.asfortranarray(F1[:, :, None]),
                       _col(s), I3, _col(props), _col(sv), float(i), 1.0, _col(wm), ndi=ndi, nshr=nshr,
                       n_threads=1, kinematic_form=form, start=(i == 0), condense=condense,
                       tangent_output="material" if condense == "exact" else "box")
        s, sv, wm = out[0].ravel(), out[1].ravel(), out[2].ravel()
        e, F0 = e + de.ravel(), F1
    return s, e


@pytest.mark.parametrize("ndi,nshr", [(1, 0), (2, 1)])   # isotropic: the free direction stays principal
@pytest.mark.parametrize("name", list(FINITE_LAWS))
def test_finite_condensed_point_over_increments(name, ndi, nshr):
    """Elastic and coaxial: n increments land where one does. The start and end F of each increment
    take their free stretches from etot (the caller's F lacks them): a J0, J1 or F1 built from the
    caller's F alone would drift, and the written-back free strain would be counted twice."""
    _, C = _kept(ndi, nshr)
    s1, e1 = _finite_path(name, ndi, nshr, 1)
    s8, e8 = _finite_path(name, ndi, nshr, 8)
    np.testing.assert_allclose(s8, s1, rtol=1e-9, atol=1e-9 * np.abs(s1).max())
    np.testing.assert_allclose(e8[C], e1[C], rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("condense", ["approximate", "exact"])
@pytest.mark.parametrize("ndi,nshr", [(1, 0), (2, 1)])
@pytest.mark.parametrize("name", list(FINITE_LAWS))
def test_finite_condensed_point_over_increments_total_form(name, ndi, nshr, condense):
    """As above under the total form, approximate and exact (there Detot is written back on all six
    components, the law's own increment H1 - H0): n increments land where one does."""
    _, C = _kept(ndi, nshr)
    s1, e1 = _finite_path(name, ndi, nshr, 1, "total", condense)
    s8, e8 = _finite_path(name, ndi, nshr, 8, "total", condense)
    np.testing.assert_allclose(s8, s1, rtol=1e-8, atol=1e-8 * np.abs(s1).max())
    np.testing.assert_allclose(e8[C], e1[C], rtol=1e-8, atol=1e-11)


@pytest.mark.parametrize("ndi,nshr", [(1, 0), (2, 1)])
def test_exact_rate_form_over_increments(ndi, nshr):
    """condense='exact' under the rate form feeds the law the log_R increment of the rebuilt F0 -> F1:
    NEOHC (a function of F) follows the total form at every n; MODUL carries the log_R integration
    error of those increments, which vanishes with the increment."""
    s_rate, _ = _finite_path("NEOHC", ndi, nshr, 8, "rate", "exact")
    s_total, _ = _finite_path("NEOHC", ndi, nshr, 8, "total", "exact")
    np.testing.assert_allclose(s_rate, s_total, rtol=1e-9, atol=1e-9 * np.abs(s_total).max())
    ref, _ = _finite_path("MODUL", ndi, nshr, 1, "total", "exact")
    gaps = [np.abs(_finite_path("MODUL", ndi, nshr, n, "rate", "exact")[0] - ref).max() for n in (4, 8, 16)]
    assert gaps[1] < 0.3 * gaps[0] and gaps[2] < 0.3 * gaps[1]


def test_umat_T_finite_condensed_dSdT_follows_the_free_stretch():
    """Plane stress under finite strain: theta moves eps_33, hence J; dSdT is the derivative of the
    Cauchy stress, J included (finite differences on DT)."""
    props = [RHO, C_P, E, NU, ALPHA]
    lnv = np.array([0.05, -0.02, 0.0, 0.03, 0.0, 0.0])
    F1 = np.asfortranarray(expm(sim.v2t_strain(lnv))[:, :, None])

    def run(DT):
        out = sim.umat_T("ELISO", _col(np.zeros(6)), _col(lnv), _col(np.zeros(6)), I3, _col(props), _col([0.0]),
                         0.0, 1.0, _col(np.zeros(4)), _col(np.zeros(3)), np.array([293.15]), np.array([DT]),
                         ndi=2, n_threads=1, start=True, F0=I3, F1=F1)
        return out[0].ravel(), out[6].ravel()
    s, dSdT = run(20.0)
    h = 1e-3
    fd = (run(20.0 + h)[0] - run(20.0 - h)[0]) / (2 * h)
    K = [0, 1, 3]
    np.testing.assert_allclose(dSdT[K], fd[K], rtol=1e-6, atol=1e-9 * np.abs(fd).max())


# --- the total form, and condense='exact' ----------------------------------------------------
# Under kinematic_form='total' the law works on H = ln U and its dual T; the caller passes
# etot = H0, Detot = H1 - H0 on the kept components and F for its rotation. 'approximate' makes
# the law's own stress vanish on C (T_C = 0); 'exact' the second Piola-Kirchhoff S_C, i.e. the
# traction on the free material faces, on the Green-Lagrange E_C (tangent dS/dE).

def _T_of(tau, F):
    """The dual stress T of the total form: T_ab = (R^T tau R)_ab / kappa_ab in the eigenbasis of C."""
    c, N = np.linalg.eigh(F.T @ F)
    lam = np.sqrt(c)
    Rm = F @ N @ np.diag(1.0 / lam) @ N.T
    t = np.log(lam[:, None] / lam[None, :])
    small = np.abs(t) < 1e-12
    kappa = np.where(small, 1.0, t / np.sinh(np.where(small, 1.0, t)))
    return N @ ((N.T @ Rm.T @ tau @ Rm @ N) / kappa) @ N.T


def _S_of(sigma, F):
    Fi = np.linalg.inv(F)
    return sim.t2v_stress(np.linalg.det(F) * Fi @ sim.v2t_stress(sigma) @ Fi.T).ravel()


# name: (props, nstatev, strain scale[, tangent_mode])
TOTAL_LAWS = {
    "NEOHC": ([1000.0, 10000.0], 1, 1.0),
    "MODUL": (list(_NEO_MOD.props), _NEO_MOD.nstatev, 1.0),
    "ELISO": ([E, NU, 0.0], 1, 0.02),
    "EPICP": (EPICP, 8, 0.05),
    "NEOHI": (HYPER["NEOHI"], 1, 1.0),
    "SNTVE": (HYPER["SNTVE"], 1, 0.02),
    "MOORI": (HYPER["MOORI"], 1, 1.0),
    "YEOHH": (HYPER["YEOHH"], 1, 1.0),
    "ELIST": (LAWS["ELIST"][0], 1, 0.02),
    "ELORT": (ELORT, 1, 0.02),
    "EPJCK": (EPJCK, 9, 0.05),
    "EPKCP": (EPKCP, 14, 0.05, 3),
    "EPCHA": (EPCHA, 33, 0.05, 3),
    "EPHIL": (EPHIL, 33, 0.05, 3),
}
# condense='exact': the total-form laws (the rate form: test_exact_rate_form_over_increments; at small
# strain it is 'approximate', test_small_strain_exact_is_approximate)
EXACT_LAWS = list(TOTAL_LAWS)
H_PATH = np.array([0.25, -0.08, 0.0, 0.12, -0.07, 0.05])


def _total(name, H, ndi=3, nshr=None, condense="approximate", Rm=None):
    """One total-form increment from the virgin state to H (kept components; the free ones are the
    caller's guess, 0 here). Returns sigma (Cauchy), Lt, the written-back Detot, the end F."""
    props, nstatev, _, *tm = TOTAL_LAWS[name]
    Rm = R_ if Rm is None else Rm
    Detot = _col(H)
    out = sim.umat(name, _col(np.zeros(6)), Detot, I3, np.asfortranarray((Rm @ expm(sim.v2t_strain(H)))[:, :, None]),
                   _col(np.zeros(6)), I3, _col(props), _col(np.zeros(nstatev)), 0.0, 1.0, _col(np.zeros(4)),
                   ndi=ndi, nshr=nshr, n_threads=1, kinematic_form="total", start=True, condense=condense,
                   tangent_output="material" if condense == "exact" else "box", tangent_mode=(tm or [2])[0])
    H1 = Detot.ravel()
    return out[0].ravel(), out[3][:, :, 0], H1, Rm @ expm(sim.v2t_strain(H1))


@pytest.mark.parametrize("ndi,nshr", STATES)
@pytest.mark.parametrize("name", list(TOTAL_LAWS))
def test_total_condensed_is_the_3D_total_law(name, ndi, nshr):
    """T_C = 0 around the 3D total law (Newton on H_C), with the exact Schur tangent dT_K/dH_K."""
    K, C = _kept(ndi, nshr)
    scale = TOTAL_LAWS[name][2]
    H = np.zeros(6)
    H[K] = scale * H_PATH[K]
    s, Lt, H1, F1 = _total(name, H.copy(), ndi, nshr)
    h = H.copy()   # the 3D reference: Newton on H_C so that T_C = 0
    for _ in range(60):
        s3, L3, _, F = _total(name, h.copy())
        T3 = sim.t2v_stress(_T_of(np.linalg.det(F) * sim.v2t_stress(s3), F)).ravel()
        if np.abs(T3[C]).max() <= 1e-11 * np.abs(T3).max():
            break
        h[C] -= np.linalg.solve(L3[np.ix_(C, C)], T3[C])
    np.testing.assert_allclose(s, s3, rtol=1e-8, atol=1e-8 * np.abs(s3).max())
    np.testing.assert_allclose(H1[C], h[C], rtol=1e-8, atol=1e-11)
    T = sim.t2v_stress(_T_of(np.linalg.det(F1) * sim.v2t_stress(s), F1)).ravel()
    assert np.abs(T[C]).max() < 1e-9 * np.abs(T).max()
    np.testing.assert_array_equal(np.delete(Lt, K, axis=0), 0.0)
    d = 1e-5 * scale   # out of the round-off of the nearly incompressible laws (MOORI: kappa/mu ~ 4e4)
    fd = np.zeros((len(K), len(K)))
    for j, k in enumerate(K):
        Tpm = []
        for sgn in (1, -1):
            Hp = H.copy()
            Hp[k] += sgn * d
            sp, _, _, Fp = _total(name, Hp, ndi, nshr)
            Tpm.append(sim.t2v_stress(_T_of(np.linalg.det(Fp) * sim.v2t_stress(sp), Fp)).ravel()[K])
        fd[:, j] = (Tpm[0] - Tpm[1]) / (2 * d)
    assert np.abs(Lt[np.ix_(K, K)] - fd).max() < 1e-5 * np.abs(fd).max()


@pytest.mark.parametrize("ndi,nshr", [(1, 0), (1, 1), (2, 1), (1, 2), (2, 3)])
@pytest.mark.parametrize("name", EXACT_LAWS)
def test_exact_condensation_frees_the_material_faces(name, ndi, nshr):
    """condense='exact': S_C = 0, i.e. P e_3 = 0 on the free face of a shell (and of the planar
    beam), P e_2 = P e_3 = 0 for the truss; E_K stays the caller's. It coincides with
    'approximate' when the free directions are principal (truss, plane stress, isotropic)."""
    K, C = _kept(ndi, nshr)
    H = np.zeros(6)
    H[K] = TOTAL_LAWS[name][2] * H_PATH[K]
    s, _, _, F1 = _total(name, H.copy(), ndi, nshr, "exact")
    S = _S_of(s, F1)
    assert np.abs(S[C]).max() < 1e-9 * np.abs(S).max()   # the Newton stops at 1e-10 |S| + 1e-14 |L_CC|
    P = F1 @ sim.v2t_stress(S)
    faces = {(1, 0): [1, 2], (1, 1): [2], (2, 1): [2]}.get((ndi, nshr), [])
    for f in faces:
        assert np.abs(P[:, f]).max() < 1e-10 * np.abs(P).max()
    Fc = R_ @ expm(sim.v2t_strain(H))
    E_of = lambda F: sim.t2v_strain(0.5 * (F.T @ F - np.eye(3))).ravel()
    np.testing.assert_allclose(E_of(F1)[K], E_of(Fc)[K], atol=1e-14)
    if (ndi, nshr) in [(1, 0), (2, 1)]:
        s_a, _, _, _ = _total(name, H.copy(), ndi, nshr)
        np.testing.assert_allclose(s, s_a, rtol=1e-9, atol=1e-9 * np.abs(s).max())


@pytest.mark.parametrize("ndi,nshr", [(2, 1), (1, 2), (2, 3)])
def test_exact_condensation_is_the_same_under_rate_and_total(ndi, nshr):
    """NEOHC depends on F only: the exact condition (on S, at the F built from E) gives the same
    point whatever form carries the law's strain."""
    K, _ = _kept(ndi, nshr)
    H = np.zeros(6)
    H[K] = H_PATH[K]
    F1c = np.asfortranarray((R_ @ expm(sim.v2t_strain(H)))[:, :, None])
    s_t, Lt_t, _, _ = _total("NEOHC", H.copy(), ndi, nshr, "exact")
    out = sim.umat("NEOHC", _col(np.zeros(6)), _col(np.zeros(6)), I3, F1c, _col(np.zeros(6)), I3,
                   _col([1000.0, 10000.0]), _col([0.0]), 0.0, 1.0, _col(np.zeros(4)), ndi=ndi, nshr=nshr,
                   n_threads=1, corate=3, start=True, condense="exact", tangent_output="material")
    np.testing.assert_allclose(out[0].ravel(), s_t, rtol=1e-9, atol=1e-9 * np.abs(s_t).max())
    np.testing.assert_allclose(out[3][:, :, 0], Lt_t, rtol=1e-7, atol=1e-7 * np.abs(Lt_t).max())


@pytest.mark.parametrize("ndi,nshr", [(2, 1), (1, 2)])
@pytest.mark.parametrize("name", EXACT_LAWS)
def test_exact_tangent_is_dS_dE_on_the_kept_components(name, ndi, nshr):
    K, _ = _kept(ndi, nshr)
    H = np.zeros(6)
    H[K] = TOTAL_LAWS[name][2] * H_PATH[K]
    _, Lt, _, F1 = _total(name, H.copy(), ndi, nshr, "exact")
    E_of = lambda F: sim.t2v_strain(0.5 * (F.T @ F - np.eye(3))).ravel()
    E = E_of(F1)
    np.testing.assert_array_equal(np.delete(Lt, K, axis=0), 0.0)
    d = 1e-5 * TOTAL_LAWS[name][2]   # out of the round-off of the nearly incompressible laws
    fd = np.zeros((len(K), len(K)))
    for j, k in enumerate(K):
        Spm = []
        for sgn in (1, -1):
            Ep = E.copy()
            Ep[k] += sgn * d
            Hp = sim.t2v_strain(0.5 * logm(np.eye(3) + 2 * sim.v2t_strain(Ep)).real).ravel()
            sp, _, _, Fp = _total(name, Hp, ndi, nshr, "exact")
            Spm.append(_S_of(sp, Fp)[K])
        fd[:, j] = (Spm[0] - Spm[1]) / (2 * d)
    assert np.abs(Lt[np.ix_(K, K)] - fd).max() < 1e-5 * np.abs(fd).max()


def test_approximate_defect_is_second_order_in_transverse_shear():
    """(2, 3): T_33 = 0 leaves S_33 = O(gamma_t^2) (about -mu gamma_t^2, compressive); halving the
    transverse shear divides it by 4. 'exact' removes it."""
    def S33(g, condense):
        H = np.array([0.1, -0.03, 0.0, 0.05, g, 0.0])
        s, _, _, F1 = _total("NEOHC", H, 2, 3, condense, Rm=np.eye(3))
        return _S_of(s, F1)[2]
    a, b = S33(0.2, "approximate"), S33(0.1, "approximate")
    assert a < 0.0 and b < 0.0
    np.testing.assert_allclose(a / b, 4.0, rtol=0.1)
    np.testing.assert_allclose(a / (1000.0 * 0.2**2), -1.0, rtol=0.3)
    assert abs(S33(0.2, "exact")) < 1e-10 * abs(a)


def test_umat_T_total_and_exact_condensed():
    """umat_T, ELISO, plane stress: the total form condensed (T_C = 0) and condense='exact'
    (S_C = 0, coinciding here: the free direction is principal); under 'exact' dSdE, dSdT, drdE
    are dS/dE, dS/dtheta, dr/dE on the kept components (finite differences)."""
    props = [RHO, C_P, E, NU, ALPHA]
    K = [0, 1, 3]

    def run(H, DT, condense, Rm=R_):
        Detot = _col(H)
        F1 = np.asfortranarray((Rm @ expm(sim.v2t_strain(H)))[:, :, None])
        out = sim.umat_T("ELISO", _col(np.zeros(6)), Detot, _col(np.zeros(6)), I3, _col(props), _col([0.0]),
                         0.0, 1.0, _col(np.zeros(4)), _col(np.zeros(3)), np.array([293.15]), np.array([DT]),
                         ndi=2, n_threads=1, start=True, F0=I3, F1=F1, kinematic_form="total", condense=condense)
        H1 = Detot.ravel()
        return dict(s=out[0].ravel(), r=out[4].ravel()[0], L=out[5][:, :, 0], dT=out[6].ravel(),
                    rE=out[7].ravel(), rT=out[8].ravel()[0], F=Rm @ expm(sim.v2t_strain(H1)))
    H = np.array([0.04, -0.015, 0.0, 0.02, 0.0, 0.0])
    a, x = run(H, 20.0, "approximate"), run(H, 20.0, "exact")
    np.testing.assert_allclose(x["s"], a["s"], rtol=1e-9, atol=1e-9 * np.abs(a["s"]).max())
    S = _S_of(x["s"], x["F"])
    assert abs(S[2]) < 1e-9 * np.abs(S).max()
    # finite differences of S_K and r in E_K and theta
    E_of = lambda F: sim.t2v_strain(0.5 * (F.T @ F - np.eye(3))).ravel()
    E1 = E_of(x["F"])
    d, hT = 1e-7, 1e-3
    for j, k in enumerate(K):
        pm = []
        for sgn in (1, -1):
            Ep = E1.copy()
            Ep[k] += sgn * d
            Hp = sim.t2v_strain(0.5 * logm(np.eye(3) + 2 * sim.v2t_strain(Ep)).real).ravel()
            pm.append(run(Hp, 20.0, "exact"))
        fdS = (_S_of(pm[0]["s"], pm[0]["F"]) - _S_of(pm[1]["s"], pm[1]["F"]))[K] / (2 * d)
        np.testing.assert_allclose(x["L"][K, k], fdS, rtol=1e-5, atol=1e-5 * np.abs(fdS).max())
        np.testing.assert_allclose(x["rE"][k], (pm[0]["r"] - pm[1]["r"]) / (2 * d), rtol=1e-5)
    p, m = run(H, 20.0 + hT, "exact"), run(H, 20.0 - hT, "exact")
    fdT = (_S_of(p["s"], p["F"]) - _S_of(m["s"], m["F"]))[K] / (2 * hT)
    np.testing.assert_allclose(x["dT"][K], fdT, rtol=1e-6, atol=1e-9 * np.abs(fdT).max())
    np.testing.assert_allclose(x["rT"], (p["r"] - m["r"]) / (2 * hT), rtol=1e-6)
    # approximate, total: dSdT is the derivative of the Cauchy stress, J(theta) included
    p, m = run(H, 20.0 + hT, "approximate"), run(H, 20.0 - hT, "approximate")
    np.testing.assert_allclose(a["dT"][K], ((p["s"] - m["s"]) / (2 * hT))[K], rtol=1e-6,
                               atol=1e-9 * np.abs(a["dT"]).max())


# --- thermomechanical: the stress, the heat and their four derivatives -----------------------

RHO, C_P, ALPHA = 7.8e-9, 4.6e8, 1.2e-5
THERMO = {"ELISO": ([RHO, C_P, E, NU, ALPHA], 1, 1.0, 1e-8),
          "ZENER": ([4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.], 8, 10.0, 1e-8),
          "EPICP": ([RHO, C_P, E, NU, ALPHA, SY, H, 1.0], 8, 1.0, 1e-3)}   # drdE: flow direction held fixed


def _umat_T(name, props, e, de, s, sv, wm, wt, T, DT, start, ndi=3, nshr=None):
    out = sim.umat_T(name, _col(e), _col(de), _col(s), I3, _col(props), _col(sv), 0.0, 0.05, _col(wm), _col(wt),
                     np.array([T]), np.array([DT]), ndi=ndi, nshr=nshr, n_threads=1, start=start)
    return dict(s=out[0].ravel(), sv=out[1].ravel(), wm=out[2].ravel(), wt=out[3].ravel(), r=out[4].ravel()[0],
                L=out[5][:, :, 0], dT=out[6].ravel(), rE=out[7].ravel(), rT=out[8].ravel()[0])


def _umat_T_once(de, **kwargs):
    out = _umat_T("ELISO", [RHO, C_P, E, NU, 0.0], np.zeros(6), de, np.zeros(6), [0.0], np.zeros(4), np.zeros(3),
                  293.15, 0.0, True, **kwargs)
    return out["s"]


def test_umat_T_honours_nshr():
    de = np.array([1e-3, 0, 0, 5e-4, -3e-4, 0])
    np.testing.assert_allclose(_umat_T_once(de, ndi=1, nshr=2), [E * 1e-3, 0, 0, G * 5e-4, -G * 3e-4, 0],
                               rtol=1e-12, atol=1e-12)


def test_umat_T_invalid_nshr():
    with pytest.raises(ValueError):
        _umat_T_once(np.zeros(6), ndi=1, nshr=4)


@pytest.mark.parametrize("ndi,nshr", [(1, 0), (2, 1)])
@pytest.mark.parametrize("name", list(THERMO))
def test_umat_T_condensed_heat_and_derivatives(name, ndi, nshr):
    props, nstatev, scale, tol_rE = THERMO[name]
    K, C = _kept(ndi, nshr)
    de = np.zeros(6)
    de[K] = scale * PATH[K]
    T, DT = 293.15, 0.5
    e = np.zeros(6)
    cn = [np.zeros(6), np.zeros(nstatev), np.zeros(4), np.zeros(3)]
    ref = [np.zeros(6), np.zeros(6), np.zeros(nstatev), np.zeros(4), np.zeros(3)]
    dC = np.zeros(len(C))
    for i in range(12):
        start_state, T_last, e_last = [x.copy() for x in cn], T, e.copy()
        o = _umat_T(name, props, e, de, *cn, T, DT, i == 0, ndi, nshr)
        d3 = de.copy()
        d3[C] = dC
        for _ in range(60):
            o3 = _umat_T(name, props, ref[0], d3, *ref[1:], T, DT, i == 0)
            if np.abs(o3["s"][C]).max() <= 1e-11 * np.abs(o3["s"]).max():
                break
            d3[C] -= np.linalg.solve(o3["L"][np.ix_(C, C)], o3["s"][C])
        dC = d3[C]
        cn = [o["s"], o["sv"], o["wm"], o["wt"]]
        e = e + de
        ref = [ref[0] + d3, o3["s"], o3["sv"], o3["wm"], o3["wt"]]
        T += DT
    for key in ("s", "sv", "wm", "wt", "r"):
        np.testing.assert_allclose(o[key], o3[key], rtol=1e-9, atol=1e-12)
    for key in ("dT", "rE"):
        np.testing.assert_array_equal(o[key][C], 0.0)

    def f(d, dt):
        return _umat_T(name, props, e_last, d, *start_state, T_last, dt, False, ndi, nshr)
    h, hT = 1e-8, 1e-5
    fdL = np.zeros((len(K), len(K)))
    fdr = np.zeros(len(K))
    for j, k in enumerate(K):
        dp, dm = de.copy(), de.copy()
        dp[k] += h
        dm[k] -= h
        op, om = f(dp, DT), f(dm, DT)
        fdL[:, j] = (op["s"][K] - om["s"][K]) / (2 * h)
        fdr[j] = (op["r"] - om["r"]) / (2 * h)
    op, om = f(de, DT + hT), f(de, DT - hT)
    assert np.abs(o["L"][np.ix_(K, K)] - fdL).max() < 1e-6 * np.abs(fdL).max()
    assert np.abs(o["dT"][K] - (op["s"][K] - om["s"][K]) / (2 * hT)).max() < 1e-6 * np.abs(o["dT"][K]).max()
    assert np.abs(o["rE"][K] - fdr).max() < tol_rE * np.abs(fdr).max()
    fdrT = (op["r"] - om["r"]) / (2 * hT)
    assert abs(o["rT"] - fdrT) < max(1e-6, 0.1 * tol_rE) * abs(fdrT)


# --- the capability query: sim.condensation_serves answers from the laws validated above ---------

NOT_VALIDATED = ["HYPOO", "OGDEN", "HOLZA", "SMRDI", "SMAMO", "LLDM0", "MIMTN", "PYEXT", "UMEXT"]


def test_condensation_serves_exactly_the_validated_laws():
    """True exactly for the laws whose condensed point the tests above check against the 3D law
    (rate form: LAWS, HYPER, FINITE_LAWS; total: TOTAL_LAWS; condense='exact' the same); always
    True for ndi = 3, where nothing is condensed."""
    rate, total, exact = set(LAWS) | set(HYPER) | set(FINITE_LAWS), set(TOTAL_LAWS), set(EXACT_LAWS)
    for name in rate | total | set(NOT_VALIDATED):
        for ndi, nshr in STATES + [(1, 1)]:
            assert sim.condensation_serves(name, ndi, nshr) == (name in rate), name
            assert sim.condensation_serves(name, ndi, nshr, kinematic_form="rate") == (name in rate), name
            assert sim.condensation_serves(name, ndi, nshr, kinematic_form="total") == (name in total), name
        assert sim.condensation_serves(name, 3) and sim.condensation_serves(name, 3, 1, "total")
    assert all(sim.total_formulation_serves(n) for n in total)
    assert exact == total


def test_condensation_serves_refusals():
    assert not sim.condensation_serves("EPICP", 2, 0)                          # an invalid state
    assert not sim.condensation_serves("EPICP", 0)
    assert not sim.condensation_serves("EPICP", 2, kinematic_form="jaumann")   # validated under log_R only
    assert sim.condensation_serves("EPICP", 2, kinematic_form="logarithmic_R")
    with pytest.raises(ValueError):
        sim.condensation_serves("EPICP", 2, kinematic_form="totl")


@pytest.mark.parametrize("name", list(LAWS))
def test_small_strain_exact_is_approximate(name):
    """Without F there is one stress: condense='exact' gives the 'approximate' point bit for bit."""
    props, nstatev, T, scale, tm = LAWS[name]
    out = [sim.umat(name, _col(np.zeros(6)), _col(scale * PATH), np.empty(0), np.empty(0), _col(np.zeros(6)), I3,
                    _col(props), _col(np.zeros(nstatev)), 0.0, 0.05, _col(np.zeros(4)), temp=np.array([T]), ndi=2,
                    n_threads=1, start=True, tangent_mode=tm, condense=c) for c in ("approximate", "exact")]
    for a, b in zip(*out):
        np.testing.assert_array_equal(a, b)
