"""sim.umat(ndi=..., nshr=...): condensed stress states.

With ndi < 3 the components outside the (ndi, nshr) state are stress-free and the elastic
prediction condenses them away: (1, 0) truss, sigma = E eps (not the oedometric lambda + 2 mu);
(1, 2) beam fibre, 11, 12, 13 kept. Only the kept strain components are read, so they are the
only ones a coupler passes. The plastic state is the exact uniaxial / tension-torsion J2 one.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.modular import IsotropicElasticity, LinearIsotropicHardening, ModularMaterial, Plasticity

E, NU, SY = 70000.0, 0.3, 200.0
G = E / (2.0 * (1.0 + NU))
H = 1000.0
EPICP = [E, NU, 0.0, SY, H, 1.0]   # linear hardening: R(p) = H p


def _col(a):
    return np.asfortranarray(np.asarray(a, dtype=float).reshape(-1, 1))


def _drive(name, props, nstatev, de, n_inc, **kwargs):
    etot, sigma = np.zeros(6), _col(np.zeros(6))
    statev, wm = _col(np.zeros(nstatev)), _col(np.zeros(4))
    DR = np.asfortranarray(np.eye(3)[:, :, None])
    for i in range(n_inc):
        sigma, statev, wm, Lt = sim.umat(name, _col(etot), _col(de), np.empty(0), np.empty(0), sigma, DR,
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
    np.testing.assert_allclose(np.sqrt(s * s + 3.0 * t * t), SY + H * p, rtol=1e-6)
    np.testing.assert_allclose([s, t], [E * (etot[0] - ep[0]), G * (etot[3] - ep[3])], rtol=1e-9)


def test_beam_fibre_MODUL_matches_EPICP():
    mat = ModularMaterial(elasticity=IsotropicElasticity(C1=E, C2=NU),
                          mechanisms=[Plasticity(sigma_Y=SY, isotropic_hardening=LinearIsotropicHardening(H=H))])
    de = np.array([1e-3, 0, 0, 1e-3, 0, 0])
    s_ref, _, _ = _drive("EPICP", EPICP, 8, de, 20, ndi=1, nshr=2)
    s_mod, _, _ = _drive("MODUL", mat.props, mat.nstatev, de, 20, ndi=1, nshr=2)
    np.testing.assert_allclose(s_mod, s_ref, rtol=1e-5, atol=1e-9)


def test_invalid_nshr():
    with pytest.raises(ValueError):
        _drive("EPICP", EPICP, 8, np.zeros(6), 1, ndi=1, nshr=4)


def _umat_T_once(de, **kwargs):
    props = _col([7.8e-9, 4.6e8, E, NU, 0.0])   # ELISO_T: rho, c_p, E, nu, alpha
    out = sim.umat_T("ELISO", _col(np.zeros(6)), _col(de), _col(np.zeros(6)), np.asfortranarray(np.eye(3)[:, :, None]),
                     props, _col([0.0]), 0.0, 1.0, _col(np.zeros(4)), _col(np.zeros(3)), np.full(1, 293.15),
                     np.zeros(1), n_threads=1, start=True, **kwargs)
    return out[0].ravel()


def test_umat_T_honours_nshr():
    de = np.array([1e-3, 0, 0, 5e-4, -3e-4, 0])
    np.testing.assert_allclose(_umat_T_once(de, ndi=1, nshr=2), [E * 1e-3, 0, 0, G * 5e-4, -G * 3e-4, 0],
                               rtol=1e-12, atol=1e-12)


def test_umat_T_invalid_nshr():
    with pytest.raises(ValueError):
        _umat_T_once(np.zeros(6), ndi=1, nshr=4)
