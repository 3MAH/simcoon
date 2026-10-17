"""sim.umat_T under finite strain (F0, F1 given): ELISO, rate form (log_R) or total form.

Cauchy stress in and out; r, drdE, drdT per unit reference volume. For isotropic thermoelasticity
the two forms are the same law: tau = L : (ln V - alpha (T - T_init) I) under log_R, its dual
T = R^T tau R under the total form, and the same entropy eta = c_0 ln(T/T_init) + alpha tr(tau).
"""

import numpy as np
import pytest

import simcoon as sim

RHO, C_P, E, NU, ALPHA = 7.8e-9, 4.6e8, 210000., 0.3, 1.2e-5
C0 = RHO * C_P
PROPS = np.array([RHO, C_P, E, NU, ALPHA])
T0 = 293.15


def _col(a):
    return np.asfortranarray(np.asarray(a, float).reshape(-1, 1))


def _cube(a):
    return np.asfortranarray(np.asarray(a, float).reshape(3, 3, 1))


def _v2t_strain(v):
    return np.array([[v[0], v[3] / 2, v[4] / 2], [v[3] / 2, v[1], v[5] / 2], [v[4] / 2, v[5] / 2, v[2]]])


def _t2v_strain(m):
    return np.array([m[0, 0], m[1, 1], m[2, 2], 2 * m[0, 1], 2 * m[0, 2], 2 * m[1, 2]])


def _lnU(F):
    c, N = np.linalg.eigh(F.T @ F)
    return N @ np.diag(0.5 * np.log(c)) @ N.T


def _F(t):
    c, s = np.cos(0.6 * t), np.sin(0.6 * t)
    R = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    return R @ np.array([[1. + 0.3 * t, 0.4 * t, 0.], [0., 1. - 0.1 * t, 0.], [0., 0., 1. + 0.05 * t]])


def _step(form, etot, de, DR, F0, F1, sigma, statev, Wm, Wt, T, DT, start):
    out = sim.umat_T("ELISO", _col(etot), _col(de), _col(sigma), _cube(DR), _col(PROPS), _col(statev),
                     0.0 if start else 1.0, 0.1, _col(Wm), _col(Wt), np.array([T]), np.array([DT]),
                     n_threads=1, start=start, F0=_cube(F0), F1=_cube(F1), kinematic_form=form)
    return [np.asarray(o) for o in out]


def _increments(form, F0, F1):
    """Strain at the start (in the end frame) and increment of the form, and its DR."""
    if form == "total":
        H0, H1 = _t2v_strain(_lnU(F0)), _t2v_strain(_lnU(F1))
        return H0, H1 - H0, np.eye(3)
    de, _, DR, _ = sim.objective_rate("logarithmic_R", F0, F1, 0.1, return_de=True)
    return None, np.asarray(de).ravel(), np.asarray(DR)


def _run(form, n=10, DT=2.0):
    sigma, statev, Wm, Wt = np.zeros(6), np.array([T0]), np.zeros(4), np.zeros(3)
    etot, T, hist = np.zeros(6), T0, []
    for i in range(n):
        F0, F1 = _F(i / n), _F((i + 1) / n)
        H0, de, DR = _increments(form, F0, F1)
        e_start = H0 if form == "total" else _t2v_strain(DR @ _v2t_strain(etot) @ DR.T)
        s, sv, wm, wt, r, dSdE, dSdT, drdE, drdT = _step(form, e_start, de, DR, F0, F1, sigma, statev,
                                                         Wm, Wt, T, DT, i == 0)
        sigma, statev, Wm, Wt = s.ravel(), sv.ravel(), wm.ravel(), wt.ravel()
        etot, T = e_start + de, T + DT
        hist.append((sigma.copy(), r.ravel()[0], drdE.ravel(), drdT.ravel()[0], dSdT.ravel()))
    return hist


def _gap(n):
    """Largest relative gap between the rate (log_R) and the total form along the path."""
    rate, total = _run("rate", n=n, DT=20. / n), _run("total", n=n, DT=20. / n)
    dev = max(np.abs(a[0] - b[0]).max() / np.abs(b[0]).max() for a, b in zip(rate, total))
    heat = max(abs(a[1] - b[1]) / abs(b[1]) for a, b in zip(rate, total))
    dSdT = max(np.abs(a[4] - b[4]).max() / np.abs(b[4]).max() for a, b in zip(rate, total))
    return dev, heat, dSdT


def test_rate_form_converges_to_the_total_form_for_isotropic_thermoelasticity():
    """Same law, two integrations: the total form is exact at any step, the log_R increments
    integrate ln V (and ln J) only asymptotically -- first order on the stress, second order on
    the heat, which only sees tr(tau). dSdT = -3 K alpha I / J1 is the same in both."""
    dev10, heat10, dSdT10 = _gap(10)
    dev40, heat40, _ = _gap(40)
    assert dev10 < 1e-3 and heat10 < 1e-3
    assert dev40 < dev10 / 3.
    assert heat40 < heat10 / 10.
    assert dSdT10 < 1e-14


def test_dSdT_is_the_temperature_derivative_of_the_cauchy_stress():
    F0, F1 = _F(0.0), _F(1.0)
    H0, de, DR = _increments("total", F0, F1)

    def sig(DT):
        return _step("total", H0, de, DR, F0, F1, np.zeros(6), [T0], np.zeros(4), np.zeros(3), T0, DT, True)

    out = sig(5.0)
    h = 1e-3
    fd = (sig(5.0 + h)[0] - sig(5.0 - h)[0]).ravel() / (2 * h)
    np.testing.assert_allclose(out[6].ravel(), fd, rtol=0, atol=1e-7 * np.abs(fd).max())


def test_adiabatic_finite_path_follows_the_isentrope():
    """Solve r = 0 for DT at every increment (Newton on drdT): c_0 ln(T/T_init) + alpha tr(tau) = 0,
    tau = J sigma, along the whole stretch-shear-rotation path."""
    sigma, statev, Wm, Wt, T = np.zeros(6), np.array([T0]), np.zeros(4), np.zeros(3), T0
    n = 10
    for i in range(n):
        F0, F1 = _F(i / n), _F((i + 1) / n)
        H0, de, DR = _increments("total", F0, F1)
        DT = 0.0
        for _ in range(20):
            out = _step("total", H0, de, DR, F0, F1, sigma, statev, Wm, Wt, T, DT, i == 0)
            r, drdT = out[4].ravel()[0], out[8].ravel()[0]
            if abs(r) < 1e-12:
                break
            DT -= r / drdT
        sigma, statev, Wm, Wt = out[0].ravel(), out[1].ravel(), out[2].ravel(), out[3].ravel()
        T += DT
        tr_tau = np.linalg.det(F1) * sigma[:3].sum()
        np.testing.assert_allclose(T, T0 * np.exp(-ALPHA * tr_tau / C0), rtol=1e-12)


@pytest.mark.parametrize("kwargs, match", [
    (dict(umat="EPICP"), "ELISO"),
    (dict(form="jaumann"), "rate"),
    (dict(only_F0=True), "both F0 and F1"),
    (dict(no_F=True, form="total"), "needs F0 and F1"),
])
def test_finite_requests_outside_the_scope_are_refused(kwargs, match):
    F = _cube(_F(0.5))
    args = dict(F0=F, F1=F, kinematic_form=kwargs.get("form"))
    if kwargs.get("only_F0"):
        del args["F1"]
    if kwargs.get("no_F"):
        del args["F0"], args["F1"]
    props = PROPS if kwargs.get("umat", "ELISO") == "ELISO" else np.r_[PROPS, 300., 500., 0.4]
    with pytest.raises(ValueError, match=match):
        sim.umat_T(kwargs.get("umat", "ELISO"), _col(np.zeros(6)), _col(np.zeros(6)), _col(np.zeros(6)),
                   _cube(np.eye(3)), _col(props), _col(np.zeros(8)), 0., 0.1, _col(np.zeros(4)),
                   _col(np.zeros(3)), np.array([T0]), np.array([0.]), **args)


def test_rate_form_work_matches_umat():
    """Under log_R the work is the true stress power (Delta_work_conjugacy), as sim.umat: with no
    temperature change the thermal ELISO is the mechanical one, Wm included."""
    props_M = np.array([E, NU, ALPHA])
    sigma, sigma_M, statev, Wm, Wm_M, Wt = np.zeros(6), np.zeros(6), np.array([T0]), np.zeros(4), np.zeros(4), np.zeros(3)
    etot, statev_M = np.zeros(6), np.zeros(1)
    for i in range(10):
        F0, F1 = _F(i / 10), _F((i + 1) / 10)
        _, de, DR = _increments("rate", F0, F1)
        e_start = _t2v_strain(DR @ _v2t_strain(etot) @ DR.T)
        out = _step("rate", e_start, de, DR, F0, F1, sigma, statev, Wm, Wt, T0, 0.0, i == 0)
        sigma, statev, Wm, Wt = out[0].ravel(), out[1].ravel(), out[2].ravel(), out[3].ravel()
        s_M, sv_M, wm_M, _ = sim.umat("ELISO", _col(e_start), _col(de), _cube(F0), _cube(F1), _col(sigma_M), _cube(DR),
                                      _col(props_M), _col(statev_M), 0.0 if i == 0 else 1.0, 0.1, _col(Wm_M),
                                   temp=np.array([T0]), n_threads=1, start=(i == 0), corate=3)
        sigma_M, statev_M, Wm_M = s_M.ravel(), sv_M.ravel(), wm_M.ravel()
        etot = e_start + de
    np.testing.assert_allclose(sigma, sigma_M, rtol=1e-12, atol=1e-12 * np.abs(sigma_M).max())
    np.testing.assert_allclose(Wm, Wm_M, rtol=1e-12, atol=1e-15)
