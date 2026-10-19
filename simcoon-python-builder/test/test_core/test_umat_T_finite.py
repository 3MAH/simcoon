"""sim.umat_T under finite strain (F0, F1 given): every thermomechanical law, rate form (log_R)
or total form.

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
    with pytest.raises(ValueError, match=match):
        sim.umat_T("ELISO", _col(np.zeros(6)), _col(np.zeros(6)), _col(np.zeros(6)),
                   _cube(np.eye(3)), _col(PROPS), _col(np.zeros(1)), 0., 0.1, _col(np.zeros(4)),
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


SMA = [5.5, 0.5, 0.5, 0, 67538., 67538., 0.349, 0.349, 1e-6, 1e-6, 0., 0.0418, 0.021, 0., 10., 10.,
       300., 290., 295., 305., 0.2, 0.2, 0.2, 0.2, 300., 0.2, 2., 1e-6, 1e-3, 1., 1e8]
# name: (props, nstatev, T_init), as test_thermo_heat_contract
LAWS = {
    "ELISO": (PROPS, 1, T0),
    "ELIST": ([RHO, C_P, 1, 140000., 9000., 0.3, 0.4, 4500., 1e-6, 3e-5], 1, T0),
    "ELORT": ([RHO, C_P, 140000., 9000., 9500., 0.3, 0.32, 0.4, 4500., 4600., 3000., 1e-6, 3e-5, 2.5e-5], 1, T0),
    "EPICP": ([RHO, C_P, E, NU, ALPHA, 300., 500., 0.4], 8, T0),
    "EPKCP": ([RHO, C_P, E, NU, ALPHA, 300., 500., 0.4, 20000.], 14, T0),
    "EPJCK": ([7.85e-9, 4.75e8, 200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0], 9, T0),
    "ZENER": ([4.4, 0.656, 3000., 0.35, 1e-4, 1500., 0.35, 3000., 1200.], 8, T0),
    "PRONK": ([4.4, 0.656, 3000., 0.35, 1e-4, 1., 1500., 0.35, 3000., 1200.], 14, T0),
    "SMADI": (SMA, 17, 320.),
}


def _v2t_stress(v):
    return np.array([[v[0], v[3], v[4]], [v[3], v[1], v[5]], [v[4], v[5], v[2]]])


def _t2v_stress(m):
    return np.array([m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[0, 2], m[1, 2]])


def _spin(t):
    """Superposed rigid rotation, about two axes, 1.7 pi over the path."""
    c, s = np.cos(1.7 * np.pi * t), np.sin(1.7 * np.pi * t)
    return np.array([[1., 0, 0], [0, c, -s], [0, s, c]]) @ np.array([[c, 0, s], [0, 1., 0], [-s, 0, c]])


def _F_small(t):
    return np.eye(3) + 0.1 * (_F(t) - np.eye(3))


def _path(name, form, n=20, DT=0.5, Q=None, F=_F_small):
    """Run the law along F (times Q); Cauchy stress, Wm, Wt and the heat of every increment."""
    props, nstatev, T_init = LAWS[name]
    sigma, statev, Wm, Wt = np.zeros(6), np.zeros(nstatev), np.zeros(4), np.zeros(3)
    etot, T, r = np.zeros(6), T_init, []
    for i in range(n):
        F0, F1 = F(i / n), F((i + 1) / n)
        if Q is not None:
            F0, F1 = Q(i / n) @ F0, Q((i + 1) / n) @ F1
        H0, de, DR = _increments(form, F0, F1)
        e_start = H0 if form == "total" else _t2v_strain(DR @ _v2t_strain(etot) @ DR.T)
        out = sim.umat_T(name, _col(e_start), _col(de), _col(sigma), _cube(DR), _col(props), _col(statev),
                         0.0 if i == 0 else 1.0, 0.1, _col(Wm), _col(Wt), np.array([T]), np.array([DT]),
                         n_threads=1, start=(i == 0), F0=_cube(F0), F1=_cube(F1), kinematic_form=form)
        sigma, statev, Wm, Wt = (np.asarray(o).ravel() for o in out[:4])
        etot, T = e_start + de, T + DT
        r.append(np.asarray(out[4]).ravel()[0])
    return sigma, Wm, Wt, np.array(r)


@pytest.mark.parametrize("form", ["rate", "total"])
@pytest.mark.parametrize("name", list(LAWS))
def test_every_law_balances_its_heat_under_finite_strain(name, form):
    """The heat of the increment r = Delta gamma - theta_{n+1} Delta eta sums to Wm_d - Wt along a
    stretch-shear-rotation path with heating, per reference volume, for every law and both forms."""
    sigma, Wm, Wt, r = _path(name, form)
    assert np.abs(sigma).max() > 1.0 and np.abs(r).max() > 0.0
    np.testing.assert_allclose(r.sum(), Wm[3] - Wt[0], rtol=1e-12, atol=1e-12 * np.abs(r).sum())


@pytest.mark.parametrize("name", list(LAWS))
def test_total_form_is_objective_for_every_law(name):
    """A superposed rigid rotation Q(t) leaves ln U, T, r and the works unchanged: the Cauchy stress
    turns with Q, to round-off, anisotropic laws included (ln U lives in the reference frame)."""
    s0, Wm0, Wt0, r0 = _path(name, "total")
    s1, Wm1, Wt1, r1 = _path(name, "total", Q=_spin)
    Q = _spin(1.0)
    np.testing.assert_allclose(s1, _t2v_stress(Q @ _v2t_stress(s0) @ Q.T), rtol=0, atol=1e-12 * np.abs(s0).max())
    np.testing.assert_allclose(r1, r0, rtol=0, atol=1e-11 * np.abs(r0).max())
    np.testing.assert_allclose(np.r_[Wm1, Wt1], np.r_[Wm0, Wt0], rtol=0, atol=1e-12 * np.abs(np.r_[Wm0, Wt0]).max())


@pytest.mark.parametrize("name", ["ELISO", "EPICP", "PRONK", "SMADI"])
def test_rate_form_is_objective_to_second_order(name):
    """Under log_R the discrete increment of a superposed spin is objective only asymptotically:
    the stress and heat gaps fall as the square of the step (isotropic laws: the kernels run in the
    lab frame transported by DR, so anisotropy axes are not convected)."""
    def gap(n):
        s0, _, _, r0 = _path(name, "rate", n=n, DT=10. / n)
        s1, _, _, r1 = _path(name, "rate", n=n, DT=10. / n, Q=_spin)
        Q = _spin(1.0)
        return (np.abs(s1 - _t2v_stress(Q @ _v2t_stress(s0) @ Q.T)).max() / np.abs(s0).max(),
                abs(r1.sum() - r0.sum()) / abs(r0.sum()))
    (ds20, dr20), (ds80, dr80) = gap(20), gap(80)
    assert ds20 < 0.05 and dr20 < 0.05
    assert ds80 < ds20 / 10. and dr80 < dr20 / 10.


MECH_PROPS = {"ELISO": ([E, NU, ALPHA], 1), "EPICP": ([E, NU, ALPHA, 300., 500., 0.4], 8),
              "EPJCK": ([200000.0, 0.33, 1.0e-5, 792.0, 510.0, 0.26, 0.014, 1.0, 1.03, 293.0, 1793.0], 9)}


@pytest.mark.parametrize("form", ["rate", "total"])
@pytest.mark.parametrize("name", list(MECH_PROPS))
def test_isothermal_finite_law_is_the_mechanical_one(name, form):
    """At constant temperature the thermomechanical kernel under finite strain is the mechanical one
    of sim.umat on the same path: stress and works, bit for bit."""
    sigma_T, Wm_T, _, _ = _path(name, form, DT=0.0)
    props, nstatev = MECH_PROPS[name]
    sigma, statev, Wm, etot = np.zeros(6), np.zeros(nstatev), np.zeros(4), np.zeros(6)
    n = 20
    for i in range(n):
        F0, F1 = _F_small(i / n), _F_small((i + 1) / n)
        H0, de, DR = _increments(form, F0, F1)
        e_start = H0 if form == "total" else _t2v_strain(DR @ _v2t_strain(etot) @ DR.T)
        out = sim.umat(name, _col(e_start), _col(de), _cube(F0), _cube(F1), _col(sigma), _cube(DR), _col(props),
                       _col(statev), 0.0 if i == 0 else 1.0, 0.1, _col(Wm), temp=np.array([T0]), n_threads=1,
                       start=(i == 0), kinematic_form=form)
        sigma, statev, Wm = (np.asarray(o).ravel() for o in out[:3])
        etot = e_start + de
    np.testing.assert_allclose(sigma_T, sigma, rtol=0, atol=1e-12 * np.abs(sigma).max())
    np.testing.assert_allclose(Wm_T, Wm, rtol=0, atol=1e-12 * np.abs(Wm).max())


@pytest.mark.parametrize("name", ["EPICP", "PRONK", "SMADI"])
def test_finite_law_tends_to_the_small_strain_one(name):
    """F = I + h grad u: the finite answer (Cauchy stress, heat) differs from the small-strain call
    by O(h) relative."""
    def gap(h):
        F = lambda t: np.eye(3) + h * (_F(t) - np.eye(3))
        s_f, _, _, r_f = _path(name, "total", F=F)
        props, nstatev, T_init = LAWS[name]
        sigma, statev, Wm, Wt, etot, T, r = np.zeros(6), np.zeros(nstatev), np.zeros(4), np.zeros(3), np.zeros(6), T_init, []
        n = 20
        for i in range(n):
            de = _t2v_strain(0.5 * ((F((i + 1) / n) - F(i / n)) + (F((i + 1) / n) - F(i / n)).T))
            out = sim.umat_T(name, _col(etot), _col(de), _col(sigma), _cube(np.eye(3)), _col(props), _col(statev),
                             0.0 if i == 0 else 1.0, 0.1, _col(Wm), _col(Wt), np.array([T]), np.array([0.5]),
                             n_threads=1, start=(i == 0))
            sigma, statev, Wm, Wt = (np.asarray(o).ravel() for o in out[:4])
            etot, T = etot + de, T + 0.5
            r.append(np.asarray(out[4]).ravel()[0])
        return np.abs(s_f - sigma).max() / np.abs(sigma).max(), abs(r_f.sum() - sum(r)) / abs(sum(r))
    (ds1, dr1), (ds2, dr2) = gap(0.02), gap(0.005)
    assert ds2 < ds1 / 2.5 and dr2 < dr1 / 2.5
    assert ds2 < 0.05 and dr2 < 0.05
