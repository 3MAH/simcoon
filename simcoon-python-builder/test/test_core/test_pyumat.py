"""Tests for constitutive laws written in Python (PythonUMAT / 'PYEXT').

A numpy law integrated by the C++ solver through the PYEXT callback must
reproduce the equivalent built-in kernel (ELISO, EPICP) to round-off, under
strain, stress and mixed control, including the tangent that feeds the Newton
loop, step cuts and the batch entry point ``sim.umat``.
"""

import numpy as np
import pytest

import simcoon as sim
from simcoon.solver import Block, StepMeca, solve
from solver_harness import call_pyumat_batch

E, NU = 70000.0, 0.3
ELISO_PROPS = [E, NU, 1.0e-5]
K_HARD = 1000.0
SIGMA_Y = 300.0
EPICP_LIN_PROPS = [E, NU, 1.0e-5, SIGMA_Y, K_HARD, 1.0]  # m = 1: linear isotropic hardening
EPICP_NSTATEV = 8

UNIAXIAL = ["strain"] + ["stress"] * 5


def _voigt_stress_to_strain_form(v):
    """Stress-like Voigt vector -> strain-like (engineering shear) Voigt vector."""
    return v * np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])


# ---------------------------------------------------------------------------
# Python laws used by the tests
# ---------------------------------------------------------------------------

class LinearElasticPy(sim.PythonUMAT):
    """Linear isotropic elasticity; Wm bookkeeping identical to ELISO."""

    nstatev = 0

    def __init__(self, E=E, nu=NU):
        self.L = sim.L_iso([E, nu], "Enu")
        self.calls = 0

    def integrate(self, *, Etot, DEtot, sigma, statev, Wm, ndi, **kw):
        self.calls += 1
        eps = Etot + DEtot
        L = self.L
        if ndi == 2:
            # plane stress: static condensation on the out-of-plane direction (index 2),
            # same algebra as el_pred(ndi=2) in constitutive.cpp
            F = [0, 1, 3]
            Q = np.zeros((6, 6))
            Q[np.ix_(F, F)] = L[np.ix_(F, F)] - np.outer(L[F, 2], L[2, F]) / L[2, 2]
            L = Q
        stress = L @ eps
        Wm = np.array(Wm, dtype=float)
        Wm[0] += 0.5 * (sigma + stress) @ DEtot
        Wm[1] = Wm[0]
        return stress, L, statev, Wm


class J2LinearPy(sim.PythonUMAT):
    """J2 plasticity, linear isotropic hardening R = sigma_Y + k p, radial return with the
    Simo-Hughes consistent tangent. statev = [p, EP(6)]."""

    nstatev = 7

    def __init__(self, E=E, nu=NU, sigma_y=SIGMA_Y, k=K_HARD):
        self.L = sim.L_iso([E, nu], "Enu")
        self.G = E / (2.0 * (1.0 + nu))
        self.K = E / (3.0 * (1.0 - 2.0 * nu))
        self.sigma_y = sigma_y
        self.k = k
        # projectors acting on STRAIN-like vectors (stiffness convention); the stress
        # deviator below is taken directly (a stiffness projector would halve shear stresses)
        self.Idev = sim.Tensor4.deviatoric("stiffness").mat
        self.Ivol = sim.Tensor4.volumetric("stiffness").mat

    def integrate(self, *, Etot, DEtot, sigma, statev, Wm, tangent_mode, **kw):
        p = statev[0]
        EP = statev[1:7]
        eps = Etot + DEtot
        L, G, k = self.L, self.G, self.k
        sig_tr = L @ (eps - EP)
        s_tr = sig_tr.copy()
        s_tr[:3] -= sig_tr[:3].mean()                  # deviatoric trial stress (stress Voigt)
        norm_s = np.sqrt(s_tr[:3] @ s_tr[:3] + 2.0 * (s_tr[3:] @ s_tr[3:]))
        q_tr = np.sqrt(1.5) * norm_s
        f = q_tr - (self.sigma_y + k * p)
        Lt = L.copy()
        stress = sig_tr
        if f > 0.0:
            dp = f / (3.0 * G + k)
            n = s_tr / norm_s                          # unit deviatoric direction (tensor norm)
            stress = sig_tr - 2.0 * G * dp * np.sqrt(1.5) * n
            dEP = dp * np.sqrt(1.5) * _voigt_stress_to_strain_form(n)
            EP = EP + dEP
            p = p + dp
            if tangent_mode != 0:
                beta = 1.0 - 3.0 * G * dp / q_tr
                gamma = 3.0 * G / (3.0 * G + k) - (1.0 - beta)
                n_eps = _voigt_stress_to_strain_form(n)
                Lt = 3.0 * self.K * self.Ivol + 2.0 * G * beta * self.Idev \
                    - 2.0 * G * gamma * np.outer(n, n_eps)
        Wm = np.array(Wm, dtype=float)
        dW = 0.5 * (sigma + stress) @ DEtot
        Wm[0] += dW
        statev = np.concatenate([[p], EP])
        return stress, Lt, statev, Wm, L


class StepCutAbove(LinearElasticPy):
    """Elastic law that refuses strain increments larger than a threshold."""

    def __init__(self, threshold):
        super().__init__()
        self.threshold = threshold
        self.cuts = 0

    def integrate(self, *, DEtot, **kw):
        if np.abs(DEtot).max() > self.threshold:
            self.cuts += 1
            raise sim.StepCut()
        return super().integrate(DEtot=DEtot, **kw)


def _assert_same(res_py, res_ref, rtol=1e-10, atol=1e-9, fields=("Stress", "Strain", "Wm")):
    assert res_py.status == 0 and res_ref.status == 0
    assert len(res_py) == len(res_ref)
    for f in fields:
        np.testing.assert_allclose(res_py[f], res_ref[f], rtol=rtol, atol=atol, err_msg=f)


# ---------------------------------------------------------------------------
# solve(): Python law vs built-in kernels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("control", [
    ["strain"] * 6,
    UNIAXIAL,
    ["stress"] * 6,
])
def test_linear_elastic_matches_eliso(control):
    value = [0.01, -0.002, 0.001, 0.004, 0.0, -0.003] if control[1] == "strain" \
        else [300.0, 0.0, 0.0, 50.0, 0.0, 0.0]
    step = StepMeca(control=control, value=value, ninc=25)
    res_py = solve(step, LinearElasticPy())
    res_ref = solve(step, "ELISO", ELISO_PROPS, 1)
    _assert_same(res_py, res_ref, fields=("Stress", "Strain", "Wm", "TangentMatrix"))
    assert not sim._core.has_python_umat()      # context manager cleaned up


def test_linear_elastic_cyclic_blocks_matches_eliso():
    load = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=20)
    unload = StepMeca(control=UNIAXIAL, value=[-0.01, 0, 0, 0, 0, 0], ninc=20)
    blocks = [Block(steps=[load, unload], ncycle=2)]
    res_py = solve(blocks, LinearElasticPy())
    res_ref = solve(blocks, "ELISO", ELISO_PROPS, 1)
    _assert_same(res_py, res_ref, fields=("Stress", "Strain", "Wm", "TangentMatrix"))


def test_linear_elastic_finite_strain_block_matches_eliso():
    """Under NLGEOM the Python law is a log-strain / Kirchhoff box like ELISO
    (kirchhoff_box membership): a missing entry would show up as a factor J."""
    step = StepMeca(control=UNIAXIAL, value=[0.2, 0, 0, 0, 0, 0], ninc=20)
    blocks = [Block(steps=[step], control_type="logarithmic")]
    res_py = solve(blocks, LinearElasticPy())
    res_ref = solve(blocks, "ELISO", ELISO_PROPS, 1)
    _assert_same(res_py, res_ref, rtol=1e-9, atol=1e-7,
                 fields=("Stress", "Strain", "Wm", "TangentMatrix", "F"))
    assert abs(res_py["F"][0, 0, -1] - np.exp(0.2)) < 1e-6


def test_j2_linear_hardening_matches_epicp():
    load = StepMeca(control=UNIAXIAL, value=[0.02, 0, 0, 0, 0, 0], ninc=40)
    rev = StepMeca(control=UNIAXIAL, value=[-0.02, 0, 0, 0, 0, 0], ninc=40)
    blocks = [Block(steps=[load, rev], ncycle=2)]
    res_py = solve(blocks, J2LinearPy())
    res_ref = solve(blocks, "EPICP", EPICP_LIN_PROPS, EPICP_NSTATEV)
    assert res_py.status == 0 and res_ref.status == 0
    assert res_ref["Stress"][0].max() > SIGMA_Y * 1.1          # plasticity actually happened
    np.testing.assert_allclose(res_py["Stress"], res_ref["Stress"], rtol=1e-8, atol=1e-6)
    np.testing.assert_allclose(res_py["Strain"], res_ref["Strain"], rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(res_py["Statev"][0], res_ref["Statev"][1], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(res_py["Statev"][1:7], res_ref["Statev"][2:8], rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(res_py["Wm"][0], res_ref["Wm"][0], rtol=1e-8, atol=1e-6)
    # consistent (Simo-Hughes) tangent
    np.testing.assert_allclose(res_py["TangentMatrix"], res_ref["TangentMatrix"], rtol=1e-6, atol=1e-3)


def test_j2_linear_hardening_matches_epicp_in_shear():
    """Pure shear exercises the deviatoric/engineering-shear bookkeeping of the law."""
    ctrl = ["stress", "stress", "stress", "strain", "stress", "stress"]
    load = StepMeca(control=ctrl, value=[0, 0, 0, 0.03, 0, 0], ninc=30)
    rev = StepMeca(control=ctrl, value=[0, 0, 0, -0.03, 0, 0], ninc=30)
    res_py = solve([Block(steps=[load, rev])], J2LinearPy())
    res_ref = solve([Block(steps=[load, rev])], "EPICP", EPICP_LIN_PROPS, EPICP_NSTATEV)
    assert res_ref["Stress"][3].max() > SIGMA_Y / np.sqrt(3) * 1.02      # plastic in shear
    np.testing.assert_allclose(res_py["Stress"], res_ref["Stress"], rtol=1e-8, atol=1e-6)
    np.testing.assert_allclose(res_py["Statev"][0], res_ref["Statev"][1], rtol=1e-8, atol=1e-12)


def test_props_and_nstatev_from_law_or_override():
    law = LinearElasticPy()
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=5)
    res = solve(step, law)                          # props/nstatev from the object
    assert res["Statev"].shape[0] == 0
    res = solve(step, law, props=[1.0, 2.0], nstatev=3)   # explicit override
    assert res["Statev"].shape[0] == 3
    with pytest.raises(TypeError, match="props and nstatev"):
        solve(step, "ELISO")


# ---------------------------------------------------------------------------
# step cut, errors
# ---------------------------------------------------------------------------

def test_step_cut_is_honored():
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=10, Dn_mini=1e-3)
    law = StepCutAbove(threshold=5.1e-4)            # 1e-3 per increment -> one cut each
    res = solve(step, law)
    assert res.status == 0
    assert law.cuts >= 10
    np.testing.assert_allclose(res["Stress"][0, -1], E * 0.01, rtol=1e-10)
    np.testing.assert_allclose(res["Strain"][1, -1], -NU * 0.01, atol=1e-12)


def test_step_cut_forever_aborts_without_inforce():
    """inforce=0: the solver refuses to go below Dn_mini and aborts with status 1 (it used
    to exit(0), killing the interpreter). The law is left unregistered afterwards."""
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=5, Dn_mini=1e-2)
    law = StepCutAbove(threshold=1e-9)              # cuts every genuine increment
    res = solve(step, law, inforce=0, raise_on_abort=False)
    assert res.status == 1 and len(res) == 0
    with pytest.raises(RuntimeError, match="aborted early"):
        solve(step, law, inforce=0)
    assert law.cuts > 0
    assert not sim._core.has_python_umat()


class _BadShape(LinearElasticPy):
    def integrate(self, **kw):
        return np.zeros(5), np.eye(6), kw["statev"], kw["Wm"]


class _BadCount(LinearElasticPy):
    def integrate(self, **kw):
        return np.zeros(6), np.eye(6), kw["statev"]


class _Raises(LinearElasticPy):
    def integrate(self, **kw):
        return 1 / 0


class _NonFinite(LinearElasticPy):
    def integrate(self, **kw):
        s = np.full(6, np.nan)
        return s, np.eye(6), kw["statev"], kw["Wm"]


@pytest.mark.parametrize("law, exc, match", [
    (_BadShape(), ValueError, "sigma"),
    (_BadCount(), ValueError, "4 or 5 return values"),
    (_Raises(), ZeroDivisionError, "division"),
    (_NonFinite(), ValueError, "non-finite"),
])
def test_errors_propagate_unchanged(law, exc, match):
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=5)
    with pytest.raises(exc, match=match):
        solve(step, law)
    assert not sim._core.has_python_umat()


def test_no_registered_law_is_a_clear_error():
    sim.pyumat.unregister()
    step = StepMeca(control=UNIAXIAL, value=[0.01, 0, 0, 0, 0, 0], ninc=5)
    with pytest.raises(RuntimeError, match="no UMAT callback"):
        solve(step, "PYEXT", props=[], nstatev=0)


def test_registered_restores_previous_law():
    a, b = LinearElasticPy(), LinearElasticPy()
    with sim.registered(a):
        assert sim.pyumat.current() is a
        with sim.registered(b):
            assert sim.pyumat.current() is b
        assert sim.pyumat.current() is a
        assert sim._core.has_python_umat()
    assert sim.pyumat.current() is None
    assert not sim._core.has_python_umat()


# ---------------------------------------------------------------------------
# batch entry point sim.umat("PYEXT", ...)
# ---------------------------------------------------------------------------

def _batch_inputs(n, seed=0):
    rng = np.random.default_rng(seed)
    etot = np.asfortranarray(rng.normal(scale=1e-3, size=(6, n)))
    Detot = np.asfortranarray(rng.normal(scale=1e-3, size=(6, n)))
    sigma = np.asfortranarray(sim.L_iso([E, NU], "Enu") @ etot)
    Wm = np.zeros((4, n), order="F")
    DR = np.empty((3, 3, n), order="F")
    DR[...] = np.eye(3)[:, :, None]
    return etot, Detot, sigma, Wm, DR


def test_batch_umat_matches_eliso_and_single_calls():
    n = 250                                          # > parallel cutoff (100): serial branch
    etot, Detot, sigma, Wm, DR = _batch_inputs(n)
    law = LinearElasticPy()
    with sim.registered(law):
        out = call_pyumat_batch(etot, Detot, sigma, DR, Wm)
    assert law.calls == n
    ref = sim.umat("ELISO", etot, Detot, np.array([]), np.array([]), sigma, DR,
                   np.asfortranarray(np.array(ELISO_PROPS)[:, None]),
                   np.zeros((1, n), order="F"), 0.5, 1.0, Wm)
    np.testing.assert_allclose(out[0], ref[0], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(out[3], ref[3], rtol=1e-12)
    assert out[1].shape == (0, n)
    # N single-point calls give the same columns
    with sim.registered(law):
        for pt in [0, 17, n - 1]:
            col = lambda a: a[..., [pt]].copy(order="F")
            one = call_pyumat_batch(col(etot), col(Detot), col(sigma), col(DR), col(Wm))
            np.testing.assert_allclose(one[0][:, 0], out[0][:, pt], rtol=1e-12, atol=1e-9)


def test_batch_umat_plane_stress_ndi2():
    n = 3
    etot, Detot, sigma, Wm, DR = _batch_inputs(n, seed=1)
    etot[:] = 0.0
    sigma[:] = 0.0
    with sim.registered(LinearElasticPy()):
        out = call_pyumat_batch(etot, Detot, sigma, DR, Wm, ndi=2)
    Q = out[3][:, :, 0]
    np.testing.assert_allclose(Q[0, 0], E / (1 - NU ** 2), rtol=1e-12)
    np.testing.assert_allclose(Q[0, 1], NU * E / (1 - NU ** 2), rtol=1e-12)
    np.testing.assert_allclose(Q[3, 3], E / (2 * (1 + NU)), rtol=1e-12)
    assert np.all(Q[2, :] == 0.0) and np.all(Q[:, 2] == 0.0)
    np.testing.assert_allclose(out[0][2], 0.0, atol=1e-12)     # sigma_33 = 0


def test_batch_umat_step_cut_is_a_clear_error():
    """The batch entry point cannot subdivide an increment: a StepCut must not be swallowed."""
    n = 4
    etot, Detot, sigma, Wm, DR = _batch_inputs(n)
    with sim.registered(StepCutAbove(threshold=-1.0)):
        with pytest.raises(RuntimeError, match="cannot subdivide"):
            call_pyumat_batch(etot, Detot, sigma, DR, Wm)


def test_batch_umat_stops_at_the_first_failing_point():
    """A failing point aborts the batch instead of calling the law for the remaining ones."""
    n = 20

    class FailsAt(LinearElasticPy):
        def __init__(self, at):
            super().__init__()
            self.at = at
            self.points = 0

        def integrate(self, **kw):
            self.points += 1
            if self.points > self.at:
                raise ZeroDivisionError("boom")
            return super().integrate(**kw)

    etot, Detot, sigma, Wm, DR = _batch_inputs(n)
    law = FailsAt(3)
    with sim.registered(law):
        with pytest.raises(ZeroDivisionError):
            call_pyumat_batch(etot, Detot, sigma, DR, Wm)
    assert law.points == 4                     # stopped at the failing point, not n


def test_batch_umat_without_registration_is_a_clear_error():
    sim.pyumat.unregister()
    n = 2
    etot, Detot, sigma, Wm, DR = _batch_inputs(n)
    with pytest.raises(RuntimeError, match="no UMAT callback"):
        call_pyumat_batch(etot, Detot, sigma, DR, Wm)


def test_probe_call_leaves_statev_unchanged():
    """The zero-increment tangent probe (DTime == 0, DEtot == 0) never advances statev."""

    class Counter(sim.PythonUMAT):
        nstatev = 1

        def integrate(self, *, statev, Wm, sigma, **kw):
            statev += 1.0                        # in place, on purpose
            return sigma + 5.0, np.eye(6), statev, Wm + 1.0

    law = Counter()
    probe = dict(Etot=np.zeros(6), DEtot=np.zeros(6), sigma=np.ones(6), DR=np.eye(3),
                 props=np.zeros(0), statev=np.zeros(1), T=293.15, DT=0.0, Time=0.0, DTime=0.0,
                 Wm=np.zeros(4), ndi=3, nshr=3, start=True, tangent_mode=2)
    out = law(**probe)
    np.testing.assert_array_equal(out[0], np.ones(6))    # stress, statev and Wm untouched
    assert out[2][0] == 0.0 and out[3][0] == 0.0
    real = dict(probe, DEtot=np.array([1e-3, 0, 0, 0, 0, 0]), DTime=0.1, statev=np.zeros(1))
    out = law(**real)
    assert out[2][0] == 1.0 and out[0][0] == 6.0 and out[3][0] == 1.0


def test_python_law_rejects_thermomechanical_blocks():
    """select_umat_T has no PYEXT entry: say so in Python, not from deep inside the C++ solve."""
    from simcoon.solver import StepThermomeca

    step = StepThermomeca(control="stress", value=[0.0] * 6, ninc=5,
                          thermal_control="heat_flux", Q=0.0)
    with pytest.raises(TypeError, match="thermomechanical"):
        solve(Block(steps=[step]), LinearElasticPy())
    with pytest.raises(TypeError, match="thermomechanical"):
        solve(step, LinearElasticPy())                 # bare step too
