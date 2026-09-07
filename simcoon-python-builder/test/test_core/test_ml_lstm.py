"""Tests for simcoon.ml: the stress LSTM, its data/loss/training helpers and the
LSTMLaw constitutive law driven by the C++ solver through PYEXT."""

import numpy as np
import pytest

torch = pytest.importorskip("torch")

import simcoon as sim  # noqa: E402
from simcoon import ml  # noqa: E402
from solver_harness import call_pyumat_batch  # noqa: E402
from simcoon.identify import calc_cost  # noqa: E402
from simcoon.solver import StepMeca, solve  # noqa: E402

E, NU = 70000.0, 0.3
ELISO_PROPS = [E, NU, 1.0e-5]
UNIAXIAL = ["strain"] + ["stress"] * 5


# ---------------------------------------------------------------------------
# losses: torch_cost == calc_cost
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("metric", ["mse", "nmse", "nmse_per_response", "rmse", "mae", "mape", "wmape"])
def test_torch_cost_matches_calc_cost(metric):
    rng = np.random.default_rng(3)
    N, T, R = 4, 12, 3
    y_true = rng.normal(size=(N, T, R)) * 50 + 10
    y_pred = y_true + rng.normal(size=(N, T, R)) * 2
    lengths = [12, 9, 12, 5]
    mask = np.zeros((N, T))
    for i, n in enumerate(lengths):
        mask[i, :n] = 1.0
    y_true[mask == 0] = 0.0
    y_pred[mask == 0] = 0.0
    w_test = rng.uniform(0.5, 2.0, size=N)
    w_resp = rng.uniform(0.5, 2.0, size=R)
    tc = ml.torch_cost(torch.as_tensor(y_true), torch.as_tensor(y_pred),
                       w_test=w_test, w_response=w_resp, w_point=mask, metric=metric).item()
    # calc_cost keeps the padded rows with zero point weight (same numbers, same denominators)
    cc = calc_cost([y_true[i] for i in range(N)], [y_pred[i] for i in range(N)],
                   w_test=w_test, w_response=[w_resp] * N,
                   w_point=[np.repeat(mask[i][:, None], R, axis=1) for i in range(N)], metric=metric)
    assert tc == pytest.approx(cc, rel=1e-10, abs=1e-12)


def test_torch_cost_is_differentiable():
    y_true = torch.randn(2, 5, 3, dtype=torch.float64)
    y_pred = torch.randn(2, 5, 3, dtype=torch.float64, requires_grad=True)
    ml.torch_cost(y_true, y_pred, metric="nmse").backward()
    assert torch.isfinite(y_pred.grad).all()


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

def test_stress_lstm_shapes_step_and_roundtrip(tmp_path):
    torch.manual_seed(0)
    m = ml.StressLSTM(components=("xx", "yy", "xy"), features=("strain", "dstrain", "dtime"),
                      hidden_size=8, num_layers=2)
    assert (m.n_in, m.n_out, m.state_size) == (7, 3, 32)      # [h, c], 2 layers x 8
    x = torch.randn(4, 10, 7)
    y, s, psi = m(x)
    assert y.shape == (4, 10, 3) and s.shape == (4, 32) and psi is None
    # sequential stepping reproduces the sequence forward
    ss = m.zero_state(4)
    ys = []
    for t in range(10):
        yt, ss = m.step(x[:, t], ss)
        ys.append(yt)
    torch.testing.assert_close(torch.stack(ys, 1), y, rtol=1e-5, atol=1e-6)
    # scalers and save/load
    m.fit_scalers(x * 3 + 1, y * 2 - 1)
    assert torch.allclose(m.x_mean, (x * 3 + 1).reshape(-1, 7).mean(0))
    p = tmp_path / "m.pt"
    m.save(p)
    m2 = ml.StressLSTM.load(p)
    assert m2.hparams == m.hparams
    torch.testing.assert_close(m2(x)[0], m(x)[0])
    # optional free-energy head
    m3 = ml.StressLSTM(components=("xx", "yy", "xy"), hidden_size=8, psi_head=True)
    assert m3(torch.randn(2, 5, 3))[2].shape == (2, 5, 1)
    # the flat state packs and unpacks the LSTM (h, c) without loss
    h, c = m.unpack_state(s)
    assert h.shape == (2, 4, 8)
    torch.testing.assert_close(m.pack_state(h, c), s)


def test_stress_lstm_rejects_bad_config():
    with pytest.raises(ValueError):
        ml.StressLSTM(features=("dstrain",))
    with pytest.raises(ValueError):
        ml.StressLSTM(components=("xx", "ab"))


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def test_random_paths_and_generate_dataset_eliso():
    targets, ninc = ml.random_strain_paths(6, n_segments=3, n_sub=5, amplitude=0.01, seed=2)
    t2, n2 = ml.random_strain_paths(6, n_segments=3, n_sub=5, amplitude=0.01, seed=2)
    np.testing.assert_array_equal(targets, t2)
    assert targets.shape == (6, 3, 6) and ninc.shape == (6, 3) and np.abs(targets).max() <= 0.01
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc, mode="3D",
                             features=("strain", "dtime"))
    assert ds.x.shape == (6, 15, 7) and ds.y.shape == (6, 15, 6) and bool(ds.mask.all())
    L = sim.L_iso([E, NU], "Enu")
    np.testing.assert_allclose(ds.y.numpy(), ds.x[..., :6].numpy() @ L.T, rtol=1e-4, atol=1e-2)
    np.testing.assert_allclose(ds.x[..., 6].numpy(), 0.2, rtol=1e-6)   # dtime = 1/5
    tr, te = ml.split_dataset(ds, 0.5, seed=0)
    assert len(tr) + len(te) == 6 and len(te) == 3
    ys = te.as_lists()
    assert len(ys) == 3 and ys[0].shape == (15, 6)
    # variable sub-steps -> padded, masked
    targets, ninc = ml.random_strain_paths(3, n_segments=2, n_sub=6, amplitude=0.01, seed=5,
                                           variable_substeps=True)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc)
    assert ds.lengths.tolist() == ninc.sum(axis=1).tolist()


def test_generate_dataset_plane_stress_and_uniaxial_constraints():
    targets, ninc = ml.random_strain_paths(2, n_segments=2, n_sub=4, components=("xx", "yy", "xy"),
                                           amplitude=0.005, seed=0)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc, mode="plane_stress",
                             stress_components=("xx", "yy", "zz", "xy"))
    np.testing.assert_allclose(ds.y[..., 2].numpy(), 0.0, atol=1e-6)      # sigma_zz = 0
    targets, ninc = ml.random_strain_paths(2, n_segments=2, n_sub=4, components=("xx",),
                                           amplitude=0.005, seed=0)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc, mode="uniaxial")
    np.testing.assert_allclose(ds.y[..., 0].numpy(), E * ds.x[..., 0].numpy(), rtol=1e-6)


def test_load_csv_stress_lstm_format(tmp_path):
    pd = pytest.importorskip("pandas")
    rows = []
    for sid in (7, 3):
        for t in range(0, 5):
            rows.append(dict(simulation_load_id=sid, timestep=t, total_strain_xx=0.001 * t * sid,
                             total_strain_yy=0.0, total_strain_xy=0.0005 * t,
                             stress_xx=70.0 * t * sid, stress_yy=1.0, stress_xy=2.0 * t))
    p = tmp_path / "d.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    ds = ml.load_csv(p)
    assert ds.x.shape == (2, 4, 3) and ds.y.shape == (2, 4, 3)       # timestep 0 dropped
    np.testing.assert_allclose(ds.x[0, :, 0].numpy(), 0.001 * np.arange(1, 5) * 3, rtol=1e-6)


# ---------------------------------------------------------------------------
# a trained LSTM as a constitutive law
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def eliso_lstm():
    """Small LSTM trained a few seconds on ELISO paths (3D)."""
    torch.manual_seed(0)
    targets, ninc = ml.random_strain_paths(48, n_segments=4, n_sub=10, amplitude=0.01, seed=1)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc)
    tr, te = ml.split_dataset(ds, 0.25, seed=0)
    model = ml.StressLSTM(hidden_size=32, num_layers=2)
    trl, vl = ml.train(model, tr, te, epochs=200, batch_size=16, lr=3e-3, verbose=False, seed=0)
    assert trl[-1] < trl[0] * 1e-2
    return model, te


def test_evaluate_reports_identification_metrics(eliso_lstm):
    model, te = eliso_lstm
    rep = ml.evaluate(model, te, metrics=("mse", "nmse", "wmape"), per_component=True, per_sequence=True)
    assert rep["nmse"] < 0.02 and rep["wmape"] < 0.15
    assert set(rep["per_component"]) == set(ml.VOIGT)
    assert rep["per_sequence"]["nmse"].shape == (len(te),)


def test_lstm_law_in_solver_strain_and_mixed_control(eliso_lstm):
    model, _ = eliso_lstm
    law = ml.LSTMLaw(model)
    assert law.nstatev == model.state_size + 6
    step = StepMeca(control=["strain"] * 6, value=[0.008, -0.002, 0.001, 0.004, 0.0, -0.003], ninc=40)
    r = solve(step, law)
    ref = solve(step, "ELISO", ELISO_PROPS, 1)
    assert r.status == 0
    rel = np.abs(r["Stress"] - ref["Stress"]).max() / np.abs(ref["Stress"]).max()
    assert rel < 0.05
    # mixed control: the autograd tangent drives the Newton loop to convergence
    step = StepMeca(control=UNIAXIAL, value=[0.008, 0, 0, 0, 0, 0], ninc=40)
    r = solve(step, law)
    assert r.status == 0
    np.testing.assert_allclose(r["Stress"][1:, -1], 0.0, atol=1e-4)
    assert abs(r["Stress"][0, -1] - E * 0.008) / (E * 0.008) < 0.05
    assert abs(r["Strain"][1, -1] + NU * 0.008) < 0.002
    # tangent_mode = 0 -> elastic operator identified at the origin
    r0 = solve(StepMeca(control=["strain"] * 6, value=[0.001, 0, 0, 0, 0, 0], ninc=2), law,
               tangent_mode="none")
    np.testing.assert_allclose(r0["TangentMatrix"][:, :, -1], law.elastic_L, rtol=1e-10)


def test_lstm_law_state_lives_in_statev(eliso_lstm):
    """Newton re-calls of one increment must all start from the same (h, c): the
    recurrent state is rolled back by the solver through statev."""
    model, _ = eliso_lstm

    class Spy(ml.LSTMLaw):
        def __init__(self, m):
            super().__init__(m)
            self.seen = {}

        def integrate(self, *, Time, DTime, statev, **kw):
            key = (round(Time, 12), round(DTime, 12))
            self.seen.setdefault(key, []).append(np.array(statev))
            return super().integrate(Time=Time, DTime=DTime, statev=statev, **kw)

    law = Spy(model)
    step = StepMeca(control=UNIAXIAL, value=[0.008, 0, 0, 0, 0, 0], ninc=8)
    r = solve(step, law)
    assert r.status == 0
    multi = [v for v in law.seen.values() if len(v) > 1]
    assert multi, "expected several Newton evaluations per increment"
    for calls in multi:
        for sv in calls[1:]:
            np.testing.assert_array_equal(sv, calls[0])
    # the zero-increment probe (Time = 0, DTime = 0) did not advance the state: the
    # first real increment starts from h = c = 0 ...
    real = sorted(k for k in law.seen if k[1] > 0.0)
    np.testing.assert_array_equal(law.seen[real[0]][0], 0.0)
    # ... and the state does evolve between increments
    assert np.abs(law.seen[real[-1]][0] - law.seen[real[0]][0]).max() > 0.0


def test_lstm_law_batch_equals_single_and_autograd_matches_fd(eliso_lstm):
    model, _ = eliso_lstm
    law = ml.LSTMLaw(model)
    rng = np.random.default_rng(0)
    N = 5
    s6 = rng.normal(scale=3e-3, size=(6, N))
    h0 = np.zeros((law.state_size, N))
    st, Lt, h1, used = law.step_batch(s6, h0)
    assert st.shape == (6, N) and Lt.shape == (6, 6, N) and h1.shape == (law.state_size, N)
    for i in range(N):
        one = law.step_batch(s6[:, [i]], h0[:, [i]])
        np.testing.assert_allclose(one[0][:, 0], st[:, i], rtol=1e-10, atol=1e-8)
        np.testing.assert_allclose(one[1][:, :, 0], Lt[:, :, i], rtol=1e-8, atol=1e-6)
        np.testing.assert_allclose(one[2][:, 0], h1[:, i], rtol=1e-10, atol=1e-12)
    # central finite differences on the first point
    d = 1e-6
    Jfd = np.zeros((6, 6))
    for j in range(6):
        ep, em = s6[:, [0]].copy(), s6[:, [0]].copy()
        ep[j] += d
        em[j] -= d
        Jfd[:, j] = (law.step_batch(ep, h0[:, [0]])[0][:, 0]
                     - law.step_batch(em, h0[:, [0]])[0][:, 0]) / (2 * d)
    np.testing.assert_allclose(Lt[:, :, 0], Jfd, rtol=1e-5, atol=1e-5 * np.abs(Jfd).max())


def test_lstm_law_through_sim_umat_batch(eliso_lstm):
    model, _ = eliso_lstm
    law = ml.LSTMLaw(model)
    n = 4
    rng = np.random.default_rng(1)
    Detot = np.asfortranarray(rng.normal(scale=2e-3, size=(6, n)))
    zeros6 = np.zeros((6, n), order="F")
    DR = np.empty((3, 3, n), order="F")
    DR[...] = np.eye(3)[:, :, None]
    with sim.registered(law):
        sig, sv, wm, Lt = call_pyumat_batch(zeros6, Detot, zeros6, DR, np.zeros((4, n), order="F"),
                                            nstatev=law.nstatev, time=0.0)
    st, Lt_b, h1, used = law.step_batch(Detot, np.zeros((law.state_size, n)),
                                        dtime=1.0, temperature=0.0)
    np.testing.assert_allclose(sig, st, rtol=1e-10, atol=1e-8)
    np.testing.assert_allclose(Lt, Lt_b, rtol=1e-8, atol=1e-6)
    np.testing.assert_allclose(sv[:law.state_size], h1, rtol=1e-10, atol=1e-12)
    assert sv.shape == (law.nstatev, n)


# ---------------------------------------------------------------------------
# ndi: plane stress / uniaxial condensation on a linear oracle
# ---------------------------------------------------------------------------

class _LinearOracle(ml.StressLSTM):
    """StressLSTM interface, exact linear elastic response (state ignored)."""

    def __init__(self, L):
        super().__init__(hidden_size=4, num_layers=1)
        self.register_buffer("Lmat", torch.as_tensor(L, dtype=torch.float64))

    def forward(self, x, s=None):
        y = x[..., :6] @ self.Lmat.to(x.dtype).T
        if s is None:
            s = self.zero_state(x.shape[0], dtype=x.dtype, device=x.device)
        return y, s, None


def test_condensation_plane_stress_and_uniaxial_on_linear_oracle():
    L = sim.L_iso([E, NU], "Enu")
    law = ml.LSTMLaw(_LinearOracle(L))
    np.testing.assert_allclose(law.elastic_L, L, rtol=1e-10)
    eps = np.zeros((6, 1))
    eps[0, 0] = 0.01
    eps[1, 0] = 0.002
    eps[3, 0] = 0.004
    # ndi = 2: sigma_zz = 0, condensed tangent = plane-stress stiffness
    st, Lt, h, used = law.step_batch(eps, np.zeros((law.state_size, 1)), ndi=2)
    assert st[2, 0] == 0.0
    np.testing.assert_allclose(Lt[0, 0, 0], E / (1 - NU ** 2), rtol=1e-10)
    np.testing.assert_allclose(Lt[0, 1, 0], NU * E / (1 - NU ** 2), rtol=1e-10)
    np.testing.assert_allclose(Lt[3, 3, 0], E / (2 * (1 + NU)), rtol=1e-10)
    assert np.all(Lt[2, :, 0][[0, 1, 3, 4, 5]] == 0.0) and np.all(Lt[:, 2, 0][[0, 1, 3, 4, 5]] == 0.0)
    np.testing.assert_allclose(used[2, 0], -NU / (1 - NU) * (0.01 + 0.002), rtol=1e-10)
    np.testing.assert_allclose(st[0, 0], (E / (1 - NU ** 2)) * (0.01 + NU * 0.002), rtol=1e-10)
    # ndi = 1: sigma_yy = sigma_zz = 0, axial stiffness E
    st, Lt, h, used = law.step_batch(eps[:, [0]] * np.array([[1], [0], [0], [0], [0], [0]]),
                                    np.zeros((law.state_size, 1)), ndi=1)
    np.testing.assert_allclose(Lt[0, 0, 0], E, rtol=1e-10)
    np.testing.assert_allclose(st[0, 0], E * 0.01, rtol=1e-10)
    np.testing.assert_allclose(used[1, 0], -NU * 0.01, rtol=1e-10)
    # through sim.umat with ndi=2 (the solver-side convention)
    zeros6 = np.zeros((6, 1), order="F")
    DR = np.eye(3).reshape(3, 3, 1).copy(order="F")
    with sim.registered(law):
        sig, sv, wm, Lt6 = call_pyumat_batch(zeros6, np.asfortranarray(eps), zeros6, DR,
                                             np.zeros((4, 1), order="F"), nstatev=law.nstatev,
                                             time=0.0, ndi=2)
    np.testing.assert_allclose(Lt6[0, 0, 0], E / (1 - NU ** 2), rtol=1e-10)
    assert sig[2, 0] == 0.0
    # tangent_mode = 0 under plane stress returns the condensed ELASTIC operator
    with sim.registered(law):
        _, _, _, Lt0 = call_pyumat_batch(zeros6, np.asfortranarray(eps), zeros6, DR,
                                         np.zeros((4, 1), order="F"), nstatev=law.nstatev,
                                         time=0.0, ndi=2, tangent_mode=0)
    np.testing.assert_allclose(Lt0[:, :, 0], law._elastic_by_ndi[2], rtol=1e-12)
    assert Lt0[2, 2, 0] == 0.0 and Lt0[0, 0, 0] == pytest.approx(E / (1 - NU ** 2))


# ---------------------------------------------------------------------------
# scikit-learn estimator
# ---------------------------------------------------------------------------

def test_sklearn_regressor_fit_predict_score_clone():
    sklearn = pytest.importorskip("sklearn")
    from sklearn.base import clone
    from sklearn.model_selection import cross_val_score

    targets, ninc = ml.random_strain_paths(24, n_segments=3, n_sub=6, amplitude=0.01,
                                           components=("xx", "yy", "xy"), seed=4)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc, mode="plane_strain")
    X, Y = ds.x.numpy(), ds.y.numpy()
    reg = ml.StressLSTMRegressor(components=("xx", "yy", "xy"), hidden_size=16, epochs=80,
                                 batch_size=8, lr=3e-3, random_state=0)
    params = reg.get_params()
    assert params["hidden_size"] == 16 and clone(reg).get_params() == params
    reg.fit(X, Y)
    assert reg.predict(X).shape == Y.shape
    assert reg.score(X, Y) > 0.9                       # R2
    assert reg.score(X, Y, metric="nmse") < 0.0        # negated cost
    rep = reg.evaluate(X, Y, metrics=("nmse",))
    assert rep["nmse"] < 0.1
    law = reg.law()
    assert isinstance(law, ml.LSTMLaw) and law.components == ("xx", "yy", "xy")
    scores = cross_val_score(clone(reg).set_params(epochs=30), X, Y, cv=2)
    assert scores.shape == (2,)


# ---------------------------------------------------------------------------
# time-discretisation robustness: random resampling (training) and committed state (inference)
# ---------------------------------------------------------------------------

def test_resample_time_is_exact_on_linear_paths_and_rebuilds_features():
    from simcoon.ml.data import resample_time
    targets, ninc = ml.random_strain_paths(3, n_segments=2, n_sub=8, amplitude=0.01, seed=3)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc,
                             features=("strain", "dstrain", "dtime"))
    for n in (5, 16, 41):
        x, y, mask, _ = resample_time(ds.x, ds.y, ds.mask, n, ds.meta["features"], 6)
        assert x.shape == (3, n, 13) and y.shape == (3, n, 6) and bool(mask.all())
        # linear elasticity is exactly preserved by linear interpolation
        L = sim.L_iso([E, NU], "Enu")
        np.testing.assert_allclose(y.numpy(), x[..., :6].numpy() @ L.T, rtol=1e-4, atol=1e-2)
        # end points kept, dstrain rebuilt as the backward difference, dtime rescaled
        np.testing.assert_allclose(x[:, -1, :6].numpy(), ds.x[:, -1, :6].numpy(), rtol=1e-6)
        np.testing.assert_allclose(x[:, 1:, 6:12].numpy(), (x[:, 1:, :6] - x[:, :-1, :6]).numpy(), atol=1e-7)
        # the resampled sequence describes the same path over the same duration
        np.testing.assert_allclose(x[..., 12].sum(1).numpy(), ds.x[..., 12].sum(1).numpy(), rtol=1e-6)


def test_train_records_median_increment_and_sets_commit_tol(eliso_lstm):
    targets, ninc = ml.random_strain_paths(16, n_segments=2, n_sub=10, amplitude=0.01, seed=9)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc)
    model = ml.StressLSTM(hidden_size=8, num_layers=1)
    tr, _ = ml.train(model, ds, epochs=3, batch_size=8, verbose=False, seed=0, clip_grad_norm=1e-3)
    assert len(tr) == 3 and all(np.isfinite(tr))
    d = (ds.x[:, 1:, :6] - ds.x[:, :-1, :6]).norm(dim=-1)
    assert float(model.median_increment) == pytest.approx(float(d.median()), rel=1e-6)
    law = ml.LSTMLaw(model)
    assert law.commit_tol == pytest.approx(0.1 * float(model.median_increment))


def test_committed_state_rule(eliso_lstm):
    model, _ = eliso_lstm
    tol = 1e-3
    law = ml.LSTMLaw(model, commit_tol=tol)
    h0 = np.zeros((law.state_size, 1))
    eps_c = np.zeros((6, 1))
    eps_c[0] = 2e-3
    # first commit from the origin
    st, Lt, h1, used = law.step_batch(eps_c, h0)
    assert np.abs(h1).max() > 0 and np.allclose(used[:, 0], eps_c[:, 0])
    # a move below the tolerance: trial stress, state and committed strain unchanged
    small = eps_c + 2e-4
    st2, Lt2, h2, used2 = law.step_batch(small, h1, strain_prev=used)
    np.testing.assert_array_equal(h2, h1)
    np.testing.assert_array_equal(used2, used)
    assert not np.allclose(st2, st)                    # but the stress does respond
    # a move above the tolerance commits
    big = eps_c + 2e-3
    st3, Lt3, h3, used3 = law.step_batch(big, h1, strain_prev=used)
    assert np.abs(h3 - h1).max() > 0 and np.allclose(used3[:, 0], big[law.idx_in, 0])
    # through the solver: statev (h, c, strain_c) is frozen across sub-tolerance increments
    law0 = ml.LSTMLaw(model, commit_tol=0.0)
    step = StepMeca(control=["strain"] * 6, value=[1e-3, 0, 0, 0, 0, 0], ninc=50)   # 2e-5 per inc
    r_tol = solve(step, law)
    r_0 = solve(step, law0)
    sv = r_tol["Statev"]
    n_commits = int((np.abs(np.diff(sv[:law.state_size], axis=1)).max(axis=0) > 0).sum())
    assert 0 < n_commits < 10                          # ~ path length / tol = 1e-3 / 1e-3 -> a few
    assert (np.abs(np.diff(r_0["Statev"][:law.state_size], axis=1)).max(axis=0) > 0).sum() == 49
    assert r_tol.status == 0 and r_0.status == 0


# ---------------------------------------------------------------------------
# point weights on the plastic increment, and the tangent-aware loss
# ---------------------------------------------------------------------------

EPICP_PROPS = [E, NU, 1.0e-5, 300.0, 1000.0, 0.3]      # E nu alpha sigma_Y k m
EPICP_NSTATEV = 8


def test_generate_dataset_records_statev_and_tangent():
    targets, ninc = ml.random_strain_paths(4, n_segments=2, n_sub=6, amplitude=0.01, seed=1)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc,
                             record=("statev", "tangent"))
    assert ds.extras["statev"].shape == (4, 12, 1)
    assert ds.extras["tangent"].shape == (4, 12, 6, 6)
    # elasticity: the recorded operator is L at every step
    L = sim.L_iso([E, NU], "Enu")
    np.testing.assert_allclose(ds.extras["tangent"].numpy(), np.broadcast_to(L, (4, 12, 6, 6)),
                               rtol=1e-5, atol=1e-2)
    with pytest.raises(ValueError, match="unknown record"):
        ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc, record=("nope",))
    # extras survive subset and resampling
    sub = ds.subset([0, 2])
    assert sub.extras["tangent"].shape == (2, 12, 6, 6)
    from simcoon.ml.data import resample_time
    x, y, mask, ex = resample_time(ds.x, ds.y, ds.mask, 7, ds.meta["features"], 6, ds.extras)
    assert ex["tangent"].shape == (4, 7, 6, 6) and ex["statev"].shape == (4, 7, 1)


def test_sequence_tangent_matches_the_law_tangent():
    """The batched block-diagonal tangent is the one the solver receives, step by step."""
    torch.manual_seed(0)
    targets, ninc = ml.random_strain_paths(3, n_segments=2, n_sub=5, amplitude=0.01, seed=2)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc)
    model = ml.StressLSTM(hidden_size=8, num_layers=2)
    model.fit_scalers(ds.x, ds.y, ds.mask)
    D = ml.sequence_tangent(model, ds.x).detach().numpy()
    assert D.shape == (3, 10, 6, 6)
    law = ml.LSTMLaw(model, commit_tol=0.0)
    for i in (0, 2):
        h = np.zeros((law.state_size, 1))
        prev = np.zeros((law.n_comp, 1))
        for t in range(4):
            eps6 = np.zeros((6, 1))
            eps6[law.idx_in, 0] = ds.x[i, t, :6].numpy()
            _, Lt, h, prev = law.step_batch(eps6, h, strain_prev=prev)
            np.testing.assert_allclose(Lt[np.ix_(law.idx_out, law.idx_in)][:, :, 0], D[i, t],
                                       rtol=1e-6, atol=1e-6 * np.abs(D).max())


@pytest.mark.parametrize("cell", ["lstm_dstrain", "lmsc"])
def test_sequence_tangent_sums_every_strain_block(cell):
    """Perturbing the end-of-increment strain moves the 'strain' and 'dstrain' features
    together: the diagnostic must add both contributions, or it reports the wrong operator
    for a model that reads the increment (and exactly zero for the LMSC, which reads only
    the increment)."""
    torch.manual_seed(0)
    targets, ninc = ml.random_strain_paths(3, n_segments=2, n_sub=5, amplitude=0.01, seed=2)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc,
                             features=("strain", "dstrain"))
    model = (ml.StressLSTM(features=("strain", "dstrain"), hidden_size=8, num_layers=2)
             if cell == "lstm_dstrain" else ml.LMSC(n_state=8, depth=2, width=16))
    model.fit_scalers(ds.x, ds.y, ds.mask)
    D = ml.sequence_tangent(model, ds.x).detach().numpy()
    assert D.shape == (3, 10, 6, 6) and np.abs(D).max() > 0.0
    law = ml.LSTMLaw(model, commit_tol=0.0)
    for i in (0, 2):
        st = np.zeros((law.state_size, 1))
        prev = np.zeros((law.n_comp, 1))
        for t in range(4):
            eps6 = np.zeros((6, 1))
            eps6[law.idx_in, 0] = ds.x[i, t, :6].numpy()
            _, Lt, st, prev = law.step_batch(eps6, st, strain_prev=prev)
            np.testing.assert_allclose(Lt[np.ix_(law.idx_out, law.idx_in)][:, :, 0], D[i, t],
                                       rtol=1e-6, atol=1e-6 * np.abs(D).max())


# ---------------------------------------------------------------------------
# LMSC: the properties the architecture guarantees by construction
# ---------------------------------------------------------------------------

def _lmsc(**kw):
    torch.manual_seed(0)
    m = ml.LMSC(n_state=kw.pop("n_state", 8), depth=kw.pop("depth", 2),
                width=kw.pop("width", 16), **kw).double()
    m.y_std.fill_(300.0)                    # a realistic stress scale
    return m


def test_lmsc_stationarity_is_exact():
    """A zero increment leaves state and stress strictly unchanged (Eq. 25 with nu = 0).

    This is the condition (Eq. 11 of Bonatti and Mohr) that a gated cell cannot satisfy.
    """
    m = _lmsc()
    chi = torch.randn(4, m.n_state, dtype=torch.float64) * 0.3
    with torch.no_grad():
        for _ in range(100):                # stationarity holds step after step
            chi2 = m.update(torch.zeros(4, 6, dtype=torch.float64), chi)
            assert float((chi2 - chi).abs().max()) == 0.0
            assert float((m.stress(chi2) - m.stress(chi)).abs().max()) == 0.0
            chi = chi2


def test_lmsc_is_rate_independent_and_self_consistent():
    m = _lmsc()
    chi = torch.randn(3, m.n_state, dtype=torch.float64) * 0.3
    d = torch.randn(3, 6, dtype=torch.float64) * 1e-3
    # no time enters the equations: the law cannot depend on the time step at all
    law = ml.LSTMLaw(m)
    assert law.commit_tol == 0.0            # no committed-state rule needed
    common = dict(Etot=np.zeros(6), sigma=np.zeros(6), statev=np.zeros(law.nstatev),
                  Wm=np.zeros(4), T=293.15, DT=0.0, ndi=3, start=True, tangent_mode=2)
    a = law.integrate(DEtot=d[0].numpy(), DTime=1.0, **common)
    b = law.integrate(DEtot=d[0].numpy(), DTime=10.0, **common)
    for x, y in zip(a[:3], b[:3]):
        np.testing.assert_array_equal(x, y)
    # self-consistency: at frozen coefficients, splitting an increment changes nothing
    with torch.no_grad():
        nu, n = m._split(d)
        alpha, beta = m._alpha_beta(chi, n)
        step = lambda c, v: c + torch.expm1(-v * alpha) * (c - beta)
        one, two = step(chi, nu), step(step(chi, nu / 2), nu / 2)
        assert float((one - two).abs().max()) < 1e-15

        # and on the full cell the refinement error decreases monotonically
        def rollout(total, k):
            c = chi.clone()
            for _ in range(k):
                c = m.update(total / k, c)
            return m.stress(c)
        ref = rollout(d * 20, 2048)
        errs = [float((rollout(d * 20, k) - ref).abs().max()) for k in (1, 2, 4, 8, 16, 32)]
    assert all(b < a for a, b in zip(errs, errs[1:]))
    assert errs[-1] < errs[0] / 20


def test_lmsc_zero_state_is_stress_free():
    """The output map has no bias (Eq. 26), so chi = 0 is exactly the virgin state."""
    m = _lmsc()
    with torch.no_grad():
        assert float(m.stress(torch.zeros(2, m.n_state, dtype=torch.float64)).abs().max()) == 0.0
    x = torch.randn(4, 5, m.n_in, dtype=torch.float64)
    m.fit_scalers(x, torch.randn(4, 5, 6, dtype=torch.float64) * 100 + 50)
    assert float(m.y_mean.abs().max()) == 0.0          # a mean would break the property
    assert float(m.y_std.min()) > 0.0


def test_lmsc_analytic_tangent_matches_finite_differences():
    m = _lmsc()
    chi = torch.randn(3, m.n_state, dtype=torch.float64) * 0.3
    d = torch.randn(3, 6, dtype=torch.float64) * 1e-3
    with torch.no_grad():
        J = m.analytic_tangent(m.build_inputs(torch.zeros(3, 6, dtype=torch.float64), d), chi)
    h = 1e-7
    fd = torch.zeros_like(J)
    with torch.no_grad():
        for j in range(6):
            e = torch.zeros(3, 6, dtype=torch.float64)
            e[:, j] = h
            fd[:, :, j] = (m.stress(m.update(d + e, chi)) - m.stress(m.update(d - e, chi))) / (2 * h)
    torch.testing.assert_close(J, fd, rtol=1e-6, atol=1e-6 * float(fd.abs().max()))
    # the law hands the solver that same operator, and falls back to the elastic
    # predictor where the increment - hence the loading direction - vanishes
    law = ml.LSTMLaw(m)
    eps6 = np.zeros((6, 3))
    eps6[:, :] = d.numpy().T
    st, Lt, s1, used = law.step_batch(eps6, chi.numpy().T)
    np.testing.assert_allclose(Lt[:, :, 0], J[0].numpy(), rtol=1e-10, atol=1e-8 * abs(J).max())
    st, Lt, s1, used = law.step_batch(np.zeros((6, 1)), chi[:1].numpy().T)
    np.testing.assert_allclose(Lt[:, :, 0], law.elastic_L, rtol=1e-12)


def test_lmsc_in_the_solver_state_and_rollback():
    m = _lmsc()
    law = ml.LSTMLaw(m)
    assert law.nstatev == m.n_state + 6                # no committed strain to carry beyond
    step = StepMeca(control=["strain"] * 6, value=[2e-3, -1e-3, 0.0, 1e-3, 0.0, 0.0], ninc=20)
    r = solve(step, law)
    assert r.status == 0 and len(r) == 20
    # replaying the recorded strain path reproduces the run exactly: all the history the
    # solver rewinds lives in statev
    s = np.zeros((law.state_size, 1))
    prev = np.zeros((law.n_comp, 1))
    for k in range(len(r)):
        eps = r["Strain"][:, k].reshape(6, 1)
        sig, _, s, prev = law.step_batch(eps, s, strain_prev=prev)
    np.testing.assert_allclose(sig[:, 0], r["Stress"][:, -1], rtol=1e-10, atol=1e-9)
    # mixed control: the analytic tangent drives the Newton loop
    r2 = solve(StepMeca(control=UNIAXIAL, value=[2e-3, 0, 0, 0, 0, 0], ninc=20), law)
    assert r2.status == 0
    np.testing.assert_allclose(r2["Stress"][1:, -1], 0.0, atol=1e-6)


@pytest.fixture(scope="module")
def eliso_lmsc():
    """Small LMSC trained a few seconds on ELISO paths (3D)."""
    torch.manual_seed(0)
    targets, ninc = ml.random_strain_paths(96, n_segments=4, n_sub=10, amplitude=0.01, seed=1)
    ds = ml.generate_dataset("ELISO", ELISO_PROPS, 1, targets=targets, ninc=ninc,
                             features=("strain", "dstrain"))
    tr, te = ml.split_dataset(ds, 0.25, seed=0)
    model = ml.LMSC(n_state=20, depth=3, width=40)
    trl, _ = ml.train(model, tr, te, epochs=500, batch_size=16, lr=5e-3, verbose=False, seed=0)
    assert trl[-1] < trl[0] * 1e-2
    return model, te


def test_lmsc_learns_and_condenses(eliso_lmsc):
    """A trained LMSC drives the solver and condenses the stress-free directions."""
    model, te = eliso_lmsc
    law = ml.LSTMLaw(model)
    r = solve(StepMeca(control=["strain"] * 6, value=[4e-3, -1e-3, 0, 2e-3, 0, 0], ninc=20), law)
    ref = solve(StepMeca(control=["strain"] * 6, value=[4e-3, -1e-3, 0, 2e-3, 0, 0], ninc=20),
                "ELISO", ELISO_PROPS, 1)
    assert r.status == 0
    assert np.abs(r["Stress"] - ref["Stress"]).max() / np.abs(ref["Stress"]).max() < 0.25
    # plane stress and uniaxial: the local Newton finds the free strains
    eps = np.array([[4e-3], [1e-3], [0.0], [2e-3], [0.0], [0.0]])
    s0 = np.zeros((law.state_size, 1))
    for ndi, free in ((2, [2]), (1, [1, 2])):
        st, Lt, _, used = law.step_batch(eps, s0, ndi=ndi)
        np.testing.assert_allclose(st[free, 0], 0.0, atol=1e-6)
        assert np.isfinite(Lt).all() and abs(Lt[0, 0, 0]) > 0.0


def test_lmsc_batch_and_persistence(tmp_path):
    m = _lmsc()
    law = ml.LSTMLaw(m)
    rng = np.random.default_rng(0)
    N = 4
    eps = rng.normal(scale=1e-3, size=(6, N))
    s0 = rng.normal(scale=0.2, size=(law.state_size, N))
    st, Lt, s1, used = law.step_batch(eps, s0)
    for i in range(N):                                  # batch == N single calls
        one = law.step_batch(eps[:, [i]], s0[:, [i]])
        np.testing.assert_allclose(one[0][:, 0], st[:, i], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(one[1][:, :, 0], Lt[:, :, i], rtol=1e-8, atol=1e-8)
    # save / load through the generic cell registry
    p = tmp_path / "lmsc.pt"
    m.save(p)
    m2 = ml.StateModel.load(p)
    assert isinstance(m2, ml.LMSC) and m2.hparams == m.hparams
    x = torch.randn(2, 4, m.n_in, dtype=torch.float64)
    torch.testing.assert_close(m2.double()(x)[0], m(x)[0])
