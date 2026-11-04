"""Options of simcoon.identify.identification (least squares, threads), calc_residuals,
and the in-memory key templates of simcoon.parameter."""

import os

import numpy as np
import pytest
from scipy.optimize import differential_evolution

from simcoon.constant import Constant
from simcoon.identify import calc_cost, calc_residuals, identification
from simcoon.parameter import KeyTemplate, Parameter

T = np.linspace(0.0, 1.0, 40)
A_REF, B_REF = 3.0, 250.0


def decay_residuals(x):
    """Residuals of y = a exp(-b t / 100) against the reference curve."""
    return x[0] * np.exp(-x[1] * T / 100.0) - A_REF * np.exp(-B_REF * T / 100.0)


def decay_params():
    return [Parameter(0, bounds=(0.1, 10.0)), Parameter(1, bounds=(1.0, 1000.0))]


class TestLeastSquares:
    def test_recovers_parameters_and_writes_back(self):
        params = decay_params()
        result = identification(decay_residuals, params, method="least_squares")
        assert result.x == pytest.approx([A_REF, B_REF], rel=1e-8)
        assert [p.value for p in params] == pytest.approx([A_REF, B_REF], rel=1e-8)
        assert np.sum(result.fun ** 2) == pytest.approx(2 * result.cost)
        assert result.cost < 1e-16
        assert np.all((result.z >= 0) & (result.z <= 1))
        assert result.nfev >= len(result.starts[0].x0) + 1
        assert len(result.starts) == 1

    def test_first_start_is_middle_of_scaled_range(self):
        params = [Parameter(0, bounds=(1.0, 100.0)), Parameter(1, bounds=(-1.0, 1.0))]
        result = identification(lambda x: x - [20.0, 0.3], params, method="least_squares")
        # log scale (positive bounds): geometric mean; linear otherwise: midpoint
        assert result.starts[0].x0 == pytest.approx([10.0, 0.0])
        assert result.x == pytest.approx([20.0, 0.3])

    def test_first_start_uses_set_value(self):
        params = decay_params()
        params[1].value = 300.0
        result = identification(decay_residuals, params, method="least_squares")
        assert result.starts[0].x0 == pytest.approx([np.sqrt(0.1 * 10.0), 300.0])

    def test_multi_start_reproducible(self):
        r1 = identification(decay_residuals, decay_params(), method="least_squares",
                            n_starts=4, seed=3)
        r2 = identification(decay_residuals, decay_params(), method="least_squares",
                            n_starts=4, seed=3)
        assert len(r1.starts) == 4
        assert [s.x0.tolist() for s in r1.starts] == [s.x0.tolist() for s in r2.starts]
        assert r1.cost == min(s.cost for s in r1.starts)
        assert r1.nfev == sum(s.nfev for s in r1.starts)
        for s in r1.starts:
            lb, ub = np.array([0.1, 1.0]), np.array([10.0, 1000.0])
            assert np.all((s.x0 >= lb) & (s.x0 <= ub))

    @pytest.mark.parametrize("scale", ["linear", "log", ["log", "linear"]])
    def test_scales(self, scale):
        result = identification(decay_residuals, decay_params(), method="least_squares",
                                scale=scale)
        assert result.x == pytest.approx([A_REF, B_REF], rel=1e-6)

    def test_threads_do_not_change_the_result(self):
        r1 = identification(decay_residuals, decay_params(), method="least_squares",
                            n_starts=3, seed=0)
        r4 = identification(decay_residuals, decay_params(), method="least_squares",
                            n_starts=3, seed=0, n_threads=4)
        assert np.array_equal(r1.x, r4.x)
        assert r1.nfev == r4.nfev

    def test_args_forwarded(self):
        def res(x, target, scale):
            return scale * (x - target)
        params = [Parameter(0, bounds=(-5.0, 5.0))]
        result = identification(res, params, method="least_squares", args=([1.5], 2.0))
        assert result.x == pytest.approx([1.5])

    def test_nan_residuals_reject_the_step(self):
        """A model failing (NaN) in part of the box does not stop the optimizer."""
        def res(x):
            if x[1] > 600.0:
                return np.full(T.size, np.nan)
            return decay_residuals(x)
        params = decay_params()
        params[1].value = 550.0
        result = identification(res, params, method="least_squares")
        assert result.x == pytest.approx([A_REF, B_REF], rel=1e-6)

    def test_failed_start_is_skipped(self):
        def res(x):
            return np.full(3, np.nan) if x[0] > 9.0 else np.array([x[0] - 2.0, 0.0, 0.0])
        params = [Parameter(0, bounds=(1.0, 10.0))]
        params[0].value = 9.5
        result = identification(res, params, method="least_squares", n_starts=2, seed=0)
        assert result.starts[0].status == -1
        assert result.x == pytest.approx([2.0])

    def test_all_starts_failing_raises(self):
        params = [Parameter(0, bounds=(1.0, 10.0))]
        with pytest.raises(RuntimeError, match="every starting point"):
            identification(lambda x: np.full(3, np.nan), params, method="least_squares")

    def test_invalid_options(self):
        with pytest.raises(ValueError, match="method"):
            identification(decay_residuals, decay_params(), method="nelder")
        with pytest.raises(ValueError, match="positive lower bound"):
            identification(decay_residuals, [Parameter(0, bounds=(-1.0, 1.0))],
                           method="least_squares", scale="log")
        with pytest.raises(ValueError, match="scale"):
            identification(decay_residuals, decay_params(), method="least_squares",
                           scale="exp")
        with pytest.raises(ValueError, match="n_starts"):
            identification(decay_residuals, decay_params(), method="least_squares",
                           n_starts=0)
        with pytest.raises(TypeError, match="jac"):
            identification(decay_residuals, decay_params(), method="least_squares",
                           jac="3-point")
        with pytest.raises(ValueError, match="n_threads"):
            identification(decay_residuals, decay_params(), n_threads=-1)


class TestDifferentialEvolutionOptions:
    @staticmethod
    def cost(x):
        return float(np.sum(decay_residuals(x) ** 2))

    def test_default_is_unchanged(self):
        """No option: the plain scipy call with polish=True."""
        params = decay_params()
        result = identification(self.cost, params, seed=1, maxiter=30)
        ref = differential_evolution(self.cost, [p.bounds for p in decay_params()],
                                     seed=1, maxiter=30, polish=True)
        assert np.array_equal(result.x, ref.x)
        assert result.nfev == ref.nfev

    def test_threads_match_deferred_serial(self):
        r1 = identification(self.cost, decay_params(), seed=1, maxiter=30,
                            updating="deferred")
        r3 = identification(self.cost, decay_params(), seed=1, maxiter=30, n_threads=3)
        assert np.array_equal(r1.x, r3.x)
        assert r1.nfev == r3.nfev

    def test_threads_and_workers_exclusive(self):
        with pytest.raises(ValueError, match="workers"):
            identification(self.cost, decay_params(), n_threads=2, workers=2)


class TestSolverIdentificationThreads:
    """End to end on the simcoon solver, whose solve releases the GIL."""

    @staticmethod
    def setup_case():
        import simcoon as sim
        from simcoon.solver import Block, StepMeca

        block = Block(steps=[StepMeca(control=["strain"] + ["stress"] * 5,
                                      value=[0.02, 0, 0, 0, 0, 0], ninc=20)])

        def stress(x):
            props = np.array([70000.0, 0.3, 0.0, x[0], x[1], 0.3])
            res = sim.solver.solve(block, "EPICP", props, 8, T_init=293.15)
            return np.asarray(res["Stress"][0])

        sigma_exp = stress([300.0, 1000.0])

        def residuals(x):
            return calc_residuals([sigma_exp], [stress(x)], metric="nmse_per_response")

        return residuals

    def test_epicp_serial_and_threads_identical(self):
        residuals = self.setup_case()
        params = [Parameter(0, bounds=(100.0, 600.0)), Parameter(1, bounds=(100.0, 5000.0))]
        r1 = identification(residuals, params, method="least_squares")
        params = [Parameter(0, bounds=(100.0, 600.0)), Parameter(1, bounds=(100.0, 5000.0))]
        r2 = identification(residuals, params, method="least_squares", n_threads=2)
        assert np.array_equal(r1.x, r2.x)
        assert r1.x == pytest.approx([300.0, 1000.0], rel=1e-4)


class TestHomogenizationIdentificationThreads:
    """L_eff releases the GIL and its quadrature tables are per (mp, np): a
    threaded identification of a self-consistent composite equals the serial one."""

    @staticmethod
    def residuals(x, c_exp, E_exp):
        import simcoon as sim
        from simcoon.solver.micromechanics import Ellipsoid

        matrix = Ellipsoid(number=0, umat_name="ELISO", save=1, nstatev=1,
                           props=np.array([2250.0, 0.19, 8.8e-5]))
        fibre = Ellipsoid(number=1, umat_name="ELISO", save=1, nstatev=1,
                          props=np.array([x[0], x[1], 0.5e-6]))
        E = []
        for c in c_exp:
            matrix.concentration, fibre.concentration = 1.0 - c, c
            L = sim.L_eff("MISCN", np.array([10.0, 10.0, 0.0]), 0, phases=[matrix, fibre])
            E.append(sim.L_iso_props(L).flatten()[0])
        return np.array(E) - E_exp

    def test_miscn_serial_and_threads_identical(self):
        c_exp = np.array([0.1, 0.2, 0.3, 0.4])
        E_exp = self.residuals([60000.0, 0.25], c_exp, np.zeros(4))
        results = []
        for nt in (1, 2):
            params = [Parameter(0, bounds=(10000.0, 200000.0)),
                      Parameter(1, bounds=(0.01, 0.45))]
            results.append(identification(self.residuals, params, method="least_squares",
                                          args=(c_exp, E_exp), n_threads=nt))
        assert np.array_equal(results[0].x, results[1].x)
        assert results[0].x[0] == pytest.approx(60000.0, rel=1e-4)


class TestCalcResiduals:
    @staticmethod
    def data():
        rng = np.random.default_rng(0)
        y_exp = [rng.normal(size=(20, 2)) * [100.0, 0.1],
                 rng.normal(size=(15, 2)) * [50.0, 0.05],
                 rng.normal(size=10)]
        y_num = [y + rng.normal(size=y.shape) * 0.1 * np.abs(y).mean(axis=0) for y in y_exp]
        weights = dict(
            w_test=np.array([1.0, 2.0, 0.5]),
            w_response=[np.array([1.0, 3.0]), np.array([2.0, 1.0]), np.array([1.0])],
            w_point=[rng.uniform(size=y.shape) for y in y_exp[:2]]
            + [rng.uniform(size=(10, 1))],
        )
        return y_exp, y_num, weights

    @pytest.mark.parametrize("metric", ["mse", "nmse", "nmse_per_response"])
    @pytest.mark.parametrize("weighted", [False, True])
    def test_sum_of_squares_is_calc_cost(self, metric, weighted):
        y_exp, y_num, weights = self.data()
        kw = weights if weighted else {}
        r = calc_residuals(y_exp, y_num, metric=metric, **kw)
        assert r.ndim == 1 and r.size == 20 * 2 + 15 * 2 + 10
        assert np.sum(r ** 2) == pytest.approx(calc_cost(y_exp, y_num, metric=metric, **kw),
                                               rel=1e-12)

    def test_sign_is_model_minus_experiment(self):
        r = calc_residuals([np.array([1.0, 2.0])], [np.array([2.0, 2.0])])
        assert r[0] > 0 and r[1] == 0

    def test_unsupported_metric(self):
        with pytest.raises(ValueError, match="sum-of-squares"):
            calc_residuals([np.ones(3)], [np.ones(3)], metric="mae")

    def test_shape_mismatch(self):
        with pytest.raises(ValueError, match="shape"):
            calc_residuals([np.ones(3)], [np.ones(4)])
        with pytest.raises(ValueError, match="tests"):
            calc_residuals([np.ones(3)], [np.ones(3), np.ones(3)])


class TestKeyTemplate:
    @staticmethod
    def setup(tmp_path):
        keys = tmp_path / "keys"
        keys.mkdir()
        (keys / "mat.json").write_text('{"props": [@E, @Et, @nu, @m]}')
        (keys / "other.dat").write_text("E=@E m=@m mu=@mu")
        params = [
            Parameter(0, bounds=(0, 1), key="@E", sim_input_files=["mat.json", "other.dat"]),
            Parameter(1, bounds=(0, 1), key="@Et", sim_input_files=["mat.json"]),
            Parameter(2, bounds=(0, 1), key="@nu", sim_input_files=["mat.json"]),
        ]
        consts = [
            Constant(0, "@m", np.array([0.3]), 0.3, ["mat.json", "other.dat"]),
            Constant(1, "@mu", np.array([0.5]), 0.5, ["other.dat"]),
        ]
        for p, v in zip(params, (70000.0, 2000.0, 0.3)):
            p.value = v
        return keys, params, consts

    def test_parameters_then_constants_per_file(self, tmp_path):
        keys, params, consts = self.setup(tmp_path)
        tpl = tmp_path / "tpl"
        tpl.mkdir()
        KeyTemplate(keys, params, consts).write(tpl)
        assert sorted(os.listdir(tpl)) == ["mat.json", "other.dat"]
        assert (tpl / "mat.json").read_text() == '{"props": [70000.0, 2000.0, 0.3, 0.3]}'
        assert (tpl / "other.dat").read_text() == "E=70000.0 m=0.3 mu=0.5"

    def test_templates_read_once(self, tmp_path):
        keys, params, consts = self.setup(tmp_path)
        template = KeyTemplate(keys, params, consts)
        (keys / "mat.json").write_text("changed on disk")
        out = tmp_path / "out"
        out.mkdir()
        params[0].value = 1.5
        template.write(out)
        assert (out / "mat.json").read_text() == '{"props": [1.5, 2000.0, 0.3, 0.3]}'
        assert template.render("other.dat") == "E=1.5 m=0.3 mu=0.5"

    def test_invalid_paths(self, tmp_path):
        keys, params, consts = self.setup(tmp_path)
        with pytest.raises(TypeError):
            KeyTemplate(123, params)
        with pytest.raises(TypeError):
            KeyTemplate(keys, params).write(123)
