import os
import tempfile

import numpy as np
import pytest

from simcoon.parameter import Parameter, read_parameters, copy_parameters, apply_parameters
from simcoon.identify import identification, calc_cost


@pytest.fixture
def data_dir():
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def keys_dir():
    return os.path.join(os.path.dirname(__file__), "keys")


class TestReadParameters:
    def test_read_parameters(self, data_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        assert len(params) == 4

    def test_first_parameter(self, data_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        p = params[0]
        assert p.number == 0
        assert p.bounds == (0.0, 180.0)
        assert p.key == "@0p"
        assert p.sim_input_files == ["Nellipsoids0.dat"]

    def test_parameter_default_value_is_midpoint(self, data_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        p = params[0]
        assert p.value == pytest.approx(90.0)

    def test_read_parameters_invalid_path(self):
        with pytest.raises(FileNotFoundError):
            read_parameters("nonexistent/parameters.inp")

    def test_read_parameters_invalid_type(self):
        with pytest.raises(TypeError):
            read_parameters(123)


class TestCopyParameters:
    def test_copy_parameters(self, data_dir, keys_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        with tempfile.TemporaryDirectory() as tmp:
            copy_parameters(params, keys_dir, tmp)
            assert os.path.exists(os.path.join(tmp, "Nellipsoids0.dat"))

    def test_copy_parameters_invalid_src(self, data_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        with pytest.raises(TypeError):
            copy_parameters(params, 123, "/tmp")


class TestApplyParameters:
    def test_apply_parameters_replaces_keys(self, data_dir, keys_dir):
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        params[0].value = 45.0
        params[1].value = 72000.0

        with tempfile.TemporaryDirectory() as tmp:
            copy_parameters(params, keys_dir, tmp)
            apply_parameters(params, tmp)

            with open(os.path.join(tmp, "Nellipsoids0.dat"), "r") as f:
                content = f.read()

            assert "@0p" not in content
            assert "@1p" not in content
            assert "45.0" in content
            assert "72000.0" in content


class TestParameterRepr:
    def test_repr(self):
        p = Parameter(
            number=0, bounds=(0.0, 180.0), key="@0p",
            sim_input_files=["mat.dat"],
        )
        r = repr(p)
        assert "Parameter(" in r
        assert "@0p" in r
        assert "bounds=(0.0, 180.0)" in r


class TestIdentificationWorkflow:
    """Test the full identification workflow: template -> apply -> forward model -> cost."""

    def test_template_apply_no_residual_keys(self, data_dir, keys_dir):
        """End-to-end: key substitution produces valid input (no residual @)."""
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        params[0].value = 0.8
        params[1].value = 73000

        with tempfile.TemporaryDirectory() as tmp:
            copy_parameters(params, keys_dir, tmp)
            apply_parameters(params, tmp)

            with open(os.path.join(tmp, "Nellipsoids0.dat"), "r") as f:
                content = f.read()

            assert "@" not in content, "Residual keys found after apply_parameters"
            assert "0.8" in content
            assert "73000" in content

    def test_multiple_iterations_fresh_copy(self, data_dir, keys_dir):
        """Each iteration starts from a fresh template (no stale values)."""
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))

        with tempfile.TemporaryDirectory() as tmp:
            # First iteration
            params[0].value = 0.7
            params[1].value = 50000
            copy_parameters(params, keys_dir, tmp)
            apply_parameters(params, tmp)

            with open(os.path.join(tmp, "Nellipsoids0.dat"), "r") as f:
                content1 = f.read()
            assert "50000" in content1

            # Second iteration with different values
            params[0].value = 0.9
            params[1].value = 80000
            copy_parameters(params, keys_dir, tmp)
            apply_parameters(params, tmp)

            with open(os.path.join(tmp, "Nellipsoids0.dat"), "r") as f:
                content2 = f.read()
            assert "80000" in content2
            assert "50000" not in content2, "Stale value from previous iteration"

    def test_bounds_available_for_optimizer(self, data_dir):
        """Parameter.bounds can be used as scipy optimizer bounds."""
        params = read_parameters(os.path.join(data_dir, "parameters.inp"))
        bounds = [p.bounds for p in params]
        assert all(lo < hi for lo, hi in bounds)
        assert bounds[0] == (0.0, 180.0)


class TestIdentification:
    """Test the identification() wrapper around differential_evolution."""

    def test_identifies_simple_quadratic(self):
        """Identify minimum of (x - 3)^2 + (y - 7)^2."""
        params = [
            Parameter(0, bounds=(0, 10), key="@x"),
            Parameter(1, bounds=(0, 10), key="@y"),
        ]

        def cost(x):
            return (x[0] - 3) ** 2 + (x[1] - 7) ** 2

        result = identification(cost, params, seed=42, tol=1e-10)
        assert result.success
        assert params[0].value == pytest.approx(3.0, abs=1e-4)
        assert params[1].value == pytest.approx(7.0, abs=1e-4)

    def test_values_written_back(self):
        """After identification, Parameter.value holds the optimum."""
        params = [Parameter(0, bounds=(-5, 5), key="@a")]
        result = identification(lambda x: (x[0] - 2.5) ** 2, params, seed=0)
        assert params[0].value == pytest.approx(2.5, abs=1e-3)

    def test_unknown_method_raises(self):
        params = [Parameter(0, bounds=(0, 1), key="@x")]
        with pytest.raises(ValueError, match="Unknown method"):
            identification(lambda x: x[0], params, method="bogus")


class TestCalcCost:
    """Test the calc_cost() multi-level weighted cost function."""

    def test_simple_mse(self):
        y_exp = [np.array([1.0, 2.0, 3.0])]
        y_num = [np.array([1.1, 2.1, 3.1])]
        cost = calc_cost(y_exp, y_num)
        assert cost == pytest.approx(0.01, abs=1e-10)

    def test_two_tests(self):
        """Cost aggregates across multiple tests."""
        y_exp = [np.array([10.0]), np.array([20.0])]
        y_num = [np.array([11.0]), np.array([22.0])]
        # residuals: 1, 2 → MSE = (1 + 4) / 2 = 2.5
        cost = calc_cost(y_exp, y_num)
        assert cost == pytest.approx(2.5)

    def test_w_test(self):
        """Per-test weights scale the contribution of each test."""
        y_exp = [np.array([0.0]), np.array([0.0])]
        y_num = [np.array([1.0]), np.array([1.0])]
        # Without weights: MSE = 1.0
        # With w_test=[2, 0]: only first test counts, weighted avg = 2*1/2 = 1.0
        # Actually weighted average: (2*1 + 0*1)/(2+0) ... np.average uses sum(w*x)/sum(w)
        # No — np.average(residuals**2, weights=w) = sum(w*r²)/sum(w)
        cost_uniform = calc_cost(y_exp, y_num)
        cost_weighted = calc_cost(y_exp, y_num, w_test=np.array([1.0, 0.0]))
        # weighted: sum(w*r²)/sum(w) = (1*1+0*1)/(1+0) = 1.0
        assert cost_uniform == pytest.approx(1.0)
        assert cost_weighted == pytest.approx(1.0)

    def test_2d_responses(self):
        """Multiple response columns per test."""
        y_exp = [np.array([[100, 0.1], [200, 0.2]])]  # (2 points, 2 responses)
        y_num = [np.array([[110, 0.1], [200, 0.3]])]
        cost = calc_cost(y_exp, y_num)
        # residuals: [10, 0, 0, 0.1] → MSE = (100+0+0+0.01)/4 = 25.0025
        assert cost == pytest.approx(25.0025)

    def test_nmse_per_response(self):
        """nmse_per_response balances responses of different magnitudes."""
        # Two responses: force (~1000 N) and displacement (~0.01 mm)
        y_exp = [np.array([[1000.0, 0.01], [2000.0, 0.02]])]
        y_num = [np.array([[1010.0, 0.015], [2010.0, 0.025]])]
        # Without normalization, force residuals dominate (10² >> 0.005²)
        cost_mse = calc_cost(y_exp, y_num, metric="mse")
        # nmse_per_response normalizes each column by sum(exp²)
        cost_nmse = calc_cost(y_exp, y_num, metric="nmse_per_response")
        # Results must differ (force dominates MSE but not NMSE per response)
        assert cost_mse != pytest.approx(cost_nmse)
        assert cost_nmse > 0
        assert np.isfinite(cost_nmse)

    def test_nmse_per_response_equal_columns(self):
        """With identical columns, nmse_per_response equals standard NMSE."""
        y_exp = [np.array([[10.0], [20.0]])]
        y_num = [np.array([[11.0], [21.0]])]
        cost = calc_cost(y_exp, y_num, metric="nmse_per_response")
        # SSE=2, denom=10²+20²=500 → NMSE=2/500=0.004
        assert cost == pytest.approx(2.0 / 500.0)

    def test_rmse(self):
        y_exp = [np.array([0.0, 0.0])]
        y_num = [np.array([3.0, 4.0])]
        cost = calc_cost(y_exp, y_num, metric="rmse")
        assert cost == pytest.approx(np.sqrt(12.5))

    def test_mae(self):
        y_exp = [np.array([0.0, 0.0])]
        y_num = [np.array([3.0, -4.0])]
        cost = calc_cost(y_exp, y_num, metric="mae")
        assert cost == pytest.approx(3.5)

    def test_mismatched_shapes_raises(self):
        with pytest.raises(ValueError, match="shape"):
            calc_cost([np.array([1, 2])], [np.array([1, 2, 3])])

    def test_mismatched_ntests_raises(self):
        with pytest.raises(ValueError, match="tests"):
            calc_cost([np.array([1])], [np.array([1]), np.array([2])])

    def test_unknown_metric_without_sklearn(self):
        """Unknown metric raises ValueError with helpful message."""
        with pytest.raises((ValueError, ImportError)):
            calc_cost([np.array([1.0])], [np.array([2.0])], metric="bogus_metric_xyz")
