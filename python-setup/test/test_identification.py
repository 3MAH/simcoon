"""
Tests for the simcoon.identification module.

These tests verify the identification/calibration functionality
including problem definition, optimizers, cost functions, and
sensitivity analysis.
"""

import pytest
import numpy as np


class TestParameterSpec:
    """Tests for ParameterSpec class."""

    def test_parameter_creation(self):
        """Test creating a parameter specification."""
        from simcoon.identification import ParameterSpec

        param = ParameterSpec(
            name='E',
            bounds=(100000, 300000),
            initial=200000,
        )

        assert param.name == 'E'
        assert param.bounds == (100000, 300000)
        assert param.initial == 200000
        assert param.fixed == False

    def test_parameter_default_initial(self):
        """Test that initial defaults to midpoint of bounds."""
        from simcoon.identification import ParameterSpec

        param = ParameterSpec(name='nu', bounds=(0.2, 0.4))
        assert param.initial == 0.3

    def test_normalization(self):
        """Test parameter normalization."""
        from simcoon.identification import ParameterSpec

        param = ParameterSpec(name='E', bounds=(100000, 200000))

        # Midpoint should normalize to 0.5
        assert param.normalize(150000) == 0.5

        # Denormalize back
        assert param.denormalize(0.5) == 150000

        # Bounds
        assert param.normalize(100000) == 0.0
        assert param.normalize(200000) == 1.0


class TestCostFunctions:
    """Tests for cost functions."""

    def test_mse(self):
        """Test mean squared error."""
        from simcoon.identification import mse

        y_true = np.array([1.0, 2.0, 3.0, 4.0])
        y_pred = np.array([1.0, 2.0, 3.0, 4.0])

        # Perfect prediction
        assert mse(y_true, y_pred) == 0.0

        # Known error
        y_pred_err = np.array([2.0, 2.0, 3.0, 4.0])
        assert mse(y_true, y_pred_err) == 0.25  # (1^2) / 4

    def test_mae(self):
        """Test mean absolute error."""
        from simcoon.identification import mae

        y_true = np.array([1.0, 2.0, 3.0, 4.0])
        y_pred = np.array([1.5, 2.5, 3.5, 4.5])

        assert mae(y_true, y_pred) == 0.5

    def test_r2(self):
        """Test R-squared (returns 1-R^2 for minimization)."""
        from simcoon.identification import r2

        y_true = np.array([1.0, 2.0, 3.0, 4.0])
        y_pred = np.array([1.0, 2.0, 3.0, 4.0])

        # Perfect prediction: R^2 = 1, so 1 - R^2 = 0
        assert r2(y_true, y_pred) == pytest.approx(0.0)

    def test_huber_loss(self):
        """Test Huber loss."""
        from simcoon.identification import huber_loss

        y_true = np.array([0.0, 0.0, 0.0])
        y_pred = np.array([0.5, 0.5, 10.0])  # Small, small, outlier

        # Huber should be robust to the outlier
        loss = huber_loss(y_true, y_pred, delta=1.0)
        assert loss > 0


class TestIdentificationProblem:
    """Tests for IdentificationProblem class."""

    def test_problem_creation(self):
        """Test creating an identification problem."""
        from simcoon.identification import IdentificationProblem

        def dummy_simulate(params):
            return {'output': params[0] * np.array([1, 2, 3])}

        problem = IdentificationProblem(
            parameters=[
                {'name': 'a', 'bounds': (0, 10)},
                {'name': 'b', 'bounds': (0, 1), 'fixed': True, 'initial': 0.5},
            ],
            simulate=dummy_simulate,
            exp_data={'output': np.array([5, 10, 15])},
        )

        # Only non-fixed parameters count
        assert problem.n_params == 1
        assert problem.parameter_names == ['a']

    def test_cost_function(self):
        """Test cost function evaluation."""
        from simcoon.identification import IdentificationProblem

        # Linear model: y = a * x
        def simulate(params):
            a = params[0]
            x = np.array([1, 2, 3, 4])
            return {'y': a * x}

        exp_data = {'y': np.array([2, 4, 6, 8])}  # True a = 2

        problem = IdentificationProblem(
            parameters=[{'name': 'a', 'bounds': (0, 5)}],
            simulate=simulate,
            exp_data=exp_data,
        )

        # Cost at true value should be ~0
        cost_true = problem.cost_function(np.array([2.0]))
        assert cost_true < 1e-10

        # Cost at wrong value should be > 0
        cost_wrong = problem.cost_function(np.array([3.0]))
        assert cost_wrong > 0

    def test_residual_vector(self):
        """Test residual vector computation."""
        from simcoon.identification import IdentificationProblem

        def simulate(params):
            return {'y': params[0] * np.array([1, 2])}

        problem = IdentificationProblem(
            parameters=[{'name': 'a', 'bounds': (0, 5)}],
            simulate=simulate,
            exp_data={'y': np.array([1, 2])},  # True a = 1
        )

        # Residuals at true value
        residuals = problem.residual_vector(np.array([1.0]))
        np.testing.assert_array_almost_equal(residuals, [0, 0])

        # Residuals at a = 2
        residuals = problem.residual_vector(np.array([2.0]))
        np.testing.assert_array_almost_equal(residuals, [1, 2])


class TestOptimizers:
    """Tests for optimization algorithms."""

    @pytest.fixture
    def simple_problem(self):
        """Create a simple quadratic optimization problem."""
        from simcoon.identification import IdentificationProblem

        # Minimize (a - 3)^2 + (b - 5)^2
        def simulate(params):
            a, b = params
            return {
                'output': np.array([a, b])
            }

        return IdentificationProblem(
            parameters=[
                {'name': 'a', 'bounds': (0, 10), 'initial': 1.0},
                {'name': 'b', 'bounds': (0, 10), 'initial': 1.0},
            ],
            simulate=simulate,
            exp_data={'output': np.array([3.0, 5.0])},
        )

    def test_levenberg_marquardt(self, simple_problem):
        """Test Levenberg-Marquardt optimizer."""
        from simcoon.identification import levenberg_marquardt

        result = levenberg_marquardt(simple_problem)

        assert result.success
        np.testing.assert_array_almost_equal(result.x, [3.0, 5.0], decimal=3)

    def test_nelder_mead(self, simple_problem):
        """Test Nelder-Mead optimizer."""
        from simcoon.identification import nelder_mead

        result = nelder_mead(simple_problem)

        assert result.success
        np.testing.assert_array_almost_equal(result.x, [3.0, 5.0], decimal=2)

    @pytest.mark.slow
    def test_differential_evolution(self, simple_problem):
        """Test differential evolution optimizer."""
        from simcoon.identification import differential_evolution

        result = differential_evolution(
            simple_problem,
            maxiter=50,
            seed=42,
        )

        assert result.success
        np.testing.assert_array_almost_equal(result.x, [3.0, 5.0], decimal=1)


class TestSensitivity:
    """Tests for sensitivity analysis."""

    @pytest.fixture
    def linear_problem(self):
        """Create a linear problem for sensitivity testing."""
        from simcoon.identification import IdentificationProblem

        # y = a * x + b
        def simulate(params):
            a, b = params
            x = np.linspace(0, 1, 10)
            return {'y': a * x + b}

        return IdentificationProblem(
            parameters=[
                {'name': 'a', 'bounds': (0, 10)},
                {'name': 'b', 'bounds': (0, 10)},
            ],
            simulate=simulate,
            exp_data={'y': np.linspace(0, 1, 10)},
        )

    def test_compute_sensitivity(self, linear_problem):
        """Test sensitivity computation."""
        from simcoon.identification import compute_sensitivity

        params = np.array([1.0, 0.5])
        sens = compute_sensitivity(linear_problem, params)

        assert 'y' in sens
        assert sens['y'].shape == (10, 2)

    def test_compute_jacobian(self, linear_problem):
        """Test Jacobian computation."""
        from simcoon.identification import compute_jacobian

        params = np.array([1.0, 0.5])
        jac = compute_jacobian(linear_problem, params)

        # Jacobian should have shape (n_residuals, n_params)
        assert jac.shape[1] == 2

    def test_correlation_matrix(self, linear_problem):
        """Test correlation matrix computation."""
        from simcoon.identification import correlation_matrix

        params = np.array([1.0, 0.5])
        corr = correlation_matrix(linear_problem, params)

        # Correlation matrix should be 2x2, symmetric
        assert corr.shape == (2, 2)
        np.testing.assert_array_almost_equal(corr, corr.T)

        # Diagonal should be 1
        np.testing.assert_array_almost_equal(np.diag(corr), [1, 1])


class TestOptimizationResult:
    """Tests for OptimizationResult class."""

    def test_from_scipy(self):
        """Test creating result from scipy output."""
        from simcoon.identification import OptimizationResult

        # Mock scipy result
        class MockResult:
            x = np.array([1.0, 2.0])
            fun = 0.5
            success = True
            message = "Converged"
            nit = 10
            nfev = 50
            jac = np.eye(2)

        result = OptimizationResult.from_scipy(MockResult(), ['a', 'b'])

        assert np.allclose(result.x, [1.0, 2.0])
        assert result.cost == 0.5
        assert result.success == True
        assert result.n_iterations == 10

    def test_repr(self):
        """Test string representation."""
        from simcoon.identification import OptimizationResult

        result = OptimizationResult(
            x=np.array([1.0, 2.0]),
            cost=0.5,
            success=True,
            message="OK",
            n_iterations=10,
            n_function_evals=50,
            parameter_names=['a', 'b'],
        )

        repr_str = repr(result)
        assert 'success=True' in repr_str
        assert 'a=' in repr_str


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
