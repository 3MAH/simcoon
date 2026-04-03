import os
import shutil
import tempfile

import pytest

from simcoon.parameter import Parameter, read_parameters, copy_parameters, apply_parameters


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

            # Keys should be replaced by values
            assert "@0p" not in content
            assert "@1p" not in content
            assert "45.0" in content
            assert "72000.0" in content


class TestParameterRepr:
    def test_repr(self):
        p = Parameter(number=0, bounds=(0.0, 180.0), key="@0p", sim_input_files=["mat.dat"])
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

    def test_deprecation_stubs(self):
        """Removed functions raise NotImplementedError with migration info."""
        import simcoon

        with pytest.raises(NotImplementedError, match="scipy"):
            simcoon.identification()

        with pytest.raises(NotImplementedError, match="removed"):
            simcoon.calc_cost()
