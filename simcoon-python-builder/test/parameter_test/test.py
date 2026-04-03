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
