#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy {
	pybind11::tuple launch_umat(const std::string& umat_name_py, const pybind11::array_t<double> &etot_py, const pybind11::array_t<double> &Detot_py, const pybind11::array_t<double> &sigma_py, const pybind11::array_t<double> &DR_py, const pybind11::array_t<double> &props_py, const pybind11::array_t<double> &statev_py, const float Time, const float DTime, const pybind11::array_t<double> &Wm_py);
}
