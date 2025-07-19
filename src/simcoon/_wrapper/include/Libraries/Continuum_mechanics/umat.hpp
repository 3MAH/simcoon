#pragma once
#include <optional>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py=pybind11;
using namespace std;

namespace simpy {
	py::tuple launch_umat(const std::string &umat_name_py, const py::array_t<double> &etot_py, const py::array_t<double> &Detot_py, const py::array_t<double> &F0_py, const py::array_t<double> &F1_py, const py::array_t<double> &sigma_py, const py::array_t<double> &DR_py, const py::array_t<double> &props_py, const py::array_t<double> &statev_py, const float Time, const float DTime, const py::array_t<double> &Wm_py, const py::object &T_py, const int &ndi, const unsigned int &n_threads);
}
