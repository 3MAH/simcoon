#pragma once

#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace simpy {

//In-memory solver: takes loading blocks as a list of dicts, returns results as a dict of numpy arrays
py::dict solver_run(const py::list &blocks_py, const double &T_init,
                    const std::string &umat_name, const py::array_t<double> &props_py,
                    const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve,
                    const int &solver_type, const int &corate_type,
                    const py::dict &params_py, const bool &record_tangent);

} //namespace simpy
