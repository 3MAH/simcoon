#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Densities of an orientation distribution function at the angles x. `peaks` is a sequence
//of dicts {number, method, mean, s_dev, width, ampl, params}; with radian = false, x and
//the peak angles are degrees.
pybind11::array_t<double> get_densities_ODF(const pybind11::array_t<double> &x, const pybind11::object &peaks, const bool &radian);

//The discretisation of a phase along an ODF is done in Python
//(simcoon.solver.micromechanics.discretize_odf) on these densities.

} //namespace simpy
