#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function returns the isochoric invariants from b
pybind11::array_t<double> isochoric_invariants(const pybind11::array_t<double> &input, const double &J=0., const bool &copy=true);

//This function returns the isochoric principale stretches from b
pybind11::array_t<double> isochoric_pstretch(const pybind11::array_t<double> &input, const std::string &input_tensor="V", const double &J=0., const bool &copy=true);

} //namespace simpy