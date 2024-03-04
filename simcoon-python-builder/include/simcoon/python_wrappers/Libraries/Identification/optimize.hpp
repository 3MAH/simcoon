#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
    double cost_solver(const pybind11::array_t<double> &);
    
} //namespace simpy