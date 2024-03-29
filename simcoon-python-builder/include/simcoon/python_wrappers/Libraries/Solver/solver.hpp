#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
    void solver(const std::string &, const pybind11::array_t<double> &, const int &, const double &, const double &, const double &, const int &, const int &, const std::string &, const std::string &, const std::string &, const std::string &);    
} //namespace simpy
