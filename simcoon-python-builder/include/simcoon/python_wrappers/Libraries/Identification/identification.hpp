#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
void identification(const std::string &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &);

pybind11::list read_constants_py(const int &, const int &);
    
void copy_constants_py(const pybind11::list &, const std::string &, const std::string &);
    
void apply_constants_py(const pybind11::list &, const std::string &);

pybind11::list read_parameters_py(const int &);

void copy_parameters_py(const pybind11::list &, const std::string &, const std::string &);

void apply_parameters_py(const pybind11::list &, const std::string &);

double calc_cost(const int &, const pybind11::list &);
    
} //namespace simpy
