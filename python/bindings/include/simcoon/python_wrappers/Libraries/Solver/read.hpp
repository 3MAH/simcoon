#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function reads material properties to prepare a simulation
pybind11::tuple read_matprops(const std::string &, const std::string &);

pybind11::tuple read_path(const std::string &, const std::string &);
    
} //namespace simpy
