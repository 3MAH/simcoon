#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <simcoon/Simulation/Identification/parameters.hpp>

namespace simpy{

simcoon::parameters build_parameters_full(const int &, const double&, const double &, const std::string&, const int &, const pybind11::list &);

pybind11::list parameters_get_input_files(simcoon::parameters &);

void parameters_set_input_files(simcoon::parameters &, const pybind11::list &);
    
} //namespace simcoon
