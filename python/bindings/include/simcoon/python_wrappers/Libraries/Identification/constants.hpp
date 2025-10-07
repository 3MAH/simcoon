
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <simcoon/Simulation/Identification/constants.hpp>

namespace simpy{

simcoon::constants build_constants_full(const int &, const double &, const pybind11::array_t<double> &, const std::string &, const int &, const pybind11::list &);
    
pybind11::array_t<double> constants_get_input_values(simcoon::constants &);

pybind11::list constants_get_input_files(simcoon::constants &);

void constants_set_input_values(simcoon::constants &, const pybind11::array_t<double> &);

void constants_set_input_files(simcoon::constants &, const pybind11::list &);
    
} //namespace simcoon
