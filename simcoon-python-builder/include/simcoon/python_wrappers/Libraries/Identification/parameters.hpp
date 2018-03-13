

///@file constants.hpp
///@brief Handle of input constants exposed in python
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <boost/python.hpp>

namespace simpy{

simcoon::parameters build_parameters_full(const int &, const double&, const double &, const std::string&, const int &, const boost::python::list &);

boost::python::list parameters_get_input_files(simcoon::parameters &);

void parameters_set_input_files(simcoon::parameters &, const boost::python::list &);
    
} //namespace simcoon
