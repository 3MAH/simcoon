

///@file constants.hpp
///@brief Handle of input constants exposed in python
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

simcoon::constants build_constants_full(const int &, const double &, const boost::python::numpy::ndarray &, const std::string &, const int &, const boost::python::list &);
    
boost::python::numpy::ndarray constants_get_input_values(simcoon::constants &);

boost::python::list constants_get_input_files(simcoon::constants &);

void constants_set_input_values(simcoon::constants &, const boost::python::numpy::ndarray &);

void constants_set_input_files(simcoon::constants &, const boost::python::list &);
    
} //namespace simcoon
