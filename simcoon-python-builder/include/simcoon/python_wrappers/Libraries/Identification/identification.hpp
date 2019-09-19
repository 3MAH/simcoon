#pragma once
#include <boost/python.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
void identification(const boost::python::str &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const int &, const boost::python::str &, const boost::python::str &, const boost::python::str &, const boost::python::str &, const boost::python::str &);

boost::python::list read_constants_py(const int &, const int &);
    
void copy_constants_py(const boost::python::list &, const std::string &, const std::string &);
    
void apply_constants_py(const boost::python::list &, const std::string &);

boost::python::list read_parameters_py(const int &);

void copy_parameters_py(const boost::python::list &, const std::string &, const std::string &);

void apply_parameters_py(const boost::python::list &, const std::string &);

double calc_cost(const int &, const std::string &);
    
} //namespace simpy