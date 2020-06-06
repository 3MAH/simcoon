#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
boost::python::numpy::ndarray get_densities_ODF(const boost::python::numpy::ndarray &, const std::string &, const std::string &, const bool &);
    
void ODF_discretization(const int &, const int &, const double &, const double &, const std::string &, const boost::python::numpy::ndarray &, const std::string &, const std::string &, const std::string &, const std::string &, const int &);
    
} //namespace simpy
