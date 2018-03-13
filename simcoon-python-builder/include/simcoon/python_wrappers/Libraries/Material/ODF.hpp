#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
boost::python::numpy::ndarray get_densities_ODF(const boost::python::numpy::ndarray &, const boost::python::str &, const boost::python::str &, const bool &);
    
void ODF_discretization(const int &, const int &, const double &, const double &, const boost::python::str &, const boost::python::numpy::ndarray &, const boost::python::str &, const boost::python::str &, const boost::python::str &, const boost::python::str &, const int &);
    
} //namespace simpy