#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
    void solver(const boost::python::str &, const boost::python::numpy::ndarray &, const int &, const double &, const double &, const double &, const int &, const boost::python::str &, const boost::python::str &, const boost::python::str &, const boost::python::str &);
 
//void solver();
    
} //namespace simpy