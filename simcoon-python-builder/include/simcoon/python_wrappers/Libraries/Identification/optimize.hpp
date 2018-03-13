#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
    double cost_solver(const boost::python::numpy::ndarray &);
    
} //namespace simpy