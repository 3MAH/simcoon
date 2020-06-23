#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
boost::python::tuple logarithmic(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);
    
boost::python::numpy::ndarray Delta_log_strain(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);

} //namespace simpy
