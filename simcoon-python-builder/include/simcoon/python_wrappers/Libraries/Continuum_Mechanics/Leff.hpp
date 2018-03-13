#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//return the elastic stiffness tensor of a composite material
    boost::python::numpy::ndarray L_eff(const boost::python::str &, const boost::python::numpy::ndarray &, const int &, const double &, const double &, const double &, const boost::python::str &);
    
} //namespace simpy