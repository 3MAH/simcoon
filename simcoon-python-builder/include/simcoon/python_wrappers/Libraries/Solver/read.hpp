#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function reads material properties to prepare a simulation
    void read_matprops(int &, boost::python::numpy::ndarray &, int &, double &, double &, double &, const boost::python::str &, const boost::python::str &);
    
} //namespace simpy