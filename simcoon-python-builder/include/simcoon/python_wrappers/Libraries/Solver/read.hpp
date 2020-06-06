#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function reads material properties to prepare a simulation
    void read_matprops(unsigned int &, boost::python::numpy::ndarray &, unsigned int &, double &, double &, double &, const std::string &, const std::string &);
    
} //namespace simpy
