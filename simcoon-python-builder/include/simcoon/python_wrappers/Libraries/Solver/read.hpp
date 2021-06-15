#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function reads material properties to prepare a simulation
boost::python::tuple read_matprops(const std::string &, const std::string &);

boost::python::tuple read_path(const std::string &, const std::string &);
    
} //namespace simpy
