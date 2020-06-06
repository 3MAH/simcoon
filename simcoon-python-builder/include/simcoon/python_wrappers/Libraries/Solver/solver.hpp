#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
    void solver(const std::string &, const boost::python::numpy::ndarray &, const int &, const double &, const double &, const double &, const int &, const std::string &, const std::string &, const std::string &, const std::string &);
 
//void solver();
    
} //namespace simpy
