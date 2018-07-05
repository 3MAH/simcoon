#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function transforms the strain Voigt vector into a 3*3 strain matrix
boost::python::numpy::ndarray v2t_strain(const boost::python::numpy::ndarray &);

//This function transforms a 3*3 strain matrix into a strain Voigt vector
boost::python::numpy::ndarray t2v_strain (const boost::python::numpy::ndarray &);

//This function transforms the stress Voigt vector into a 3*3 stress matrix
boost::python::numpy::ndarray v2t_stress(const boost::python::numpy::ndarray &);

//This function transforms a 3*3 stress matrix into a stress Voigt vector
boost::python::numpy::ndarray t2v_stress (const boost::python::numpy::ndarray &);

} //namespace simpy