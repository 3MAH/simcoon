#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function transforms the strain Voigt vector into a 3*3 strain matrix
pybind11::array_t<double> v2t_strain(const pybind11::array_t<double> &input, const bool & copy=true);

//This function transforms a 3*3 strain matrix into a strain Voigt vector
pybind11::array_t<double> t2v_strain (const pybind11::array_t<double> &input, const bool & copy=true);

//This function transforms the stress Voigt vector into a 3*3 stress matrix
pybind11::array_t<double> v2t_stress(const pybind11::array_t<double> &input, const bool & copy=true);

//This function transforms a 3*3 stress matrix into a stress Voigt vector
pybind11::array_t<double> t2v_stress (const pybind11::array_t<double> &input, const bool & copy=true);

} //namespace simpy