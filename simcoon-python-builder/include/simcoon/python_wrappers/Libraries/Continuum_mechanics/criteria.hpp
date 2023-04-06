#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function returns the Prager equivalent stress.
double Prager_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props);

//This function returns the derivative of the Prager equivalent stress.
pybind11::array_t<double> dPrager_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props, const bool &copy=true);

//This function returns the Prager equivalent stress.
double Tresca_stress(const pybind11::array_t<double> &input);

//This function returns the derivative of the Prager equivalent stress.
pybind11::array_t<double> dTresca_stress(const pybind11::array_t<double> &input, const bool &copy=true);

//Returns an anisotropic configurational tensor P in the Voigt format (6x6 numpy array), given its vector representation
pybind11::array_t<double> P_ani(const pybind11::array_t<double> &props, const bool &copy=true);

//Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion in the Voigt format (6x6 numpy array), given its vector representation
pybind11::array_t<double> P_hill(const pybind11::array_t<double> &props, const bool &copy=true);

//This function returns the Hill equivalent stress.
double Hill_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props);

//This function returns the derivative of the Hill equivalent stress.
pybind11::array_t<double> dHill_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props, const bool &copy=true);
    
//This function returns the Ani equivalent stress.
double Ani_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props);

//This function returns the derivative of the Ani equivalent stress.
pybind11::array_t<double> dAni_stress(const pybind11::array_t<double> &input, const pybind11::array_t<double> &props, const bool &copy=true);
    
//This function computes the selected equivalent stress function
double Eq_stress(const pybind11::array_t<double> &input, const std::string &criteria, const pybind11::array_t<double> &props);

//This function computes the deriavtive of the selected equivalent stress function
pybind11::array_t<double> dEq_stress(const pybind11::array_t<double> &input, const std::string &criteria, const pybind11::array_t<double> &props, const bool &copy=true);
    
} //namespace simpy