#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
pybind11::tuple logarithmic(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);
pybind11::tuple logarithmic_R(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);
pybind11::tuple objective_rate(const std::string &corate_name, const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &return_de, const unsigned int &n_threads);
pybind11::array_t<double> Delta_log_strain(const pybind11::array_t<double> &D, const pybind11::array_t<double> &Omega, const double &DTime, const bool &copy=true);
pybind11::array_t<double> Lt_convert(const pybind11::array_t<double> &Lt, const pybind11::array_t<double> &F, const pybind11::array_t<double> &stress, const std::string &converter_key);

} //namespace simpy
