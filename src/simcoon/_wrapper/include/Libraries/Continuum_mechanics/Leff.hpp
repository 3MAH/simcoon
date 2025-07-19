#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Return the elastic stiffness tensor of a composite material
pybind11::array_t<double> L_eff(const std::string &umat_name, const pybind11::array_t<double> &props, const int &nstatev, const double &psi_rve=0., const double &theta_rve=0., const double &phi_rve=0.);
    
} //namespace simpy
