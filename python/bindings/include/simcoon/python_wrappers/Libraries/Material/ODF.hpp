#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path
pybind11::array_t<double> get_densities_ODF(const pybind11::array_t<double> &x_py, const std::string &path_data_py, const std::string &peak_file_py, const bool &radian);
    
void ODF_discretization(const int &nphases_rve, const int &num_phase_disc, const double &angle_min, const double &angle_max, const std::string &umat_name_py, const pybind11::array_t<double> &props_py, const std::string &path_data_py, const std::string &peak_file_py, const std::string &rve_init_file_py, const std::string &rve_disc_file_py, const int &angle_mat);
    
} //namespace simpy
