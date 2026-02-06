#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//This function returns a rotated vector (3) according to a rotation matrix
pybind11::array_t<double> rotate_vec_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &copy=true);

//This function returns a rotated vector (3) according to an angle and an axis
pybind11::array_t<double> rotate_vec_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &copy=true);

//This function returns a rotated matrix (3x3) according to a rotation matrix
pybind11::array_t<double> rotate_mat_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &copy=true);

//This function returns a rotated matrix (3x3) according to an angle and an axis
pybind11::array_t<double> rotate_mat_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &copy=true);

//These functions rotate a 6*6 stiffness matrix
pybind11::array_t<double> rotate_stiffness_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_stiffness_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotate a 6*6 compliance matrix
pybind11::array_t<double> rotate_compliance_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_compliance_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotate a 6*6 strain concentration tensor (A)
pybind11::array_t<double> rotate_strain_concentration_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_strain_concentration_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotate a 6*6 stress concentration tensor (B)
pybind11::array_t<double> rotate_stress_concentration_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_stress_concentration_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotate strain vectors - Can be used with stack of arrays (vectorized)
pybind11::array_t<double> rotate_strain_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_strain_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotate stress vectors - Can be used with stack of arrays (vectorized)
pybind11::array_t<double> rotate_stress_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_stress_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

} //namespace simpy
