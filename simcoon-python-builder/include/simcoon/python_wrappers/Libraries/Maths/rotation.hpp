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
    
//This function returns the 3*3 rotation matrix according to an angle, an axis and depending if it is active or passive rotation
pybind11::array_t<double> fillR_angle(const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
    
//This function returns the 3*3 rotation matrix according to the three Euler angles, depending if it is active or passive rotation and the Euler convention (ex :"zxz")
pybind11::array_t<double> fillR_euler(const double &psi, const double &theta, const double &phi, const bool &active=true, const std::string &conv="zxz", const bool &copy=true);
    
//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from an angle and an axis
pybind11::array_t<double> fillQS_angle(const double &angle, const int &axis, const bool &active=true, const bool &copy=true);

//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress' from a rotation matrix
pybind11::array_t<double> fillQS_R(const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from an angle and an axis
pybind11::array_t<double> fillQE_angle(const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
    
//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain' from a rotation matrix
pybind11::array_t<double> fillQE_R(const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);

//These functions rotates a 6*6 stiffness arma::matrix (L)
pybind11::array_t<double> rotateL_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotateL_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_L(const pybind11::array_t<double> &input, const double &, const double &, const double &, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_L(const pybind11::array_t<double> &input, const double &, const double &, const double &, const bool &copy=true);
    
//These functions rotates a 6*6 compliance arma::matrix (M)
pybind11::array_t<double> rotateM_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotateM_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_M(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_M(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
    
//These functions rotates a 6*6 strain concentration (A)
pybind11::array_t<double> rotateA_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotateA_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_A(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_A(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
    
//These functions rotates a 6*6 stress concentration (B)
pybind11::array_t<double> rotateB_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotateB_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_B(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_B(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
    
//These functions rotates strain arma::vectors
pybind11::array_t<double> rotate_strain_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_strain_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_strain(const pybind11::array_t<double> &, const double &psi, const double &theta, const double &phi, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_strain(const pybind11::array_t<double> &, const double &psi, const double &theta, const double &phi, const bool &copy=true);

//These functions rotates stress arma::vectors
pybind11::array_t<double> rotate_stress_angle(const pybind11::array_t<double> &input, const double &angle, const int &axis, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_stress_R(const pybind11::array_t<double> &input, const pybind11::array_t<double> &R, const bool &active=true, const bool &copy=true);
pybind11::array_t<double> rotate_l2g_stress(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
pybind11::array_t<double> rotate_g2l_stress(const pybind11::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy=true);
    
} //namespace simpy