#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{

//This function returns a rotated vector (3) according to a rotation matrix
boost::python::numpy::ndarray rotate_vec_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);

//This function returns a rotated vector (3) according to an angle and an axis
boost::python::numpy::ndarray rotate_vec_angle(const boost::python::numpy::ndarray &, const double &, const int &);

//This function returns a rotated vector (3) according to a rotation matrix
boost::python::numpy::ndarray rotate_mat_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &);

//This function returns a rotated vector (3) according to an angle and an axis
boost::python::numpy::ndarray rotate_mat_angle(const boost::python::numpy::ndarray &, const double &, const int &);
    
//This function returns the 3*3 rotation matrix according to an angle, an axis and depending if it is active or passive rotation
boost::python::numpy::ndarray fillR_angle(const double &alpha, const int &axis, const bool &active);
    
//This function returns the 3*3 rotation matrix according to the three Euler angles, depending if it is active or passive rotation and the Euler convention (ex :"zxz")
boost::python::numpy::ndarray fillR_euler(const double &, const double &, const double &, const bool &, const std::string &);
    
//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress'
boost::python::numpy::ndarray fillQS_angle(const double &, const int &, const bool &);

//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'stress'
boost::python::numpy::ndarray fillQS_R(const boost::python::numpy::ndarray &, const bool &);

//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain'
boost::python::numpy::ndarray fillQE_angle(const double &, const int &, const bool &);
    
//This function returns the 6*6 rotation arma::matrix of a arma::vector of type 'strain'
boost::python::numpy::ndarray fillQE_R(const boost::python::numpy::ndarray &, const bool &);

//These functions rotates a 6*6 stiffness arma::matrix (L)
boost::python::numpy::ndarray rotateL_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotateL_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);

boost::python::numpy::ndarray rotate_l2g_L(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_L(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
    
//These functions rotates a 6*6 compliance arma::matrix (M)
boost::python::numpy::ndarray rotateM_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotateM_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);
boost::python::numpy::ndarray rotate_l2g_M(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_M(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
    
//These functions rotates a 6*6 strain concentration (A)
boost::python::numpy::ndarray rotateA_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotateA_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);
boost::python::numpy::ndarray rotate_l2g_A(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_A(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
    
//These functions rotates a 6*6 stress concentration (B)
boost::python::numpy::ndarray rotateB_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotateB_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);
boost::python::numpy::ndarray rotate_l2g_B(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_B(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
    
//These functions rotates strain arma::vectors
boost::python::numpy::ndarray rotate_strain_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotate_strain_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);
boost::python::numpy::ndarray rotate_l2g_strain(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_strain(const boost::python::numpy::ndarray &, const double &, const double &, const double &);

//These functions rotates stress arma::vectors
boost::python::numpy::ndarray rotate_stress_angle(const boost::python::numpy::ndarray &, const double &, const int &, const bool &);
boost::python::numpy::ndarray rotate_stress_R(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const bool &);
boost::python::numpy::ndarray rotate_l2g_stress(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
boost::python::numpy::ndarray rotate_g2l_stress(const boost::python::numpy::ndarray &, const double &, const double &, const double &);
    
} //namespace simpy