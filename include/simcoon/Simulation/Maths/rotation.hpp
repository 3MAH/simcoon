/* This file is part of simcoon.
 
 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file rotation.hpp
///@brief rotation of a Voigt tensor
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file rotation.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */

/**
 * @brief Rotates a 3D vector using a rotation matrix.
 * @param v The input 3D vector (arma::vec)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @return The rotated vector (arma::vec)
 * @details Example:
 * @code
 *      vec v = {1., 0., 0.};
 *      mat R = fillR(M_PI/4., 3);
 *      vec v_rot = rotate_vec(v, R);
 * @endcode
 */
arma::vec rotate_vec(const arma::vec &v, const arma::mat &R);

/**
 * @brief Rotates a 3D vector by an angle around a given axis.
 * @param v The input 3D vector (arma::vec)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @return The rotated vector (arma::vec)
 * @details Example:
 * @code
 *      vec v = {1., 0., 0.};
 *      vec v_rot = rotate_vec(v, M_PI/4., 3); // Rotate 45° around z-axis
 * @endcode
 */
arma::vec rotate_vec(const arma::vec &v, const double &angle, const int &axis);

/**
 * @brief Rotates a 3x3 matrix using a rotation matrix.
 * @param m The input 3x3 matrix (arma::mat)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @return The rotated matrix (arma::mat)
 * @details The rotation is performed as \f$ \mathbf{m}' = \mathbf{R} \cdot \mathbf{m} \cdot \mathbf{R}^T \f$
 * @code
 *      mat m = eye(3,3);
 *      mat R = fillR(M_PI/4., 3);
 *      mat m_rot = rotate_mat(m, R);
 * @endcode
 */
arma::mat rotate_mat(const arma::mat &m, const arma::mat &R);

/**
 * @brief Rotates a 3x3 matrix by an angle around a given axis.
 * @param m The input 3x3 matrix (arma::mat)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @return The rotated matrix (arma::mat)
 */
arma::mat rotate_mat(const arma::mat &m, const double &angle, const int &axis);

/**
 * @brief Generates a 3x3 rotation matrix for rotation around a single axis.
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active (alibi) rotation; if false, passive (alias) rotation (bool)
 * @return The 3x3 rotation matrix (arma::mat)
 * @details Example:
 * @code
 *      mat R = fillR(M_PI/4., 3); // 45° rotation around z-axis
 * @endcode
 */
arma::mat fillR(const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Generates a 3x3 rotation matrix using Euler angles.
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @param active If true (default), active (alibi) rotation; if false, passive (alias) rotation (bool)
 * @param conv The Euler angle convention: "zxz" (default), "zyz", "xyz", etc. Use "user" for custom convention (std::string)
 * @return The 3x3 rotation matrix (arma::mat)
 * @details Example:
 * @code
 *      mat R = fillR(M_PI/4., M_PI/6., M_PI/3., true, "zxz");
 * @endcode
 */
arma::mat fillR(const double &psi, const double &theta, const double &phi, const bool &active = true, const std::string &conv = "zxz");
    
/**
 * @brief Generates a 6x6 rotation matrix for stress tensors in Voigt notation.
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The 6x6 rotation matrix for stress (arma::mat)
 */
arma::mat fillQS(const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Generates a 6x6 rotation matrix for stress tensors from a 3x3 rotation matrix.
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The 6x6 rotation matrix for stress (arma::mat)
 */
arma::mat fillQS(const arma::mat &R, const bool &active = true);
    
/**
 * @brief Generates a 6x6 rotation matrix for strain tensors in Voigt notation.
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The 6x6 rotation matrix for strain (arma::mat)
 */
arma::mat fillQE(const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Generates a 6x6 rotation matrix for strain tensors from a 3x3 rotation matrix.
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The 6x6 rotation matrix for strain (arma::mat)
 */
arma::mat fillQE(const arma::mat &R, const bool &active = true);

/**
 * @brief Rotates a 6x6 stiffness matrix.
 * @param L The 6x6 stiffness matrix in Voigt notation (arma::mat)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 stiffness matrix (arma::mat)
 * @details Example:
 * @code
 *      mat L = L_iso(E, nu, "EnuG");
 *      mat L_rot = rotateL(L, M_PI/4., 3);
 * @endcode
 */
arma::mat rotateL(const arma::mat &L, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a 6x6 stiffness matrix using a 3x3 rotation matrix.
 * @param L The 6x6 stiffness matrix in Voigt notation (arma::mat)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 stiffness matrix (arma::mat)
 */
arma::mat rotateL(const arma::mat &L, const arma::mat &R, const bool &active = true);

/**
 * @brief Rotates a 6x6 compliance matrix.
 * @param M The 6x6 compliance matrix in Voigt notation (arma::mat)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 compliance matrix (arma::mat)
 */
arma::mat rotateM(const arma::mat &M, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a 6x6 compliance matrix using a 3x3 rotation matrix.
 * @param M The 6x6 compliance matrix in Voigt notation (arma::mat)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 compliance matrix (arma::mat)
 */
arma::mat rotateM(const arma::mat &M, const arma::mat &R, const bool &active = true);
    
/**
 * @brief Rotates a 6x6 strain localization matrix A.
 * @param A The 6x6 strain localization matrix (arma::mat)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 strain localization matrix (arma::mat)
 */
arma::mat rotateA(const arma::mat &A, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a 6x6 strain localization matrix A using a 3x3 rotation matrix.
 * @param A The 6x6 strain localization matrix (arma::mat)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 strain localization matrix (arma::mat)
 */
arma::mat rotateA(const arma::mat &A, const arma::mat &R, const bool &active = true);

/**
 * @brief Rotates a 6x6 stress localization matrix B.
 * @param B The 6x6 stress localization matrix (arma::mat)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 stress localization matrix (arma::mat)
 */
arma::mat rotateB(const arma::mat &B, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a 6x6 stress localization matrix B using a 3x3 rotation matrix.
 * @param B The 6x6 stress localization matrix (arma::mat)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated 6x6 stress localization matrix (arma::mat)
 */
arma::mat rotateB(const arma::mat &B, const arma::mat &R, const bool &active = true);
    
/**
 * @brief Rotates a stress vector in Voigt notation.
 * @param sigma The 6-component stress vector (arma::vec)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated stress vector (arma::vec)
 * @details Example:
 * @code
 *      vec sigma = {100., 50., 0., 25., 0., 0.};
 *      vec sigma_rot = rotate_stress(sigma, M_PI/4., 3);
 * @endcode
 */
arma::vec rotate_stress(const arma::vec &sigma, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a stress vector using a 3x3 rotation matrix.
 * @param sigma The 6-component stress vector (arma::vec)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated stress vector (arma::vec)
 */
arma::vec rotate_stress(const arma::vec &sigma, const arma::mat &R, const bool &active = true);
    
/**
 * @brief Rotates a strain vector in Voigt notation.
 * @param epsilon The 6-component strain vector (arma::vec)
 * @param angle The rotation angle in radians (double)
 * @param axis The axis of rotation: 1=x, 2=y, 3=z (int)
 * @param active If true (default), active rotation (bool)
 * @return The rotated strain vector (arma::vec)
 */
arma::vec rotate_strain(const arma::vec &epsilon, const double &angle, const int &axis, const bool &active = true);

/**
 * @brief Rotates a strain vector using a 3x3 rotation matrix.
 * @param epsilon The 6-component strain vector (arma::vec)
 * @param R The 3x3 rotation matrix (arma::mat)
 * @param active If true (default), active rotation (bool)
 * @return The rotated strain vector (arma::vec)
 */
arma::vec rotate_strain(const arma::vec &epsilon, const arma::mat &R, const bool &active = true);

/**
 * @brief Rotates a strain tensor from local to global coordinates using Euler angles.
 * @param E The 6-component strain vector in local coordinates (arma::vec)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 3x3 strain tensor in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_strain(const arma::vec &E, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a strain tensor from global to local coordinates using Euler angles.
 * @param E The 6-component strain vector in global coordinates (arma::vec)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 3x3 strain tensor in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_strain(const arma::vec &E, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a stress tensor from local to global coordinates using Euler angles.
 * @param sigma The 6-component stress vector in local coordinates (arma::vec)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 3x3 stress tensor in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_stress(const arma::vec &sigma, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a stress tensor from global to local coordinates using Euler angles.
 * @param sigma The 6-component stress vector in global coordinates (arma::vec)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 3x3 stress tensor in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_stress(const arma::vec &sigma, const double &psi, const double &theta, const double &phi);
    
/**
 * @brief Rotates a stiffness matrix from local to global coordinates using Euler angles.
 * @param L The 6x6 stiffness matrix in local coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 stiffness matrix in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_L(const arma::mat &L, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a stiffness matrix from global to local coordinates using Euler angles.
 * @param L The 6x6 stiffness matrix in global coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 stiffness matrix in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_L(const arma::mat &L, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a strain localization matrix from local to global coordinates using Euler angles.
 * @param A The 6x6 strain localization matrix in local coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 strain localization matrix in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_A(const arma::mat &A, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a strain localization matrix from global to local coordinates using Euler angles.
 * @param A The 6x6 strain localization matrix in global coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 strain localization matrix in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_A(const arma::mat &A, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a stress localization matrix from local to global coordinates using Euler angles.
 * @param B The 6x6 stress localization matrix in local coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 stress localization matrix in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_B(const arma::mat &B, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a stress localization matrix from global to local coordinates using Euler angles.
 * @param B The 6x6 stress localization matrix in global coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 stress localization matrix in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_B(const arma::mat &B, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a compliance matrix from local to global coordinates using Euler angles.
 * @param M The 6x6 compliance matrix in local coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 compliance matrix in global coordinates (arma::mat)
 */
arma::mat rotate_l2g_M(const arma::mat &M, const double &psi, const double &theta, const double &phi);

/**
 * @brief Rotates a compliance matrix from global to local coordinates using Euler angles.
 * @param M The 6x6 compliance matrix in global coordinates (arma::mat)
 * @param psi First Euler angle in radians (double)
 * @param theta Second Euler angle in radians (double)
 * @param phi Third Euler angle in radians (double)
 * @return The 6x6 compliance matrix in local coordinates (arma::mat)
 */
arma::mat rotate_g2l_M(const arma::mat &M, const double &psi, const double &theta, const double &phi);
    
    

/** @} */ // end of maths group

} //namespace simcoon
