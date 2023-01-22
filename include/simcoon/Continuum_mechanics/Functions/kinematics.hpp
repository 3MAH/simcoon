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

///@file kinematics.hpp
///@brief A set of function that allows various strain transformation (in Finite strains)
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/**
* @file constitutive.hpp
* @author Yves Chemisky 
* @section The constitutive library contains all the function helpful to write a constitutive relation of constitutive model
*/

/**
 * @brief Provides the transformation gradient, from the Green-Lagrange strain and the rotation.
 *
 * The transformation gradient \f$\mathbf{F}\f$ is related to the Green-Lagrange strain \f$\mathbf{E}\f$ and the rotation \f$\mathbf{R}\f$ by the following equations:
 * \f[
 *     \mathbf{F} = \mathbf{R} \cdot \mathbf{U} \quad \mathbf{E} = \frac{1}{2} \left( \sqrt{\mathbf{U}^2} - \mathbf{I} \right)
 * \f]
 * where \f$\mathbf{U}\f$ is the right stretch tensor 3x3 matrix and I is the 3x3 identity matrix.
 * 
 * @param E 3x3 matrix representing the Green-Lagrange strain tensor
 * @param R 3x3 matrix representing the rotation tensor
 * @return a 3x3 matrix representing the transformation gradient \f$\mathbf{F}\f$
 * 
 * @details Example: 
 * @code
 *      mat E = randu(3,3);
 *      mat R = eye(3,3);
 *      mat F = ER_to_F(E, R);
 * @endcode
*/
arma::mat ER_to_F(const arma::mat &E, const arma::mat &R);

/**
 * @brief Provides the transformation gradient, from the logarithmic strain and the rotation.
 *
 * The transformation gradient \f$\mathbf{F}\f$ is related to the logarithmic strain \f$\mathbf{e}\f$ and the rotation \f$\mathbf{R}\f$ by the following equations:
 * \f[
 *     \mathbf{F} = \mathbf{V} \cdot \mathbf{R} \quad \mathbf{e} = \textrm{ln} \mathbf{V}
 * \f]
 * 
 * where \f$\mathbf{V}\f$ the right stretch tensor 3x3 matrix
 *
 * @param e 3x3 matrix representing the logarithmic strain tensor
 * @param R 3x3 matrix representing the rotation tensor
 * @return a 3x3 matrix representing the transformation gradient
 * 
 * @details Example:  
 * @code
 *      mat e = randu(3,3);
 *      mat R = eye(3,3);
 *      mat F = eR_to_F(e, R);
 * @endcode
*/
arma::mat eR_to_F(const arma::mat &e, const arma::mat &R);

/**
 * @brief Provides the gradient of the displacement (Lagrangian) from the transformation gradient.
 *
 * The gradient of the displacement (Lagrangian) is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \nabla_X \mathbf{U} = \mathbf{F} - \mathbf{I}
 * \f]
 * where \f$\mathbf{I}\f$ is the 3x3 identity matrix.
 * 
 * @param F 3x3 matrix representing the transformation gradient
 * @return a 3x3 matrix representing the gradient of the displacement (Lagrangian) 
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat GradU = G_UdX(F);
 * @endcode
*/
arma::mat G_UdX(const arma::mat &F);

/**
 * @brief Provides the gradient of the displacement (Eulerian) from the transformation gradient.
 * 
 * The gradient of the displacement (Eulerian) is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \nabla_x \mathbf{U} =  \mathbf{I} - \mathbf{F}^{-1}
 * \f]
 * where \f$\mathbf{I}\f$ is the 3x3 identity matrix.
 *
 * @param F 3x3 matrix representing the transformation gradient
 * @return 3x3 matrix representing the gradient of the displacement (Eulerian)
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat gradU = G_Udx(F);
 * @endcode
*/
arma::mat G_Udx(const arma::mat &F);

/**
 * @brief Provides the Right Cauchy-Green tensor from the transformation gradient.
 *
 * The Right Cauchy-Green tensor \f$\mathbf{C}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{C} =  \mathbf{F}^T \cdot \mathbf{F}
 * \f]
 *
 * @param F 3x3 matrix representing the transformation gradient
 * @return 3x3 matrix representing the Right Cauchy-Green tensor
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat C = R_Cauchy_Green(F);
 * @endcode
*/
arma::mat R_Cauchy_Green(const arma::mat &F);

/**
 * @brief Provides the Left Cauchy-Green tensor from the transformation gradient.
 *
 * The Left Cauchy-Green tensor \f$\mathbf{B}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{B} =  \mathbf{F} \cdot \mathbf{F}^T
 * \f]
 *
 * @param F 3x3 matrix representing the transformation gradient
 * @return 3x3 matrix representing the Left Cauchy-Green tensor
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat B = L_Cauchy_Green(F);
 * @endcode
*/
arma::mat L_Cauchy_Green(const arma::mat &F);

/**
 * @brief Provides the RU decomposition of the transformation gradient.
 *
 * The RU decomposition of the transformation gradient \f$\mathbf{F}\f$ is related to the matrices \f$\mathbf{R}\f$ and \f$\mathbf{U}\f$ by the following equations:
 * \f[
 *     \mathbf{F} = \mathbf{R} \cdot \mathbf{U} \quad \mathbf{U} = \sqrt{\mathbf{F}^T \cdot \mathbf{F}} \quad \mathbf{R} = \mathbf{F} \cdot \mathbf{U}^{-1}
 * \f]
 *
 * @param R 3x3 matrix that is the rotation matrix \f$\mathbf{R}\f$
 * @param U 3x3 matrix that is the right stretch tensor \f$\mathbf{U}\f$
 * @param F 3x3 matrix representing the transformation gradient
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat R = zeros(3,3);
 *      mat U = zeros(3,3);
 *      RU_decomposition(R, U, F);
 * @endcode
*/
void RU_decomposition(arma::mat &R, arma::mat &U, const arma::mat &F);

/**
 * @brief Provides the VR decomposition of the transformation gradient.
 *
 * The VR decomposition of the transformation gradient \f$\mathbf{F}\f$ is related to the matrices \f$\mathbf{R}\f$ and \f$\mathbf{V}\f$ by the following equations:
 * \f[
 *     \mathbf{F} = \mathbf{V} \cdot \mathbf{R} \quad \mathbf{V} = \sqrt{\mathbf{F} \cdot \mathbf{F}^T} \quad \mathbf{R} = \mathbf{V}^{-1} \cdot \mathbf{F}
 * \f]
 *
 * @param R 3x3 matrix that is the rotation matrix \f$\mathbf{R}\f$
 * @param V 3x3 matrix that is the left stretch tensor \f$\mathbf{V}\f$
 * @param F 3x3 matrix representing the transformation gradient
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat R = zeros(3,3);
 *      mat V = zeros(3,3);
 *      VR_decomposition(R, V, F);
 * @endcode
*/
void VR_decomposition(arma::mat &R, arma::mat &V, const arma::mat &F);
    
/**
 * @brief Provides the invariants of a symmetric tensor.
 *
 * The invariants of the symmetric tensor \f$\mathbf{X}\f$ are defined as (note that this definition is not unique):
 * \f[
 *     \mathbf{I}_1 = \textrm{trace} \left( X \right) \quad \mathbf{I}_2 = \frac{1}{2} \left( \textrm{trace} \left( X \right)^2 - \textrm{trace} \left( X^2 \right) \right) \quad \mathbf{I}_3 = \textrm{det} \left( X \right)
 * \f]
 *
 * @param X 3x3 matrix representing the symmetric tensor
 * @return Vector containing the three invariants of the tensor
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat C = R_Cauchy_Green(F);
 *      vec I = Inv_X(F);
 * @endcode
*/
arma::vec Inv_X(const arma::mat &X);

/**
 * @brief Provides the Cauchy tensor from the transformation gradient.
 *
 * The Cauchy tensor \f$\mathbf{b}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{b} = \left( \mathbf{F} \cdot \mathbf{F}^T \right)^{-1}
 * \f]
 *
 * Note that this tensor \f$\mathbf{b}\f$ is the inverse of the Left Cauchy-Green tensor \f$\mathbf{B}\f$.
 * 
 * @param F 3x3 matrix representing the transformation gradient
 * @return 3x3 matrix representing the Cauchy tensor
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = Cauchy(F);
 * @endcode
*/
arma::mat Cauchy(const arma::mat &F);
    
/**
 * @brief Provides the Green-Lagrange tensor from the transformation gradient.
 *
 * The Green-Lagrange tensor \f$\mathbf{E}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{E} = \frac{1}{2} \left( \mathbf{F}^T \cdot \mathbf{F} - \mathbf{I} \right)
 * \f]
 *
 * @param F 3x3 matrix representing the transformation gradient \f$\mathbf{F}\f$
 * @return 3x3 matrix representing the Green-Lagrange tensor \f$\mathbf{E}\f$
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat E = Green_Lagrange(F);
 * @endcode
*/
arma::mat Green_Lagrange(const arma::mat &F);

/**
 * @brief Provides the Euler-Almansi tensor from the transformation gradient.
 *
 * The Euler-Almansi tensor \f$\mathbf{A}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{A} = \frac{1}{2} \left( \mathbf{I} - \left( \mathbf{F} \cdot \mathbf{F}^T \right)^T \right)
 * \f]
 *
 * @param F 3x3 matrix representing the transformation gradient \f$\mathbf{F}\f$
 * @return 3x3 matrix representing the Euler-Almansi tensor \f$\mathbf{A}\f$
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat A = Euler_Almansi(F);
 * @endcode
*/
arma::mat Euler_Almansi(const arma::mat &F);
    
/**
 * @brief Provides the logarithmic strain tensor from the transformation gradient.
 *
 * The logarithmic strain tensor \f$\mathbf{e}\f$ is related to the transformation gradient \f$\mathbf{F}\f$ by the following equation:
 * \f[
 *     \mathbf{e} = \textrm{ln} \left( V \right) = \frac{1}{2} \textrm{ln} \left( V^2 \right) = \frac{1}{2} \textrm{ln} \left( \mathbf{F} \cdot \mathbf{F}^T \right)
 * \f]
 *
 * Note that the log function is here the matrix logarithm of symmetric/hermitian positive definite matrix V^2. TSee the armadillo documentation for more details of the log matrix  implementation.
 * 
 * @param F 3x3 matrix representing the transformation gradient \f$\mathbf{F}\f$
 * @return 3x3 matrix representing the logarithmic strain tensor \f$\mathbf{e}\f$
 *
 * @details Example:
 * @code
 *      mat F = randu(3,3);
 *      mat e = Log_strain(F);
 * @endcode
*/
arma::mat Log_strain(const arma::mat &F);

/**
 * @brief Provides the approximation of the Eulerian velocity tensor from the transformation gradient at two different times.
 *
 * The Eulerian velocity tensor \f$\mathbf{L}\f$ is related to the transformation gradients \f$\mathbf{F_0}\f$ at time \f$t_0\f$, \f$\mathbf{F_1}\f$ at time \f$t_1\f$, and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 * \f[
 *     \mathbf{L} = \frac{1}{\Delta t} \left( \mathbf{F}_1 - \mathbf{F}_0 \right) \cdot \mathbf{F}_1^{-1}
 * \f]
 *
 * @param F0 3x3 matrix representing the transformation gradient \f$\mathbf{F_0}\f$ at time \f$t_0\f$
 * @param F1 3x3 matrix representing the transformation gradient \f$\mathbf{F_1}\f$ at time \f$t_1\f$
 * @param DTime double representing the time difference Delta \f$\Delta t = t_1 - t_0\f$
 * @return 3x3 matrix representing the approximation of the Eulerian velocity tensor
 *
 * @details Example:
 * @code
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      mat L = finite_L(F0, F1, DTime);
 * @endcode
*/
arma::mat finite_L(const arma::mat &F0, const arma::mat &F1, const double &DTime);

/**
 * @brief Provides the approximation of the Eulerian symmetric rate tensor from the transformation gradient at time \f$t_0\f$, the transformation gradient at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$.
 *
 * The Eulerian symmetric rate tensor \f$\mathbf{D}\f$ is related to the transformation gradient \f$\mathbf{F}_0\f$ at time \f$t_0\f$, the transformation gradient \f$\mathbf{F}_1\f$ at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 
  * The Eulerian antisymmetric spin tensor \f$\mathbf{W}\f$ is related to the transformation gradient \f$\mathbf{F}_0\f$ at time \f$t_0\f$, the transformation gradient \f$\mathbf{F}_1\f$ at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 * \f[
 *      \mathbf{D} = \frac{1}{2} \left( \mathbf{L} - \mathbf{L}^T \right)
 * \f]
 * 
 * where \f$\mathbf{L}\f$ is the Eulerian velocity tensor. \f$\mathbf{D}\f$ is commonly referred as the rate of deformation (this necessitates although a specific discussion).
 *
 * @param F0 3x3 matrix representing the transformation gradient \f$\mathbf{F_0}\f$ at time \f$t_0\f$
 * @param F1 3x3 matrix representing the transformation gradient \f$\mathbf{F_1}\f$ at time \f$t_1\f$
 * @param DTime time difference \f$\Delta t = t_1 - t_0\f$
 * @return 3x3 matrix representing the approximation of the Eulerian symmetric rate tensor \f$\mathbf{D}\f$
 *
 * 
 * @details Example:
 * @code
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      mat DTime = 0.1;
 *      mat D = finite_D(F0, F1, DTime);
 * @endcode
*/
arma::mat finite_D(const arma::mat &F0, const arma::mat &F1, const double &DTime);

/**
 * @brief Provides the approximation of the Eulerian antisymmetric spin tensor from the transformation gradient at time \f$t_0\f$, the transformation gradient at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$.
 *
 * The Eulerian antisymmetric spin tensor \f$\mathbf{W}\f$ is related to the transformation gradient \f$\mathbf{F}_0\f$ at time \f$t_0\f$, the transformation gradient \f$\mathbf{F}_1\f$ at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 * \f[
 *      \mathbf{W} = \frac{1}{2} \left( \mathbf{L} - \mathbf{L}^T \right)
 * \f]
 * where \f$\mathbf{L}\f$ is the Eulerian velocity tensor. \f$\mathbf{W}\f$ corresponds to the Jaumann corotationnal rate.
 *
 * @param F0 3x3 matrix representing the transformation gradient \f$\mathbf{F_0}\f$ at time \f$t_0\f$
 * @param F1 3x3 matrix representing the transformation gradient \f$\mathbf{F_1}\f$ at time \f$t_1\f$
 * @param DTime time difference \f$\Delta t = t_1 - t_0\f$
 * @return 3x3 matrix representing the approximation of the Eulerian antisymmetric spin tensor \f$\mathbf{W}\f$
 * 
 * @details Example:
 * @code
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      mat DTime = 0.1;
 *      mat W = finite_W(F0, F1, DTime);
 * @endcode
*/
arma::mat finite_W(const arma::mat &F0, const arma::mat &F1, const double &DTime);
    
//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the RU or VR polar decomposition of the deformation gradient F (F0,F1,DTime)

/**
 * @brief Provides the approximation of the Eulerian rigid-body rotation spin tensor from the transformation gradient at time \f$t_0\f$, the transformation gradient at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$.
 *
 * The Eulerian rigid-body rotation spin tensor \f$\mathbf{\Omega}\f$ is related to the transformation gradient \f$\mathbf{F}_0\f$ at time \f$t_0\f$, the transformation gradient \f$\mathbf{F}_1\f$ at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 * 
 * \f[
 *      \mathbf{\Omega} = \frac{1}{\Delta t} \left( \mathbf{R}_1 - \mathbf{R}_0 \right) \cdot \mathbf{R}_1^{T}
 * \f]
 * 
 * where \f$\mathbf{R}_0\f$ and \f$\mathbf{R}_1\f$ are the rotation matrices obtained from RU decompositions of the transformation gradients \f$\mathbf{F}_0\f$ and \f$\mathbf{F}_1\f$, respectively. This corresponds to the Green-Naghdi corotationnal rate.
 * 
 * @param F0 3x3 matrix representing the transformation gradient \f$\mathbf{F_0}\f$ at time \f$t_0\f$
 * @param F1 3x3 matrix representing the transformation gradient \f$\mathbf{F_1}\f$ at time \f$t_1\f$
 * @param DTime time difference \f$\Delta t = t_1 - t_0\f$
 * @return 3x3 matrix representing the approximation of the Eulerian rigid-body rotation spin tensor (Green-Naghdi corotational rate).
 * 
 * @details Example:
 * @code
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      mat DTime = 0.1;
 *      mat Omega = finite_Omega(F0, F1, DTime);
 * @endcode
*/
arma::mat finite_Omega(const arma::mat &F0, const arma::mat &F1, const double &DTime);
    
/**
 * @brief Provides the Hughes-Winget approximation of a increment of rotation or transformation from the spin/velocity at time \f$t_0\f$, the spin/velocity at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$.
 *
 * The Hughes-Winget approximation of a increment of rotation or transformation \f$\mathbf{\Delta Q}\f$ is related to the spin/velocity \f$\mathbf{\Omega}_0\f$ at time \f$t_0\f$, the spin/velocity \f$\mathbf{\Omega}_1\f$ at time \f$t_1\f$ and the time difference \f$\Delta t = t_1 - t_0\f$ by the following equation:
 * \f[
 *      \mathbf{\Delta Q} = \left( \mathbf{I} + \frac{1}{2} \Delta t \, \mathbf{\Omega}_0 \right) \cdot \left( \mathbf{I} - \frac{1}{2} \Delta t \, \mathbf{\Omega}_1 \right)^{-1}
 * \f]
 * 
 * where \f$\mathbf{I}\f$ is the 3x3 identity matrix.
 *
 * @param Omega0 spin/velocity at time \f$t_0\f$
 * @param Omega1 spin/velocity at time \f$t_1\f$
 * @param DTime time difference \f$\Delta t = t_1 - t_0\f$
 * @return 3x3 matrix representing the Hughes-Winget approximation of a increment of rotation or transformation
 * 
 * @details Example:
 * @code
 *      mat Omega0 = randu(3,3);
 *      mat Omega1 = randu(3,3);
 *      mat DTime = 0.1;
 *      mat DQ = finite_DQ(Omega0, Omega1, DTime);
 * @endcode
*/
arma::mat finite_DQ(const arma::mat &Omega0, const arma::mat &Omega1, const double &DTime);
    
} //namespace simcoon
