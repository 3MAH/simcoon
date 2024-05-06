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

/**
* @file hyperelastic.hpp
* @author Yves Chemisky 
* @section A set of function that allows various strain transformation (in Finite strains)
*/

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @brief Provides the isochoric strain invariants, from the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$ .
 *
 * \f[
    \begin{align}
    \bar{I}_1 = \textrm{tr} \bar{\mathbf{b}} \\
    \bar{I}_2 = \frac{1}{2} \left( \left(\textrm{tr} \bar{\mathbf{b}} \right)^2 - \textrm{tr} \bar{\mathbf{b}}^2 \right) \\
    \bar{I}_3 = \textrm{det} \bar{\mathbf{b}} = 1
    \end{align}
 * \f]
 * 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return a column vector of dimension 3 that contains the three isochoric invariants
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec I_bar = isochoric_invariants(b,J);
 * @endcode
*/
arma::vec isochoric_invariants(const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric strain invariants, from the left Cauchy-Green principal stretches \f$ \lambda^2_1, \lambda^2_2 and lambda^2_3\f$ Note that principal stretches \f$ \lambda_1, \lambda_2 and lambda_3\f$ are the ones of the Eulerian stretch tensor \f$ \mathbf{v} \f$.
 *
 * \f[
    \begin{align}
    \bar{I}_1 = \bar{\lambda}_1^2 + \bar{\lambda}_2^2 + \bar{\lambda}_3^2 \\
    \bar{I}_2 = \bar{\lambda}_1^{-2} + \bar{\lambda}_2^{-2} + \bar{\lambda}_3^{-2} \\
    \bar{I}_3 = \bar{\lambda}_1^2 \bar{\lambda}_2^2 \bar{\lambda}_3^2 = 1
    \end{align}
 * \f]
 * where \f$ \bar{\lambda}^2_i = J^{-1/3} \lambda^2_i \f$ is the i-th isochoric principal stretch from a principal decomposition of the isochoric part of \f$\mathbf{b}\f$.
 * 
 * @param lambda a column vector of dimension 3 that contains the three principal stretches \lambda_1, \lambda_2 and lambda_3 of the Eulerian stretch tensor \f$ \mathbf{v} \f$..
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return a column vector of dimension 3 that contains the three isochoric invariants
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas = eigen_sym(sqrtmat(b));
 *      vec I_bar = isochoric_invariants(lambdas,J);
 * @endcode
*/
arma::vec isochoric_invariants(const arma::vec &barlambda, const double &mJ = 0.);

/**
 * @brief Provides the coeficients \f$ a_i \f$ for the computation of invariants-based Kirchoff stress tensor
 * see (Connolly et al. Computational Mechanics (2019) 64:1273–1288 : https://doi.org/10.1007/s00466-019-01707-1) for more details
 * 
 * @param dWdI_1_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the second isochoric invariant. 
 * @param I_bar a column vector of dimension 3 that contains the three isochoric invariants
 * @return a column vector of dimension 2 that contains the two coefficients \f$ a_1 and a_2 \f$
 * @details Example: 
 * @code
 *      double dWdI_1_bar;
 *      double dWdI_2_bar; 
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec I_bar = isochoric_invariants(b,J);
 *      vec a_coefs = a_coefs(dWdI_1_bar, dWdI_2_bar, I_bar);
 * @endcode
*/
arma::vec a_coefs(const double &dWdI_1_bar, const double &dWdI_2_bar, const arma::vec &I_bar);

/**
 * @brief Provides the coeficients \f$ b_i \f$ for the computation of hyperlastic tangent modulus
 * see (Connolly et al. Computational Mechanics (2019) 64:1273–1288 : https://doi.org/10.1007/s00466-019-01707-1) for more details
 * 
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dW2dI_11_bar The second derivative of the isochoric strain energy with respect to the first isochoric invariant. 
 * @param dW2dI_12_bar The second derivative of the isochoric strain energy with respect to the first and second isochoric invariant.  
 * @param dW2dI_22_bar The second derivative of the isochoric strain energy with respect to the second isochoric invariant.  
 * @param I_bar a column vector of dimension 3 that contains the three isochoric invariants
 * @return a column vector of dimension 4 that contains the four coefficients \f$ b_1, b_2, b_3 and b_4 \f$.
 * @details Example: 
 * @code
 *      double dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar; 
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec I_bar = isochoric_invariants(b,J);
 *      vec b_coefs = b_coefs(dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar, I_bar);
 * @endcode
*/
arma::vec b_coefs(const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const arma::vec &I_bar);

/**
 * @brief Provides the coeficients \f$ delta_i \f$ for the computation of hyperlastic tangent modulus
 * see (Connolly et al. Computational Mechanics (2019) 64:1273–1288 : https://doi.org/10.1007/s00466-019-01707-1) for more details
 * 
 * @param a_coefs a column vector of dimension 2 that contains the two coefficients \f$ a_1 and a_2 \f$
 * @param b_coefs  column vector of dimension 4 that contains the four coefficients \f$ b_1, b_2, b_3 and b_4 \f$
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @return a column vector of dimension 8 that contains the eight coefficients \f$ \delta_1, \dots, \delta_8 \f$
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      ...
 *      vec a_coefs = a_coefs(dWdI_1_bar, dWdI_2_bar, I_bar);
 *      vec b_coefs = b_coefs(dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar, I_bar); 
 *      vec delta_coefs = delta_coefs(a_coefs, b_coefs, b);
 * @endcode
*/
arma::vec delta_coefs(const arma::vec &a_coefs, const arma::vec &b_coefs, const arma::mat &b);

/**
 * @brief Provides the isochoric part of the Kirchoff stress tensor.
 * 
 * The isochoric part of the Kirchoff stress tensor is defined as:
 * \f[
        \mathbf{\tau}_{\textrm{iso}} = 2. \frac{\partial \bar{W} }{\partial \bar{I}_1 } \textrm{dev} \bar{\mathbf{b}} 
        + 2 \frac{\partial \bar{W} }{\partial \bar{I}_2 } \left( \textrm{tr} \bar{\mathbf{b}} \textrm{dev} \bar{\mathbf{b}}
        - \textrm{dev} \bar{\mathbf{b}}^2 \right)
         * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{I}_1 } \f$ and \f$ \frac{\partial \bar{W} }{\partial \bar{I}_2 } \f$ arethe derivatives of the isochoric strain energy, \f$\mathbf{b}\f$ is the
 * left Cauchy-Green deformation tensor and \f$\mathbf{I}\f$ is the 3x3 identity matrix.
 *
 * @param dWdI_1_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the second isochoric invariant. 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      double dWdI_1_bar;
 *      double dWdI_2_bar; 
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      mat m_tau_iso = tau_iso(dWdI_1_bar, dWdI_2_bar, b, J);
 * @endcode
*/
arma::mat tau_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the volumetric part of the Kirchoff stress tensor.
 * 
 * The volumetric part of the Kirchoff stress tensor is related to the derivative of the volumetric strain energy \f$ U \f$:
 * \f[
        \mathbf{\tau}_{\textrm{vol}} = J \frac{\partial U}{\partial J} \, \mathbf{I}
 * \f]
 *
 * @param dUdJ the derivative of the volumetric strain energy with respect to \f$ J \f$ 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$ 
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      double dUdJ;
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      mat m_tau_vol = tau_vol(dUdJ, J);
 * @endcode
*/
arma::mat tau_vol_hyper_invariants(const double &dUdJ, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric part of the Cauchy stress tensor.
 * 
 * The isochoric part of the Kirchoff stress tensor is defined as:
 * \f[
        \mathbf{\tau}_{\textrm{iso}} = \frac{1}{J} \left[ 2. \frac{\partial \bar{W} }{\partial \bar{I}_1 } \textrm{dev} \bar{\mathbf{b}} 
        + 2 \frac{\partial \bar{W} }{\partial \bar{I}_2 } \left( \textrm{tr} \bar{\mathbf{b}} \textrm{dev} \bar{\mathbf{b}}
        - \textrm{dev} \bar{\mathbf{b}}^2 \right) \right]
         * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{I}_1 } \f$ and \f$ \frac{\partial \bar{W} }{\partial \bar{I}_2 } \f$ arethe derivatives of the isochoric strain energy, \f$\mathbf{b}\f$ is the
 * left Cauchy-Green deformation tensor and \f$\mathbf{I}\f$ is the 3x3 identity matrix.
 *
 * @param dWdI_1_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the second isochoric invariant. 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      double dWdI_1_bar;
 *      double dWdI_2_bar; 
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      mat m_sigma_iso = sigma_iso(dWdI_1_bar, dWdI_2_bar, b, J);
 * @endcode
*/
arma::mat sigma_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the volumetric part of the Kirchoff stress tensor.
 * 
 * The volumetric part of the Cauchy stress tensor is related to the derivative of the volumetric strain energy \f$ U \f$:
 * \f[
        \mathbf{\sigma}_{\textrm{vol}} = \frac{\partial U}{\partial J} \, \mathbf{I}
 * \f]
 *
 * @param dUdJ the derivative of the volumetric strain energy with respect to \f$ J \f$ 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$ 
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      double dUdJ;
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      mat m_sigma_vol = sigma_vol(dUdJ, J);
 * @endcode
*/
arma::mat sigma_vol_hyper_invariants(const double &dUdJ, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric part of the hyperelastic tangent modulus
 * 
 * The isochoric part of the hyperelastic tangent modulus is defined as:
\f[ 
    \begin{align}    
        \mathbf{L}^t_{\textrm{iso}} &= \delta_1 \left( \bar{\mathbf{b}} \otimes \bar{\mathbf{b}} \right) + \delta_2 \left[ \left( \bar{\mathbf{b}} \otimes \bar{\mathbf{b}}^2 \right) + \left( \bar{\mathbf{b}}^2 \otimes \bar{\mathbf{b}} \right) \right] \\
                &+ \delta_3 \left[ \left( \bar{\mathbf{b}} \otimes \mathbf{I} \right) + \left(\mathbf{I} \otimes \bar{\mathbf{b}} \right) \right] + \delta_4 \left( \bar{\mathbf{b}}^2 \times \bar{\mathbf{b}}^2 \right) \\
                &+ \delta_5 \left[ \left( \bar{\mathbf{b}}^2 \otimes \mathbf{I} \right) + \left(\mathbf{I} \otimes \bar{\mathbf{b}}^2 \right) \right]  + \delta_6 \left( \mathbf{I} \otimes \mathbf{I} \right) \\
                &+ \delta_7 \left( \mathbf{I} \odot \mathbf{I} \right)  + \delta_8 \left( \bar{\mathbf{b}} \odot \bar{\mathbf{b}} \right) \\
    \end{align}
\f]
 * where \f$ delta_i \f$ depends on the derivatives of the isochoric strain energy, \f$\mathbf{b}\f$ is the
 * left Cauchy-Green deformation tensor and \f$ J \f$ is the determinant of the transformation gradient
 *
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dW2dI_11_bar The second derivative of the isochoric strain energy with respect to the first isochoric invariant. 
 * @param dW2dI_12_bar The second derivative of the isochoric strain energy with respect to the first and second isochoric invariant.  
 * @param dW2dI_22_bar The second derivative of the isochoric strain energy with respect to the second isochoric invariant.  
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 6x6 matrix representing the isochoric part of the hyperelastic tangent modulus.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdI_1_bar, dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar;
 *      mat L_iso = L_iso_hyper_invariants(delta_coefs, b, J);
 * @endcode
*/
arma::mat L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const arma::mat &b, const double &mJ=0.);

/**
 * @brief Provides the volumetric part of the hyperelastic tangent modulus
 * 
 * The volumetric part of the hyperelastic tangent modulus is defined as:
\f[ 
    \mathbf{L}^t_{\textrm{vol}} = J \left( \frac{\partial U}{\partial J} + \frac{\partial^2 U}{\partial J^2 \, J} \right) \left( \mathbf{I} \otimes \mathbf{I} \right) - 2 \frac{\partial U}{\partial J} \, J \left( \mathbf{I} \odot \mathbf{I} \right)
\f]
 * where U is the volumetric strain energy and \f$ J \f$ is the determinant of the transformation gradient
 *
 * @param dUdJ the derivative of the volumetric strain energy with respect to \f$ J \f$ 
 * @param dU2dJ2 the second derivative of the volumetric strain energy with respect to \f$ J \f$ 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 6x6 matrix representing the volumetric part of the hyperelastic tangent modulus.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dUdJ, dU2dJ2;
 *      mat L_vol = L_vol_hyper_invariants(dUdJ, dU2dJ2, J);
 * @endcode
*/
arma::mat L_vol_hyper_invariants(const double &dUdJ, const double &dU2dJ2, const arma::mat &b, const double &mJ = 0.);

} //namespace simcoon
