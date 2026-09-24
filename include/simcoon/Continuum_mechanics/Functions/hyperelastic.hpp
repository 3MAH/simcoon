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
 * @brief A set of functions for hyperelastic material models.
 */

#pragma once
#include <string>
#include <vector>
#include <armadillo>

namespace simcoon{

/** @addtogroup hyperelastic
 *  @{
 */

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
 * @brief Provides the isochoric strain invariants, from the left Cauchy-Green principal stretches \f$ \lambda^2_1, \lambda^2_2 and \lambda^2_3\f$ Note that principal stretches \f$ \lambda_1, \lambda_2 and \lambda_3\f$ are the ones of the Eulerian stretch tensor \f$ \mathbf{v} \f$.
 *
 * \f[
    \begin{align}
    \bar{I}_1 = \bar{\lambda}_1^2 + \bar{\lambda}_2^2 + \bar{\lambda}_3^2 \\
    \bar{I}_2 = \bar{\lambda}_1^{-2} + \bar{\lambda}_2^{-2} + \bar{\lambda}_3^{-2} \\
    \bar{I}_3 = \bar{\lambda}_1^2 \bar{\lambda}_2^2 \bar{\lambda}_3^2 = 1
    \end{align}
 * \f]
 * where \f$ \bar{\lambda}_i = J^{-1/3} \lambda_i \f$ is the i-th isochoric principal stretch from a principal decomposition of the isochoric part of \f$\mathbf{b}\f$.
 * 
 * @param lambda a column vector of dimension 3 that contains the three principal stretches \f$ \lambda_1 \f$, \f$ \lambda_2 \f$ and \f$ \lambda_3 \f$ of the Eulerian stretch tensor \f$ \mathbf{v} \f$.
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
arma::vec isochoric_invariants(const arma::vec &lambda, const double &mJ = 0.);

/**
 * @brief Provides the isochoric principal stretches \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$ , from the eulerian stretch tensor \f$ \mathbf{v} \f$.
 *
 *  \f$ \lambda_1, \lambda_2 \f$ and \f$ \lambda_3 \f$ are the principal stretches of the Eulerian stretch tensor \f$ \mathbf{v} \f$ and:
 * \f[
       \bar{\lambda}_i = J^{-1/3} \lambda_i
 * \f]
 * 
 * @param V 3x3 matrix representing the Eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return a column vector of dimension 3 that contains the three isochoric principal stretches
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat V = zeros(3,3);
 *      mat R = zeros(3,3); 
 *      VR_decomposition(V, R, F);
 *      double J = det(F);
 *      vec lambdas_bar = isochoric_pstretch_from_V(V,J);
 * @endcode
*/
arma::vec isochoric_pstretch_from_V(const arma::mat &V, const double &mJ = 0.);

/**
 * @brief Provides the isochoric principal stretches \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$ , from the from the left Cauchy-Green tensor \f$ \mathbf{b} \f$.
 *
 *  \f$ \lambda^2_1, \lambda^2_2 \f$ and \f$ \lambda^2_3 \f$ are the principal components of the left Cauchy-Green tensor \f$ \mathbf{b} \f$ and:
 * \f[
       \bar{\lambda}_i = J^{-1/3} \lambda_i
 * \f]
 * 
 * @param b 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$.
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return a column vector of dimension 3 that contains the three isochoric principal stretches
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas_bar = isochoric_pstretch_from_b(b,J);
 * @endcode
*/
arma::vec isochoric_pstretch_from_b(const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric principal stretches \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$ , from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param input 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param input_tensor a string ("b" for the left Cauchy-Green tensor or "V" for the eulerian stretch tensor ) representing the selected input
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return a column vector of dimension 3 that contains the three isochoric principal stretches
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas_bar = isochoric_pstretch(b, "b", J);
 * @endcode
*/
arma::vec isochoric_pstretch(const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief Principal stretches \f$ \lambda^2_1, \lambda^2_2 \f$ and \f$ \lambda^2_3\f$ and principal directions, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda a column vector of dimension 3 that will contain the three principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param input 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param input_tensor a string ("b" for the left Cauchy-Green tensor or "V" for the eulerian stretch tensor ) representing the selected input
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas;
 *      mat n;
 *      pstretch(lambdas, n, b, "b", J);
 * @endcode
*/
void pstretch(arma::vec &lambda, arma::mat &n_pvector, const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief Principal stretches \f$ \lambda^2_1, \lambda^2_2 \f$ and \f$ \lambda^2_3\f$ and principal directions, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda a column vector of dimension 3 that will contain the three principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonal projector corresponding to a principal vector 
 * @param input 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param input_tensor a string ("b" for the left Cauchy-Green tensor or "V" for the eulerian stretch tensor ) representing the selected input
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas;
 *      mat n;
 *      std::vector<mat> N;
 *      pstretch(lambdas, n, N, b, "b", J);
 * @endcode
*/
void pstretch(arma::vec &lambda, arma::mat &n_pvector, std::vector<arma::mat> &N_projectors, const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief isochoric principal stretches \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$ and principal directions, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda_bar  a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param input 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param input_tensor a string ("b" for the left Cauchy-Green tensor or "V" for the eulerian stretch tensor ) representing the selected input
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas_bar;
 *      mat n;
 *      isochoric_pstretch(lambdas_bar, n, "b", J);
 * @endcode
*/
void isochoric_pstretch(arma::vec &lambda_bar, arma::mat &n_pvector, const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief isochoric principal stretches,  principal directions \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 \f$ and \f$ \bar{\lambda}^2_3\f$, principal directions and principal orthogonal projectors, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda_bar a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonal projector corresponding to a principal vector
 * @param input 3x3 matrix representing the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * @param input_tensor a string ("b" for the left Cauchy-Green tensor or "V" for the eulerian stretch tensor ) representing the selected input
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec lambdas_bar;
 *      mat n;
 *      std::vector<mat> N(3);
 *      isochoric_pstretch(lambdas_bar, n, N, b, "b", J);
 * @endcode
*/
void isochoric_pstretch(arma::vec &lambda_bar, arma::mat &n_pvector, std::vector<arma::mat> &N_projectors, const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief Provides the coeficients \f$ beta_i \f$ for the computation of Kirchoff stress using isochoric principal stretch models
 * see (Connolly et al. Computational Mechanics (2019) 64:1273–1288 : https://doi.org/10.1007/s00466-019-01707-1) for more details
 * 
 * @param dWdlambda_bar a column vector of dimension 3 that contains the three derivatives of the strain energy with respect to the isochoric principal stretches
 * @param lambda_bar a column vector of dimension 3 that contains the isochoric principal stretches \f$ \bar{lambda}_i \f$
 * @return a column vector of dimension 3 that contains the three coefficients \f$ \beta_1, \beta_2, \beta_3 \f$
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      vec beta_coefs = beta_coefs(dWdlambda_bar, lambda_bar);
 * @endcode
*/
arma::vec beta_coefs(const arma::vec &dWdlambda_bar, const arma::vec &lambda_bar);

/**
 * @brief Provides the coeficients \f$ gamma_{ij} \f$ for the computation of Kirchoff stress using isochoric principal stretch models
 * see (Connolly et al. Computational Mechanics (2019) 64:1273–1288 : https://doi.org/10.1007/s00466-019-01707-1) for more details
 * and the Ph.D Thesis of V. Le Sault https://theses.hal.science/file/index/docid/542506/filename/Manuscrit_final.pdf
 * 
 * @param dWdlambda_bar a column vector of dimension 3 that contains the three derivatives of the strain energy with respect to the isochoric principal stretches
 * @param dW2dlambda_bar2 a matrix of dimension 3 that contains the nine (6 independant) second derivatives of the strain energy with respect to the isochoric principal stretches 
 * @param lambda_bar a column vector of dimension 3 that contains the isochoric principal stretches \f$ \bar{lambda}_i \f$
 * @return a 3x3 matrix that contains the nine coefficients \f$ \gamma_{ij} \f$
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      mat dW2dlambda_bar2 = ...
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      mat gamma_coefs = gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
 * @endcode
*/
arma::mat gamma_coefs(const arma::vec &dWdlambda_bar, const arma::mat &dW2dlambda_bar2, const arma::vec &lambda_bar);

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
    \begin{align}
    \mathbf{\tau}_{\textrm{iso}} = \sum_{i=1}^3 \beta_i \left(\underline{n}_i \otimes \underline{n}_i \right) \\
    \beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{\lambda}_i} \f$ is the derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \f$ \bar{\lambda}_i \f$.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energy with respect to the isochoric principal stretches
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      mat m_tau_iso = tau_iso_hyper_pstretch(dWdlambda_bar, b, J);
 * @endcode
*/
arma::mat tau_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::mat &b, const double &mJ=0.);

/**
 * @brief Provides the isochoric part of the Kirchoff stress tensor.
 * 
 * The isochoric part of the Kirchoff stress tensor is defined as:
 * \f[
    \begin{align}
    \mathbf{\tau}_{\textrm{iso}} = \sum_{i=1}^3 \beta_i \left(\underline{n}_i \otimes \underline{n}_i \right) \\
    \beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{\lambda}_i} \f$ is the derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \f$ \bar{\lambda}_i \f$.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energy with respect to the isochoric principal stretches
 * @param lambda_bar a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonal projector corresponding to a principal vector 
 * @return 3x3 matrix representing the isochoric part of the Kirchoff stress tensor.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      vec lambdas_bar;
 *      mat n;
 *      std::vector<double> N(3);
 *      isochoric_pstretch(lambdas_bar, n, N, "b", J);
 *      mat m_tau_iso = tau_iso_hyper_pstretch(dWdlambda_bar, lambda_bar, N_projectors);
 * @endcode
*/
arma::mat tau_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, arma::vec &lambda_bar, std::vector<arma::mat> &N_projectors);

/**
 * @brief Provides the isochoric part of the Kirchoff stress tensor.
 * 
 * The isochoric part of the Kirchoff stress tensor is defined as:
 * \f[
        \mathbf{\tau}_{\textrm{iso}} = 2. \frac{\partial \bar{W} }{\partial \bar{I}_1 } \textrm{dev} \bar{\mathbf{b}} 
        + 2 \frac{\partial \bar{W} }{\partial \bar{I}_2 } \left( \textrm{tr} \bar{\mathbf{b}} \textrm{dev} \bar{\mathbf{b}}
        - \textrm{dev} \bar{\mathbf{b}}^2 \right)
         * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{I}_1 } \f$ and \f$ \frac{\partial \bar{W} }{\partial \bar{I}_2 } \f$ are the derivatives of the isochoric strain energy, \f$\mathbf{b}\f$ is the
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
 *      mat m_tau_iso = tau_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, b, J);
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
arma::mat tau_vol_hyper(const double &dUdJ, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric part of the Cauchy stress tensor.
 * 
 * The isochoric part of the Cauchy stress tensor is defined as:
 * \f[
    \begin{align}
    \mathbf{\sigma}_{\textrm{iso}} = \frac{1}{J} \sum_{i=1}^3 \beta_i \left(\underline{n}_i \otimes \underline{n}_i \right) \\
    \beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{\lambda}_i} \f$ is the derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \f$ \bar{\lambda}_i \f$.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energy with respect to the isochoric principal stretches
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 3x3 matrix representing the isochoric part of the Cauchy stress tensor.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      mat m_sigma_iso = sigma_iso_hyper_pstretch(dWdlambda_bar, b, J);
 * @endcode
*/
arma::mat sigma_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::mat &b, const double &mJ=0.);

/**
 * @brief Provides the isochoric part of the Cauchy stress tensor.
 * 
 * The isochoric part of the Cauchy stress tensor is defined as:
 * \f[
    \begin{align}
    \mathbf{\sigma}_{\textrm{iso}} = \frac{1}{J} \sum_{i=1}^3 \beta_i \left(\underline{n}_i \otimes \underline{n}_i \right) \\
    \beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{\lambda}_i} \f$ is the derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \f$ \bar{\lambda}_i \f$.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energy with respect to the isochoric principal stretches
 * @param lambda_bar  a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonal projector corresponding to a principal vector 
 * @param J the determinant of the transformation gradient \f$\mathbf{F}\f$ (mandatory since left Cauchy-Green deformation tensor \f$\mathbf{b}\f$ is not provided)
 * @return 3x3 matrix representing the isochoric part of the Cauchy stress tensor.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      vec lambdas_bar = isochoric_pstretch_from_b(b, J);
 *      mat n;
 *      std::vector<double> N(3);
 *      isochoric_pstretch(lambdas_bar, n, N, "b", J);
 *      mat m_tau_iso = sigma_iso_hyper_pstretch(dWdlambda_bar, lambda_bar, N_projectors, J);
 * @endcode
*/
arma::mat sigma_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::vec &lambda_bar, const std::vector<arma::mat> &N_projectors, const double &J);

/**
 * @brief Provides the isochoric part of the Cauchy stress tensor.
 * 
 * The isochoric part of the Kirchoff stress tensor is defined as:
 * \f[
        \mathbf{\tau}_{\textrm{iso}} = \frac{1}{J} \left[ 2. \frac{\partial \bar{W} }{\partial \bar{I}_1 } \textrm{dev} \bar{\mathbf{b}} 
        + 2 \frac{\partial \bar{W} }{\partial \bar{I}_2 } \left( \textrm{tr} \bar{\mathbf{b}} \textrm{dev} \bar{\mathbf{b}}
        - \textrm{dev} \bar{\mathbf{b}}^2 \right) \right]
         * \f]
 * where \f$ \frac{\partial \bar{W} }{\partial \bar{I}_1 } \f$ and \f$ \frac{\partial \bar{W} }{\partial \bar{I}_2 } \f$ are the derivatives of the isochoric strain energy, \f$\mathbf{b}\f$ is the
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
arma::mat sigma_vol_hyper(const double &dUdJ, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Provides the isochoric part of the hyperelastic tangent modulus, considering principal stretches
 * 
 * The isochoric part of the hyperelastic tangent modulus is defined as:
\f[ 
    \begin{align}    
        \mathbf{L}^t_{\textrm{iso}} &= \displaystyle \sum_{a,b = 1}^3 \left( \gamma_{ab} - \delta_{ab} \beta_a \right) 
               \left( \mathbf{n}_a \otimes \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \right) \\
               & + \displaystyle \sum_{a,b=1, \, a \neq b } \frac{\beta_b \lambda_a^2 - \beta_a \lambda_b^2}{\lambda_a^2 - \lambda_b^2}
               \left(\mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_a \otimes \mathbf{n}_b + \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \otimes \mathbf{n}_a \right)
    \end{align}
\f]
 * where \f$ \beta_{ij} \f$ and \f$ \gamma_{ij} \f$ depend on the derivatives of the isochoric strain energy with respect to principal stretches and \f$ \mathbf{n}_a \f$ is the a-th principal vector.
 *
 * @param dWdlambda_bar The derivative of the isochoric strain energy with respect to the isochoric principal stretches.
 * @param dW2dlambda_bar2 The second derivative of the isochoric strain energy with respect to the isochoric principal stretches. 
 * @param lambda_bar The isochoric principal stretches.  
 * @param n_pvectors The principal vectors, stacked as columns of a 3x3 matrix.  
 * @param J the determinant of the transformation gradient \f$\mathbf{F}\f$
 * @return 6x6 matrix representing the isochoric part of the hyperelastic tangent modulus.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      mat dW2dlambda_bar = {{ dW2dlambda_bar2_11, dW2dlambda_bar2_21, dW2dlambda_bar2_31}, { dW2dlambda_bar2_12, dW2dlambda_bar2_22, dW2dlambda_bar2_32}, { dW2dlambda_bar2_13, dW2dlambda_bar2_23, dW2dlambda_bar2_33}}; 
 *      vec lambdas_bar = isochoric_pstretch_from_b(b, J);
 *      mat n;
 *      std::vector<double> N(3);
 *      isochoric_pstretch(lambdas_bar, n_pvectors, N, "b", J);
 *      mat L_iso = L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, n_pvectors, J);
 * @endcode
*/
arma::mat L_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::mat &dW2dlambda_bar2, const arma::vec &lambda_bar, const arma::mat &n_pvectors, const double &J);

/**
 * @brief Provides the isochoric part of the hyperelastic tangent modulus, considering principal stretches
 * 
 * The isochoric part of the hyperelastic tangent modulus is defined as:
\f[ 
    \begin{align}    
        \mathbf{L}^t_{\textrm{iso}} &= \displaystyle \sum_{a,b = 1}^3 \left( \gamma_{ab} - \delta_{ab} \beta_a \right) 
               \left( \mathbf{n}_a \otimes \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \right) \\
               & + \displaystyle \sum_{a,b=1, \, a \neq b } \frac{\beta_b \lambda_a^2 - \beta_a \lambda_b^2}{\lambda_a^2 - \lambda_b^2}
               \left(\mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_a \otimes \mathbf{n}_b + \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \otimes \mathbf{n}_a \right)
    \end{align}
\f]
 * where \f$ \beta_{ij} \f$ and \f$ \gamma_{ij} \f$ depend on the derivatives of the isochoric strain energy with respect to principal stretches and \f$ \mathbf{n}_a \f$ is the a-th principal vector.
 *
 * @param dWdlambda_bar The derivative of the isochoric strain energy with respect to the isochoric principal stretches.
 * @param dW2dlambda_bar2 The second derivative of the isochoric strain energy with respect to the isochoric principal stretches. 
 * @param b 3x3 matrix representing the left Cauchy-Green deformation tensor \f$\mathbf{b}\f$ 
 * @param mJ the determinant of the transformation gradient \f$\mathbf{F}\f$ (optional)
 * @return 6x6 matrix representing the isochoric part of the hyperelastic tangent modulus.
 * 
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat b = L_Cauchy_Green(F); 
 *      double J = det(F);
 *      double dWdlambda_bar_1;
 *      double dWdlambda_bar_2;
 *      double dWdlambda_bar_3; 
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3};
 *      mat dW2dlambda_bar = {{ dW2dlambda_bar2_11, dW2dlambda_bar2_21, dW2dlambda_bar2_31}, { dW2dlambda_bar2_12, dW2dlambda_bar2_22, dW2dlambda_bar2_32}, { dW2dlambda_bar2_13, dW2dlambda_bar2_23, dW2dlambda_bar2_33}}; 
 *      mat L_iso = L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, b, J);
 * @endcode
*/
arma::mat L_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::mat &dW2dlambda_bar2, const arma::mat &b, const double &mJ);

/**
 * @brief Provides the isochoric part of the hyperelastic tangent modulus, considering Invariants
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
 * @param dWdI_1_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dWdI_2_bar The derivative of the isochoric strain energy with respect to the second isochoric invariant.
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
arma::mat L_vol_hyper(const double &dUdJ, const double &dU2dJ2, const arma::mat &b, const double &mJ = 0.);

/**
 * @brief Derivatives of an isochoric-invariant hyperelastic potential.
 *
 * The scalars every potential of the form
 * \f$ W(\bar{I}_1, \bar{I}_2) + U(J) \f$ hands to the stress and tangent
 * builders. Zero-initialised, so a potential only writes the terms it has.
 *
 * An anisotropic potential additionally writes the two vector members, one entry
 * per fibre family, holding the derivatives with respect to the fibre
 * pseudo-invariant \f$ \bar{I}^{*}_{4,i} \f$ (see structure_tensors_push_forward).
 * They stay empty for an isotropic potential, which is how the builders tell the
 * two cases apart.
 */
struct hyper_invariants_dW {
    double dWdI_1_bar = 0.;    ///< \f$ \partial W / \partial \bar{I}_1 \f$
    double dWdI_2_bar = 0.;    ///< \f$ \partial W / \partial \bar{I}_2 \f$
    double dW2dI_11_bar = 0.;  ///< \f$ \partial^2 W / \partial \bar{I}_1^2 \f$
    double dW2dI_12_bar = 0.;  ///< \f$ \partial^2 W / \partial \bar{I}_1 \partial \bar{I}_2 \f$
    double dW2dI_22_bar = 0.;  ///< \f$ \partial^2 W / \partial \bar{I}_2^2 \f$
    double dUdJ = 0.;          ///< \f$ \partial U / \partial J \f$
    double dU2dJ2 = 0.;        ///< \f$ \partial^2 U / \partial J^2 \f$
    arma::vec dWdI_a_bar;      ///< \f$ \partial W / \partial \bar{I}^{*}_{4,i} \f$, one per fibre family
    arma::vec dW2dI_aa_bar;    ///< \f$ \partial^2 W / \partial \bar{I}^{*\,2}_{4,i} \f$, one per fibre family
};

/**
 * @brief Isochoric-invariant potentials of hyper_potential_derivatives.
 *
 * Shared by the standalone UMAT (umat_generic_hyper_invariants, which maps its
 * 5-letter names onto these values) and by the modular composition, which
 * stores the value in its props: adding a potential here serves both.
 */
enum class HyperPotential {
    NEOHC = 0,  ///< compressible neo-Hookean, props [mu, kappa]
    MOORI = 1,  ///< Mooney-Rivlin, props [C10, C01, kappa]
    YEOHH = 2,  ///< Yeoh, props [C10, C20, C30, kappa]
    ISHAH = 3,  ///< Isihara, props [C10, C20, C01, kappa]
    GETHH = 4,  ///< Gent-Thomas, props [c1, c2, kappa]
    SWANH = 5,  ///< Swanson, props [N, kappa, (A, B, alpha, beta) x N]
    HOLZA = 6   ///< Gasser-Ogden-Holzapfel, props [C10, k1, k2, kappa_d, n_fam, (a0x, a0y, a0z) x n_fam, kappa]
};
// Every potential above, and OGDEN, may carry ONE more prop after those listed: the
// volumetric potential (VolumetricPotential, 0 when absent).

/**
 * @brief Volumetric part \f$ U(J) \f$ of a hyperelastic potential.
 *
 * Both share the ground-state bulk modulus \f$ U''(1) = \kappa \f$.
 */
enum class VolumetricPotential {
    LOG_J = 0,      ///< \f$ U = \kappa (J \ln J - J + 1) \f$ (default)
    QUADRATIC = 1   ///< \f$ U = \frac{\kappa}{2} (J - 1)^2 \f$
};

/**
 * @brief The volumetric potential selected by an optional trailing prop.
 *
 * @param props the props of the law
 * @param n_used how many of them the isochoric potential consumed; props(n_used), when
 *        present, selects the VolumetricPotential (0 or 1)
 */
VolumetricPotential volumetric_potential_of(const arma::vec &props, const arma::uword n_used);

/**
 * @brief First and second derivatives of \f$ U(J) \f$.
 */
void volumetric_derivatives(const VolumetricPotential &vol, const double &kappa, const double &J, double &dUdJ, double &dU2dJ2);

/**
 * @brief The fibre anisotropy carried by a hyperelastic potential's props.
 *
 * @see hyper_potential_anisotropy, which extracts it, and
 *      structure_tensors_push_forward, which consumes it.
 */
struct hyper_anisotropy {
    arma::mat a0;        ///< 3 x n_fam, one UNIT reference fibre direction \f$ \mathbf{a}_{0,i} \f$ per column; empty for an isotropic potential
    double kappa_d = 0.; ///< the Gasser-Ogden-Holzapfel dispersion \f$ \kappa_d \in [0, 1/3] \f$
};

/**
 * @brief The fibre directions and dispersion an anisotropic potential declares in its props.
 *
 * Isotropic potentials return an empty @c a0 and \f$ \kappa_d = 0 \f$. This is the
 * single place that knows where the directions sit in a potential's props, so
 * neither the standalone UMAT nor the modular block duplicates the layout.
 *
 * The reference directions are read as three direction cosines per family and
 * normalised defensively: the Python API expresses them as a
 * @c simcoon.Rotation applied to \f$ \mathbf{e}_1 \f$, so no Euler triplet (and
 * hence no gimbal lock) sits anywhere on the path from the user to the kernel.
 * They are expressed in the LOCAL material frame; the solver's material
 * orientation places them globally, exactly as for ELIST/ELORT.
 *
 * @param potential the potential (see HyperPotential for its props)
 * @param props the potential's own parameters, starting at index 0
 * @return the reference fibre directions and the dispersion
 * @throw std::invalid_argument if \f$ \kappa_d \notin [0, 1/3] \f$, if the family
 *        count is not strictly positive, or if a direction has zero norm
 */
hyper_anisotropy hyper_potential_anisotropy(const HyperPotential &potential, const arma::vec &props);

/**
 * @brief The pushed-forward isochoric structure tensors of a dispersed fibre family.
 *
 * With \f$ \bar{\mathbf{a}}_i = J^{-1/3} \mathbf{F} \mathbf{a}_{0,i} \f$ the
 * isochoric push-forward of the i-th reference direction, the Gasser-Ogden-Holzapfel
 * generalised structure tensor becomes, in the spatial configuration,
 * \f[
    \mathbf{A}_i = \kappa_d \, \bar{\mathbf{b}} + (1 - 3 \kappa_d) \, \bar{\mathbf{a}}_i \otimes \bar{\mathbf{a}}_i
 * \f]
 * whose trace is the fibre pseudo-invariant the potential is written in:
 * \f[
    \bar{I}^{*}_{4,i} = \textrm{tr} \, \mathbf{A}_i = \kappa_d \, \bar{I}_1 + (1 - 3 \kappa_d) \, \bar{I}_{4,i},
    \qquad \bar{I}_{4,i} = \mathbf{a}_{0,i} \cdot \bar{\mathbf{C}} \, \mathbf{a}_{0,i}
 * \f]
 * \f$ \kappa_d = 0 \f$ gives perfectly aligned fibres (the Holzapfel-Gasser-Ogden
 * 2000 model), \f$ \kappa_d = 1/3 \f$ an isotropic distribution, for which
 * \f$ \mathbf{A}_i = \bar{\mathbf{b}}/3 \f$ and the fibre term degenerates to a
 * function of \f$ \bar{I}_1 \f$ alone.
 *
 * Returning \f$ \mathbf{A}_i \f$ rather than the bare \f$ \bar{\mathbf{a}}_i \f$
 * is what makes the anisotropic stress and tangent reuse the isotropic
 * \f$ \bar{I}_1 \f$ machinery: \f$ \bar{\mathbf{b}} \f$ and
 * \f$ \bar{\mathbf{a}}_i \otimes \bar{\mathbf{a}}_i \f$ share the convected rate
 * form \f$ \mathbf{l}\mathbf{A} + \mathbf{A}\mathbf{l}^T - \frac{2}{3}
 * \textrm{tr}(\mathbf{d}) \mathbf{A} \f$, hence so does their combination, hence
 * \f$ \dot{\bar{I}}^{*}_{4,i} = 2 \, \textrm{dev} \mathbf{A}_i : \mathbf{d} \f$ —
 * the very relation \f$ \bar{I}_1 \f$ satisfies with \f$ \bar{\mathbf{b}} \f$.
 *
 * @warning The directions are carried from the REFERENCE configuration and pushed
 *          forward by @p F. When @p F is an elastic stretch that differs from the
 *          total one -- a modular composition whose mechanism subtracts an inelastic
 *          strain (plasticity, viscoelasticity) -- the fibres should first be convected
 *          into the intermediate configuration, and are not: the inelastic strain does
 *          not reorient them. The composition remains well posed and converges, and the
 *          approximation is exact while the inelastic strain is small or leaves the
 *          fibre directions fixed; its error grows with how much that strain reorients
 *          them. Representing it exactly would need the convected directions as state,
 *          which the additive corotational kinematics of the modular UMAT cannot express
 *          (there is no plastic deformation gradient to convect with). Damage is
 *          unaffected: it contributes no inelastic strain, so the elastic stretch is the
 *          total one and the push-forward is exact.
 *
 * @note When @p F is \f$ \mathbf{V}^{el} \f$ (the modular block), the push-forward is the
 *       exact corotated one, \f$ \mathbf{U}\mathbf{a}_0 = \mathbf{R}^T\mathbf{F}\mathbf{a}_0 \f$,
 *       ONLY because MODUL under finite strain rejects any corate other than 3 (log_R) --
 *       see select_umat_M_finite in umat_smart.cpp. Under another corate the same code would
 *       silently convect the fibres with a different spin, which is still objective but is a
 *       different model. That guard and this function must stay in step.
 *
 * @param F deformation gradient \f$ \mathbf{F} \f$ (or \f$ \mathbf{V}^{el} \f$ for an elastic state, as in hyper_invariants_response)
 * @param a0 3 x n_fam matrix of unit reference directions, one per column (empty gives an empty result)
 * @param kappa_d the dispersion \f$ \kappa_d \f$
 * @param mJ the determinant of \f$ \mathbf{F} \f$ (optional)
 * @return one 3x3 symmetric \f$ \mathbf{A}_i \f$ per fibre family
 *
 * @details Example:
 * @code
 *      mat F = randu(3,3);
 *      mat a0 = {{1.},{0.},{0.}};
 *      std::vector<mat> A = structure_tensors_push_forward(F, a0, 0.1, det(F));
 *      double I4_star = trace(A[0]);
 * @endcode
*/
std::vector<arma::mat> structure_tensors_push_forward(const arma::mat &F, const arma::mat &a0, const double &kappa_d, const double &mJ = 0.);

/**
 * @brief Derivatives of an isochoric-invariant potential.
 *
 * @param potential the potential (see HyperPotential for its props)
 * @param props the potential's own parameters, starting at index 0
 * @param I_bar the isochoric invariants \f$ (\bar{I}_1, \bar{I}_2, \bar{I}_3) \f$, as
 *        isochoric_invariants returns them
 * @param J determinant of the deformation gradient
 * @param A the pushed-forward structure tensors (structure_tensors_push_forward), one per
 *        fibre family. An anisotropic potential reads its pseudo-invariants off them as
 *        \f$ \bar{I}^{*}_{4,i} = \textrm{tr}\,\mathbf{A}_i \f$, so the invariant and the
 *        tensor the tangent is built from can never come from different code. Empty for an
 *        isotropic potential, which is the default.
 * @return the derivatives of the potential
 *
 * @note The framework assumes the potential is ADDITIVELY SEPARABLE in \f$ \bar{I}_1 \f$,
 *       \f$ \bar{I}_2 \f$ and each \f$ \bar{I}^{*}_{4,i} \f$: hyper_invariants_dW carries no
 *       \f$ \partial^2 W / \partial \bar{I}_1 \partial \bar{I}^{*}_4 \f$ slot, and
 *       hyper_invariants_response builds no cross term. A future coupled potential (the
 *       Holzapfel-Ogden myocardium model, for instance) needs that slot added, not just a
 *       new case here.
 */
hyper_invariants_dW hyper_potential_derivatives(const HyperPotential &potential, const arma::vec &props, const arma::vec &I_bar, const double &J, const std::vector<arma::mat> &A = {});

/**
 * @brief Cauchy stress and canonical box tangent of an invariant potential.
 *
 * Assembles \f$ \boldsymbol{\sigma} \f$ and \f$ \partial \hat{\boldsymbol{\tau}} /
 * \partial \mathbf{D}_e \f$ (Kirchhoff, no J, XBM rate — the same object the
 * small-strain boxes return) from the potential derivatives and the left
 * Cauchy-Green tensor.
 *
 * @p F is only used to move the tangent into the box convention. The standalone
 * UMAT passes the deformation gradient; a caller that has an ELASTIC state
 * rather than a total one passes \f$ \mathbf{V}^{el} = \exp(\boldsymbol{
 * \varepsilon}^{el}) \f$, whose square is @p b — the tangent is then
 * \f$ \partial \boldsymbol{\tau} / \partial \boldsymbol{\varepsilon}^{el} \f$.
 *
 * @param[in] dW potential derivatives at (@p b, @p J)
 * @param[in] b left Cauchy-Green tensor \f$ \mathbf{b} = \mathbf{F}\mathbf{F}^T \f$
 * @param[in] J \f$ \det \mathbf{F} \f$
 * @param[in] F deformation gradient (or V for an elastic state, see above)
 * @param[out] sigma Cauchy stress, 6-Voigt. Cauchy and not Kirchhoff because
 *             that is what the builders produce: a caller wanting
 *             \f$ \boldsymbol{\tau} \f$ multiplies by @p J once, rather than
 *             this function multiplying and the caller dividing back (which
 *             is not exact in floating point).
 * @param[out] Lt_box canonical box tangent, 6x6
 * @param[in] A the pushed-forward structure tensors of an anisotropic potential
 *            (structure_tensors_push_forward), one per fibre family and in the same
 *            order as @c dW.dWdI_a_bar. Empty for an isotropic potential, which is
 *            the default.
 */
void hyper_invariants_response(const hyper_invariants_dW &dW, const arma::mat &b, const double &J, const arma::mat &F, arma::vec &sigma, arma::mat &Lt_box, const std::vector<arma::mat> &A = {});

/** @} */ // end of hyperelastic group

} //namespace simcoon
