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
arma::vec isochoric_invariants(const arma::vec &lambda, const double &mJ = 0.);

/**
 * @brief Provides the isochoric principal stretches \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$ , from the eulerian stretch tensor \f$ \mathbf{v} \f$.
 *
 *  \f$ \lambda_1, \lambda_2 and lambda_3\f$ are the principal stretches of the Eulerian stretch tensor \f$ \mathbf{v} \f$ and:
 * \f[
       \f$ \bar{\lambda}_i = J^{-1/3} \lambda_i
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
 *  \f$ \lambda^2_1, \lambda^2_2 and lambda^2_3\f$ are the principal components of the left Cauchy-Green tensor \f$ \mathbf{v} \f$ and:
 * \f[
       \f$ \bar{\lambda}_i = J^{-1/3} \lambda_i
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
 * @brief Principal stretches \f$ \lambda^2_1, \lambda^2_2 and \lambda^2_3\f$ and principal directions, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
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
void pstretch(arma::vec &lambda, arma::mat &n_pvector, const arma::mat &input, const std::string &input_tensor, const double &mJ = 0.);

/**
 * @brief Principal stretches \f$ \lambda^2_1, \lambda^2_2 and \lambda^2_3\f$ and principal directions, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda_bar  a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonl projector corresponding to a principal vector 
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
 * @brief isochoric principal stretches,  principal directions \f$ \bar{\lambda}^2_1, \bar{\lambda}^2_2 and \bar{\lambda}^2_3\f$, principal directions and principal orthogonal projectors, from either the left Cauchy-Green tensor \f$ \mathbf{b} \f$ or the the eulerian stretch tensor \f$ \mathbf{v} \f$.
 * 
 * @param lambda _bar a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param n_pvector a 3x3 matrix, where each column is a principal direction vector.
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonl projector corresponding to a principal vector
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
 *      std::vector<double> N(3);
 *      isochoric_pstretch(lambdas_bar, n, N, "b", J);
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
    beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\bar{\lambda}_i} } \f$ is derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \bar{\lambda}_i.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energywith respect to the isochoric principal stretches
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
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3}
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      mat m_tau_iso = t_iso_hyper_pstretch(dWdlambda_bar, b, J);
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
    beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\bar{\lambda}_i} } \f$ is derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \bar{\lambda}_i.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energywith respect to the isochoric principal stretches
 * @param lambda_bar  a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonl projector corresponding to a principal vector 
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
 *      vec lambdas_bar;
 *      mat n;
 *      std::vector<double> N(3);
 *      isochoric_pstretch(lambdas_bar, n, N, "b", J);
 *      mat m_tau_iso = t_iso_hyper_pstretch(dWdlambda_bar, lambda_bar, N_projectors);
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
    beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\bar{\lambda}_i} } \f$ is derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \bar{\lambda}_i.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energywith respect to the isochoric principal stretches
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
 *      vec dWdlambda_bar = {dWdlambda_bar_1, dWdlambda_bar_2, dWdlambda_bar_3}
 *      lambda_bar = isochoric_pstretch_from_b(b, J);
 *      mat m_tau_iso = sigma_iso_hyper_pstretch(dWdlambda_bar, b, J);
 * @endcode
*/
arma::mat sigma_iso_hyper_pstretch(const arma::vec &dWdlambda_bar, const arma::mat &b, const double &mJ=0.);

/**
 * @brief Provides the isochoric part of the Cauchy stress tensor.
 * 
 * The isochoric part of the Cauchy stress tensor is defined as:
 * \f[
    \begin{align}
    \mathbf{\tau}_{\textrm{iso}} = \sum_{i=1}^3 \beta_i \left(\underline{n}_i \otimes \underline{n}_i \right) \\
    beta_i = \bar{\lambda}_i \frac{\partial W}{\partial \bar{\lambda}_i} - \frac{1}{3} \sum_{j=1}^3 \bar{\lambda}_j \frac{\partial W}{\partial \bar{\lambda}_j}
    \end{align}
 * \f]
 * where \f$ \frac{\partial \bar{W} }{\bar{\lambda}_i} } \f$ is derivative of the isochoric strain energy with respect to the i-th isochoric principal stretch \bar{\lambda}_i.
 *
 * @param dWdlambda_bar A column vector of dimension 3 that contains the derivatives of the isochoric strain energywith respect to the isochoric principal stretches
 * @param lambda_bar  a column vector of dimension 3 that will contain the three isochoric principal stretches
 * @param N_projectors a std::vector of 3x3 matrices, each one being an orthogonl projector corresponding to a principal vector 
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
        \mathbf{L}^t_{\textrm{iso}} &= \displaystyle \sum_{a,b = 1}^3 \letf( \gamma_{ab} - \delta_{ab} \beta_a \right) 
               \left( \mathbf{n}_a \otimes \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \right)
               & + \displaystyle \sum_{a,b=1 \, a \neq b } \frac{\beta_b \lambda_a^2 - \beta_a \lambda_b^2}{\lambda_a^2 - \lambda_b^2}
               \left(\mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_a \otimes \mathbf{n}_b + \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \otimes \mathbf{n}_a \right) \\
    \end{align}
\f]
 * where \f$ \beta_{ij} \f$ and \f$ \gamma_{ij} \f$ depends on the derivatives of the isochoric strain energy with respect to principal stretches and \f$ \mathbf{n}_a \f$ is the a-th principal vector
 *
 * @param dWdlambda_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dW2dlambda_bar2 The second derivative of the isochoric strain energy with respect to the first isochoric invariant. 
 * @param lambda_bar The second derivative of the isochoric strain energy with respect to the first and second isochoric invariant.  
 * @param N_projectors The second derivative of the isochoric strain energy with respect to the second isochoric invariant.  
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
        \mathbf{L}^t_{\textrm{iso}} &= \displaystyle \sum_{a,b = 1}^3 \letf( \gamma_{ab} - \delta_{ab} \beta_a \right) 
               \left( \mathbf{n}_a \otimes \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \right)
               & + \displaystyle \sum_{a,b=1 \, a \neq b } \frac{\beta_b \lambda_a^2 - \beta_a \lambda_b^2}{\lambda_a^2 - \lambda_b^2}
               \left(\mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_a \otimes \mathbf{n}_b + \mathbf{n}_a \otimes \mathbf{n}_b \otimes \mathbf{n}_b \otimes \mathbf{n}_a \right) \\
    \end{align}
\f]
 * where \f$ \beta_{ij} \f$ and \f$ \gamma_{ij} \f$ depends on the derivatives of the isochoric strain energy with respect to principal stretches and \f$ \mathbf{n}_a \f$ is the a-th principal vector
 *
 * @param dWdlambda_bar The derivative of the isochoric strain energy with respect to the first isochoric invariant.
 * @param dW2dlambda_bar2 The second derivative of the isochoric strain energy with respect to the first isochoric invariant. 
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
arma::mat L_vol_hyper(const double &dUdJ, const double &dU2dJ2, const arma::mat &b, const double &mJ = 0.);

} //namespace simcoon
