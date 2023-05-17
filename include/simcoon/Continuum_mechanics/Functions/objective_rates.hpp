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
 * @file objective_rates.hpp
 * @author Yves Chemisky 
 * @section The objective_rates library contains a set of function that help to define different quantities, depending on a selected objective rate
*/

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @brief Computes the increment of rotation, the rate of deformation and the spin using the Jaumann corotational framework.
 *
 * This function computes the increment of rotation \f$ \Delta \mathbf{R} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the spin \f$ \mathbf{W} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * ( \f$ \mathbf{F} \f$ at the beginning and end of an increment) using the Jaumann corotational framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DR 3x3 matrix representing the increment of rotation \f$ \Delta \mathbf{R} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] W 3x3 matrix representing spin rate \f$ \mathbf{W} \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DR, D, W;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      Jaumann(R, DQ, D, DTime, F0, F1);
 * @endcode
 */
void Jaumann(arma::mat &DR, arma::mat &D,  arma::mat &W, const double &DTime, const arma::mat &F0, const arma::mat &F1);
    
/**
 * @brief Computes the increment of rotation, the rate of deformation and the spin using the Green-Naghdi corotational framework.
 *
 * This function computes the increment of rotation \f$ \Delta \mathbf{R} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the spin \f$ \mathbf{\Omega} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the Green-Naghdi corotational framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DR 3x3 matrix representing the increment of rotation \f$ \Delta \mathbf{R} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \mathbf{\Omega}_{\textrm{log}}
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DR, D, Omega;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      Green_Naghdi(DR, D, Omega, DTime, F0, F1);
 * @endcode
 */
void Green_Naghdi(arma::mat &DR, arma::mat &D,  arma::mat &Omega, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes the increment of rotation, the rate of deformation and the spin using the modified Logarithmic corotational framework using the "spin" \f$ \mathbf{\Omega} \f$ from \f$ \mathbf{R} \f$.
 *
 * This function computes the increment of rotation \f$ \Delta \mathbf{R} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the spin \f$ \mathbf{\Omega}_{\textrm{log}} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the modified Logarithmic corotational framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DR 3x3 matrix representing the increment of rotation \f$ \Delta \mathbf{R} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \mathbf{\Omega}_{\textrm{log}}
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DR, D, Omega;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      logarithmic_R(DR, D, Omega, DTime, F0, F1);
 * @endcode
 */
void logarithmic_R(arma::mat &DR, arma::mat &N_1,  arma::mat &N_2, arma::mat &D,  arma::mat &Omega, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes the increment of the velocity gradient, the rate of deformation and the velocity gradient using the Truesdell rate.
 *
 * This function computes the increment of the transformation gradient \f$ \Delta \mathbf{F} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the velocity gradient \f$ \mathbf{L} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the modified Logarithmic corotational framework and the time difference \f$ \Delta t \f$
 *
 * Note that this objective rate correspond to the covariant derivative
 * 
 * @param[out] DF 3x3 matrix representing the increment of transformation gradient \f$ \Delta \mathbf{F} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \mathbf{\L}
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DF, D, L;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      Truesdell(DF, D, L, DTime, F0, F1);
 * @endcode
 */
void Truesdell(arma::mat &DF, arma::mat &D, arma::mat &L, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes the increment of rotation, the rate of deformation and the spin using the modified Logarithmic corotational framework using the "spin" L
 *
 * This function computes the increment of rotation \f$ \Delta \mathbf{R} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the spin \f$ \mathbf{\Omega}_{\textrm{log}} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the modified Logarithmic corotational framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DR 3x3 matrix representing the increment of rotation \f$ \Delta \mathbf{R} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \mathbf{\Omega}_{\textrm{log}}
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DR, D, Omega;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      logarithmic_F(DR, D, Omega, DTime, F0, F1);
 * @endcode
 */
void logarithmic_F(arma::mat &DF, arma::mat &N_1, arma::mat &N_2, arma::mat &D, arma::mat &L, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes the increment of rotation, the rate of deformation and the spin using the Logarithmic corotational framework using the "spin" \f$ \mathbf{\Omega}_{\textrm{log}} \f$.
 *
 * This function computes the increment of rotation \f$ \Delta \mathbf{R} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the spin \f$ \mathbf{\Omega}_{\textrm{log}} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the modified Logarithmic corotational framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DR 3x3 matrix representing the increment of rotation \f$ \Delta \mathbf{R} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \mathbf{\Omega}_{\textrm{log}}
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 *
 * @details Example: 
 * @code
 *      mat DR, D, Omega;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      logarithmic_R(DR, D, Omega, DTime, F0, F1);
 * @endcode
*/
void logarithmic(arma::mat &, arma::mat &, arma::mat &, const double &, const arma::mat &, const arma::mat &);

/**
 * @brief This function computes the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$ from the transformation gradient \f$ \mathbf{F} \f$.
 * 
 * \f[
 *      \mathbf{B}_i = b_i \otimes b_i\,  \quad f(z) =  \frac{ 1+ \frac{\lambda_i}{\lambda_j} }{1 - \frac{\lambda_i}{\lambda_j}} + 2\, \textrm{ln} \left( \frac{\lambda_i}{\lambda_j} \right)
 * \f]
 * 
 * \f[
 *      \mathbf{\mathcal{B}}^{\textrm{log}} = \sum_{i \neq j} f(z) \, \mathbf{\beta} \left( b_i, b_j \right)
 * \f] 
 * 
 * where \f$ b_i \f$ is the i-th eigenvector of the left Cauchy Green strain tensor \f$ \mathbf{B} = \mathbf{F} \cdot \mathbf{F}^T \f$ \cite Xiao.etal.1998,
 * and \f$ \mathbf{\beta} \left( b_i, b_j \right) \f$ is a spectral operator 
 * 
 * \f[
 *      \mathbf{\beta} \left( b_i, b_j \right) = \frac{1}{2} \left( b_i \otimes b_j \otimes b_i \otimes b_j + b_i \otimes b_j \otimes b_j \otimes b_i \right)
 * \f]
 *
 * Note that this allows to compute a correction \f$ \mathbf{\mathcal{N}}^{\textrm{log}} \f$ to the spin that depends on \f$ \mathbf{D} \f$ such that \f$ \mathbf{\mathcal{N}}^{\textrm{log}} = \mathbf{\mathcal{B}} : \mathbf{D} \f$
 * The logarithmic spin is therefore :
 * 
 * \f[
 *      \mathbf{\Omega}^{\textrm{log}} = \mathbf{W} + \mathbf{\mathcal{N}}^{\textrm{log}}	
 * \f]
 * 
 * @param B the Left Cauchy-Green strain \f$ \mathbf{B} \f$
 * @return 6x6 matrix representing the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat BBBB = get_BBBB(F);
 * @endcode
*/
arma::mat get_BBBB(const arma::mat &F);

/**
 * @brief This function computes the Green-Naghdi antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{GN}} \f$ from the transformation gradient \f$ \mathbf{F} \f$.
 * 
 * \f[
 *      \mathbf{B}_i = b_i \otimes b_i\,  \quad f(z) =  \frac{ \sqrt{\lambda_j} - \sqrt{\lambda_i} }{\sqrt{\lambda_j} + \sqrt{\lambda_i}}
 * \f]
 * 
 * (sqrt(bi(j)) - sqrt(bi(i)))/(sqrt(bi(j)) + sqrt(bi(i)));
 * 
 * \f[
 *      \mathbf{\mathcal{B}}^{\textrm{GN}} = \sum_{i \neq j} f(z) \, \mathbf{\beta} \left( b_i, b_j \right)
 * \f] 
 * 
 * where \f$ b_i \f$ is the i-th eigenvector of the left Cauchy Green strain tensor \f$ \mathbf{B} = \mathbf{F} \cdot \mathbf{F}^T \f$ \cite Xiao.etal.1998,
 * and \f$ \mathbf{\beta} \left( b_i, b_j \right) \f$ is a spectral operator 
 * 
 * \f[
 *      \mathbf{\beta} \left( b_i, b_j \right) = \frac{1}{2} \left( b_i \otimes b_j \otimes b_i \otimes b_j + b_i \otimes b_j \otimes b_j \otimes b_i \right)
 * \f]
 *
 * Note that this allows to compute a correction \f$ \mathbf{\mathcal{N}}^{\textrm{GN}} \f$ to the spin that depends on \f$ \mathbf{D} \f$ such that \f$ \mathbf{\mathcal{N}}^{\textrm{GN}} = \mathbf{\mathcal{B}} : \mathbf{D} \f$
 * The logarithmic spin is therefore :
 * 
 * \f[
 *      \mathbf{\Omega}^{\textrm{GN}} = \mathbf{W} + \mathbf{\mathcal{N}}^{\textrm{GN}}	
 * \f]
 * 
 * @param B the Left Cauchy-Green strain \f$ \mathbf{B} \f$
 * @return 6x6 matrix representing the Green-Naghdi antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{GN}} \f$
 *
 * @details Example: 
 * @code
 *      mat F = randu(3,3);
 *      mat BBBB = get_BBBB(F);
 * @endcode
*/
arma::mat get_BBBB_GN(const arma::mat &F);

/**
 * @brief Computes the logarithmic strain increment
 *
 * This function takes in two matrices representing the deformation gradient at two different times, \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$ and \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
 * the time difference \f$ \Delta t = t_1 - t_0 \f$
 * It returns the matrix representing the logarithmic strain increment.
 *
 * The logarithmic strain increment is calculated using the following equation:
 * \f[
 *  \Delta \epsilon^{\text{log}} = \frac{1}{2}\left( \ln\left(F_2^TF_2\right) - \ln\left(F_1^TF_1\right) \right) 
 * \f]
 *
 * where \f$F_1\f$ and \f$F_2\f$ are the deformation gradient at the first and second times, respectively.
 * 
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time $t_0$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time $t_1$
 * @return The matrix representing the logarithmic strain increment \f$ \Delta \mathbf{e} \f$
*/
arma::mat Delta_log_strain(const arma::mat &F0, const arma::mat &F1, const double &DTime);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}}, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DtauDe_2_DSDE(const arma::mat &Lt, const arma::mat &BBBB, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DtauDe_JaumannDD_2_DSDE(const arma::mat &Lt, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}}, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_2_DSDE(const arma::mat &Lt, const arma::mat &BBBB, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_2_DSDE(const arma::mat &Lt, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_JaumannDD_2_DSDE(const arma::mat &, const arma::mat &, const arma::mat &);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame) from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
 *
 * This function takes in tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame) and 
 * the jacobian of the transformation \f$ J = \textrm{det} \mathbf{F} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
 * 
 * \f[
 *      \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} = \frac{1}{J} \frac{\partial \eth \mathbf{\tau}}{\partial \eth \mathbf{e}}  
 * \f]
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
*/
arma::mat DtauDe_2_DsigmaDe(const arma::mat &Lt, const double &J);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame) from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
 * 
 * This function takes in tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame) and 
 * the jacobian of the transformation \f$ J = \textrm{det} \mathbf{F} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
 * 
 * \f[
 *      \frac{\partial \eth \mathbf{\tau}}{\partial \eth \mathbf{e}} = J \frac{\partial \eth \mathbf{\sigma}}{\partial \eth \mathbf{e}} 
 * \f]
 * 
 * @param[in] Lt (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] J (double) the jacobian of the transformation \f$ J = \textrm{det} \mathbf{F} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
*/
arma::mat DsigmaDe_2_DtauDe(const arma::mat &Lt, const double &J);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
*/
arma::mat DSDE_2_DtauDe(const arma::mat &Lt, const arma::mat &BBBB, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
*/
arma::mat DSDE_2_DsigmaDe(const arma::mat &Lt, const arma::mat &BBBB, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
*/
arma::mat DSDE_2_Dtau_LieDD(const arma::mat &, const arma::mat &);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
*/
arma::mat DSDE_2_Dsigma_LieDD(const arma::mat &, const arma::mat &);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat DSDE_2_Dtau_JaumannDD(const arma::mat &Lt, const arma::mat &F, const arma::mat &tau);

//arma::mat DSDE_2_Dtau_GreenNaghdiDD(const arma::mat &Lt, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat DSDE_2_Dsigma_JaumannDD(const arma::mat &Lt, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\mathbf{W})} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \tau (p,s) \delta (i,r) - \frac{1}{2} \tau (r,s) \delta (i,p) \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat Dtau_LieDD_Dtau_JaumannDD(const arma::mat &Lt, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{log}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \tau (p,s) \delta (i,r) - \frac{1}{2} \tau (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
*/
arma::mat Dtau_LieDD_Dtau_logarithmicDD(const arma::mat &Lt, const arma::mat &BBBB, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\mathbf{W})} = \frac{1}{J} \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\mathbf{W})}
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * 
 * @note : The return tensor is the tangent modulus utilized by Abaqus
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat Dtau_LieDD_Dsigma_JaumannDD(const arma::mat &Dtau_JaumannDD, const double &J);
 
} //namespace simcoon
