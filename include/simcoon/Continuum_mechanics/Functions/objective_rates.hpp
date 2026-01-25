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
 * @brief The objective_rates library contains a set of functions that help to define different quantities, depending on a selected objective rate.
 */

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup objective_rates
 *  @{
 */

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
 * @param[out] Omega 3x3 matrix representing spin rate \f$ \mathbf{\Omega}_{\textrm{log}} \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
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
 * @param[out] Omega 3x3 matrix representing spin rate \f$ \mathbf{\Omega}_{\textrm{log}} \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
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
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
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
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
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
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
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
void logarithmic(arma::mat &DR, arma::mat &D, arma::mat &Omega, const double &DTime, const arma::mat &F0, const arma::mat &F1);

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
 * @param F the transformation gradient \f$ \mathbf{F} \f$
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
 * @param F the transformation gradient \f$ \mathbf{F} \f$
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
 * where \f$ F_1 \f$ and \f$ F_2 \f$ are the deformation gradient at the first and second times, respectively.
 * 
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @return The matrix representing the logarithmic strain increment \f$ \Delta \mathbf{e} \f$
*/
arma::mat Delta_log_strain(const arma::mat &F0, const arma::mat &F1, const double &DTime);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] DtauDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DtauDe_2_DSDE(const arma::mat &DtauDe, const arma::mat &BBBB, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis, 
 * the transformation gradient \f$ \mathbf{F} \f$
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Dtau_LieDD (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis, 
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat Dtau_LieDD_2_DSDE(const arma::mat &DtauDe, const arma::mat &F);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] DtauDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DtauDe_JaumannDD_2_DSDE(const arma::mat &DtauDe, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ integrated using the logarithmic spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] DsigmaDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$
 * @param[in] BBBB (6x6 arma::mat) logarithmic antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_2_DSDE(const arma::mat &DsigmaDe, const arma::mat &BBBB, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] DsigmaDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_2_DSDE(const arma::mat &DsigmaDe, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis, 
 * the transformation gradient \f$ \mathbf{F} \f$
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] Dtau_LieDD (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis, 
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat Dsigma_LieDD_2_DSDE(const arma::mat &DtauDe, const arma::mat &F);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * 
 * @param[in] DsigmaDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_JaumannDD_2_DSDE(const arma::mat &DsigmaDe, const arma::mat &F, const arma::mat &sigma);

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
 * @param[in] DtauDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] J (double) the jacobian of the transformation \f$ J = \textrm{det} \mathbf{F} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
*/
arma::mat DtauDe_2_DsigmaDe(const arma::mat &DtauDe, const double &J);

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
 * @param[in] DsigmaDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * @param[in] J (double) the jacobian of the transformation \f$ J = \textrm{det} \mathbf{F} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and logarithmic strain \f$ \mathbf{e} \f$ (or its approximation using an integration in an rotating frame)
*/
arma::mat DsigmaDe_2_DtauDe(const arma::mat &DsigmaDe, const double &J);

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
arma::mat DSDE_2_DtauDe(const arma::mat &DSDE, const arma::mat &BBBB, const arma::mat &F, const arma::mat &tau);

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
arma::mat DSDE_2_DsigmaDe(const arma::mat &DSDE, const arma::mat &BBBB, const arma::mat &F, const arma::mat &sigma);

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
arma::mat DSDE_2_Dtau_LieDD(const arma::mat &DSDE, const arma::mat &F);

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
arma::mat DSDE_2_Dsigma_LieDD(const arma::mat &DSDE, const arma::mat &F);

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
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$. 
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat DSDE_2_Dtau_JaumannDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &tau);

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
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$. 
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat DSDE_2_Dsigma_JaumannDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &sigma);

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
 * @param[in] Dtau_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat Dtau_LieDD_Dtau_JaumannDD(const arma::mat &Dtau_LieDD, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using an objective spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{obj}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \tau (p,s) \delta (i,r) - \frac{1}{2} \tau (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the objective antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$, 
 * and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * To determine the objective antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$, one can use the functions 
 * get_BBBB for logarithmic spin and get_BBBB_GN for Green-Naghdi spin
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using an objective spin
 * 
 * @param[in] Dtau_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] BBBB (6x6 arma::mat) objective antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the objective spin
*/
arma::mat Dtau_LieDD_Dtau_objectiveDD(const arma::mat &DSDE, const arma::mat &BBBB, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{GN}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \tau (p,s) \delta (i,r) - \frac{1}{2} \tau (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
 * 
 * @param[in] Dtau_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
*/
arma::mat Dtau_LieDD_Dtau_GreenNaghdiDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{log}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \tau (p,s) \delta (i,r) - \frac{1}{2} \tau (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
 * 
 * @param[in] Dtau_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
*/
arma::mat Dtau_LieDD_Dtau_logarithmicDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\mathbf{W})} = \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\mathbf{W})}
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
 * 
 * @note : The return tensor is the tangent modulus utilized by Abaqus
 * 
 * @param[in] Dsigma_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] sigma (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin
*/
arma::mat Dsigma_LieDD_Dsigma_JaumannDD(const arma::mat &Dsigma_LieDD, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using an objective spin from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{obj}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \sigma (p,s) \delta (i,r) - \frac{1}{2} \sigma (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the objective antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$, 
 * and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * 
 * To determine the objective antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$, one can use the functions 
 * get_BBBB for logarithmic spin and get_BBBB_GN for Green-Naghdi spin
 * 
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using an objective spin
 * 
 * @param[in] Dsigma_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] BBBB (6x6 arma::mat) objective antisymmetric tensor-valued function  \f$ \mathbf{\mathcal{B}}^{\textrm{obj}} \f$
 * @param[in] sigma (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the objective spin
*/
arma::mat Dsigma_LieDD_Dsigma_objectiveDD(const arma::mat &DSDE, const arma::mat &BBBB, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{GN}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \sigma (p,s) \delta (i,r) - \frac{1}{2} \sigma (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
 * 
 * @param[in] Dsigma_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
*/
arma::mat Dsigma_LieDD_Dsigma_GreenNaghdiDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis
 *
 * \f[
 *      \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} \mathbf{\Omega}_{\textrm{log}}} (i,s,p,r) = \left( \frac{\partial \mathbf{\eth \sigma}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} (i,s,p,r) - \frac{1}{2} \sigma (p,s) \delta (i,r) - \frac{1}{2} \sigma (r,s) \delta (i,p) - \frac{1}{2} \tau (i,r) \delta (s,p) - \frac{1}{2} \tau(i,p) \delta(s,r)
 * \f]
 * 
 * This function takes in the tangent modulus \f$ \left( \frac{\partial \mathbf{\eth \tau}}{\partial \eth \mathbf{e}} \right)_{\mathcal{R} (\vec{g}_i)} \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated in the natural covariant vector basis,
 * the logarithmic antisymmetric tensor-valued function \f$ \mathbf{\mathcal{B}}^{\textrm{log}} \f$, 
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\tau} \f$.
 * 
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
 * 
 * @param[in] Dsigma_LieDD (6x6 arma::mat) Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the natural covariant vector basis
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin
*/
arma::mat Dsigma_LieDD_Dsigma_logarithmicDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Biot stress tensor \f$ \mathbf{T} \f$ and rate of right stretch \f$ \mathbf{U} \f$ from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ 
 * 
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * the right stretch tensor \f$ \mathbf{U} \f$, 
 * and the Piola-Kirchoff II stress tensor  \f$ \mathbf{S} \f$.
 * 
 * It returns the tangent modulus that links the Biot stress tensor \f$ \mathbf{T} \f$ and rate of right stretch \f$ \mathbf{U} \f$
 * 
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] U(3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] S (3x3 arma::mat) Piola-Kirchoff II stress \f$ \mathbf{S} \f$ .
 * @return (6x6 arma::mat) the tangent modulus that links the Biot stress tensor \f$ \mathbf{T} \f$ and rate of right stretch \f$ \mathbf{U} \f$
*/
arma::mat DSDE_DBiotStressDU(const arma::mat &DSDE, const arma::mat &U, const arma::mat &S);

/** @} */ // end of objective_rates group

} //namespace simcoon
