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
 * @param[out] N_1 3x3 matrix, spectral spin correction \f$ \sum_{i\neq j} f(b_i/b_j)\,\mathbf{B}_i\mathbf{D}\mathbf{B}_j \f$ from the eigenprojections \f$ \mathbf{B}_i \f$ of \f$ \mathbf{B}=\mathbf{F}\mathbf{F}^T \f$
 * @param[out] N_2 3x3 matrix, second spectral correction \f$ \sum_{i\neq j} g(b_i/b_j)\,\mathbf{B}_i\mathbf{D}\mathbf{B}_j \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] Omega 3x3 matrix representing spin rate \f$ \mathbf{\Omega}_{\textrm{log}} \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
 *
 * @details Example:
 * @code
 *      mat DR, N_1, N_2, D, Omega;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      logarithmic_R(DR, N_1, N_2, D, Omega, DTime, F0, F1);
 * @endcode
 */
void logarithmic_R(arma::mat &DR, arma::mat &N_1,  arma::mat &N_2, arma::mat &D,  arma::mat &Omega, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes the increment of the transformation gradient, the rate of deformation and the velocity gradient using the Truesdell rate.
 *
 * This function computes the increment of the transformation gradient \f$ \Delta \mathbf{F} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the velocity gradient \f$ \mathbf{L} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$ 
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the Truesdell rate and the time difference \f$ \Delta t \f$
 *
 * Note that this objective rate correspond to the covariant derivative
 * 
 * @param[out] DF 3x3 matrix representing the increment of transformation gradient \f$ \Delta \mathbf{F} \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] L 3x3 matrix representing the velocity gradient \f$ \mathbf{L} \f$
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
 * @brief Computes the increment of the transformation gradient, the rate of deformation, the spectral corrections and the velocity gradient using the convected logarithmic (log_F) framework.
 *
 * This function computes the increment of the transformation gradient \f$ \Delta \mathbf{F} \f$, the rate of deformation \f$ \mathbf{D} \f$ and the velocity gradient \f$ \mathbf{L} \f$ depending on \f$ \mathbf{F}_0 \f$ and \f$ \mathbf{F}_1 \f$
 * (\f$ \mathbf{F} \f$ at the beginning and end of an increment) using the convected logarithmic (log_F) framework and the time difference \f$ \Delta t \f$
 *
 * @param[out] DF 3x3 matrix representing the increment of transformation gradient \f$ \Delta \mathbf{F} \f$
 * @param[out] N_1 3x3 matrix, spectral spin correction \f$ \sum_{i\neq j} f(b_i/b_j)\,\mathbf{B}_i\mathbf{D}\mathbf{B}_j \f$ from the eigenprojections \f$ \mathbf{B}_i \f$ of \f$ \mathbf{B}=\mathbf{F}\mathbf{F}^T \f$
 * @param[out] N_2 3x3 matrix, second spectral correction \f$ \sum_{i\neq j} g(b_i/b_j)\,\mathbf{B}_i\mathbf{D}\mathbf{B}_j \f$
 * @param[out] D 3x3 matrix representing the rate of deformation \f$ \mathbf{D} \f$
 * @param[out] L 3x3 matrix representing the velocity gradient \f$ \mathbf{L} \f$
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @param[in] F0 transformation gradient \f$ \mathbf{F}_0 \f$ at time \f$ t_0 \f$
 * @param[in] F1 transformation gradient \f$ \mathbf{F}_1 \f$ at time \f$ t_1 \f$
 *
 * @details Example:
 * @code
 *      mat DF, N_1, N_2, D, L;
 *      mat F0 = randu(3,3);
 *      mat F1 = randu(3,3);
 *      double DTime = 0.1;
 *      logarithmic_F(DF, N_1, N_2, D, L, DTime, F0, F1);
 * @endcode
 */
void logarithmic_F(arma::mat &DF, arma::mat &N_1, arma::mat &N_2, arma::mat &D, arma::mat &L, const double &DTime, const arma::mat &F0, const arma::mat &F1);

/**
 * @brief Computes an incremental rotation (or deformation) using the Hughes-Winget formula.
 *
 * This function computes an incremental rotation \f$ \Delta \mathbf{R} \f$ from a spin (or velocity gradient) tensor \f$ \mathbf{\Omega} \f$ and a time increment \f$ \Delta t \f$
 * using the Hughes-Winget midpoint formula:
 *
 * \f[
 *      \Delta \mathbf{R} = \left( \mathbf{I} - \frac{1}{2} \Delta t \, \mathbf{\Omega} \right)^{-1} \left( \mathbf{I} + \frac{1}{2} \Delta t \, \mathbf{\Omega} \right)
 * \f]
 *
 * @param[in] Omega 3x3 matrix representing the spin tensor \f$ \mathbf{\Omega} \f$ (or velocity gradient)
 * @param[in] DTime time increment \f$ \Delta t \f$
 * @return 3x3 matrix representing the incremental rotation \f$ \Delta \mathbf{R} \f$
 *
 * @details Example:
 * @code
 *      mat Omega = randu(3,3);
 *      double DTime = 0.1;
 *      mat DR = Hughes_Winget(Omega, DTime);
 * @endcode
*/
arma::mat Hughes_Winget(const arma::mat &Omega, const double &DTime);

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
 *      mat BBBB = get_BBBB_GN(F);
 * @endcode
*/
arma::mat get_BBBB_GN(const arma::mat &F);

/**
 * @brief Strain-concentration tensor \f$ \mathbf{A}^{R} \f$ for the rotated (log_R) frame.
 *
 * Maps the rate of deformation to the \f$ \mathbf{R} \f$-corotational rate of the spatial Hencky
 * strain, \f$ \mathbf{D}_e=\mathbf{A}^{R}\!:\!\mathbf{D} \f$. In the eigenbasis of
 * \f$ \mathbf{B}=\mathbf{F}\mathbf{F}^T \f$ (eigenvalues \f$ b_a=\lambda_a^2 \f$), with
 * \f$ t=\ln(\lambda_i/\lambda_j) \f$, the spectral coefficients are \f$ A^{R}_{ij}=t/\sinh t \f$
 * (\f$ \to 1 \f$ on the diagonal) -- the geometric-mean logarithmic Daleckii-Krein kernel,
 * strictly positive at every stretch so \f$ \mathbf{A}^{R} \f$ is always invertible (Hoger's
 * tangent pushed forward by \f$ \mathbf{R} \f$). Returned in the engineering strain-concentration
 * Voigt convention (\f$ \mathbf{A}^{R}(\mathbf{I})=\mathbf{I}_6 \f$, rotates as
 * \f$ v_e\,\mathbf{A}^{R}\,v_s^T \f$ like a strain-concentration tensor), so apply as
 * \f$ \mathbf{D}_e=\mathrm{v2t\_strain}(\mathbf{A}^{R}\,\mathrm{t2v\_strain}(\mathbf{D})) \f$ and
 * invert as one (the stress dual is \f$ (\mathbf{A}^{R})^{T} \f$).
 * @param[in] F deformation gradient
 * @return the 6x6 (Voigt) strain-concentration tensor
*/
arma::mat A_R(const arma::mat &F);

/**
 * @brief Strain-concentration tensor \f$ \mathbf{A}^{F} \f$ for the convected (log_F) frame.
 *
 * Spectral construction on the eigenprojections of \f$ \mathbf{B}=\mathbf{F}\mathbf{F}^T \f$ with the
 * \f$ t\coth t \f$ kernel, \f$ t=\tfrac12\ln(b_i/b_j) \f$: \f$ A^{F}_{ij}=t\coth t \f$, \f$ A^{F}_{ii}=1 \f$.
 * The kernel is positive-definite, reduces to \f$ \mathbf{I}_6 \f$ at small strain and recovers
 * \f$ \ln V \f$ like @ref A_R. Same engineering strain-concentration convention and application as @ref A_R.
 * @param[in] F deformation gradient
 * @return the 6x6 (Voigt) strain-concentration tensor
*/
arma::mat A_F(const arma::mat &F);

/**
 * @brief Corotational logarithmic-strain increment by midpoint integration of the rate of deformation.
 *
 * Builds the incremental rotation \f$ \Delta \mathbf{R} \f$ from the spin \f$ \mathbf{\Omega} \f$ via the
 * Hughes-Winget midpoint formula and returns the corotational midpoint integral of the rate of
 * deformation \f$ \mathbf{D} \f$ over the step:
 * \f[
 *  \Delta \boldsymbol{\epsilon} = \tfrac12\left( \mathbf{D} + \Delta\mathbf{R}\,\mathbf{D}\,\Delta\mathbf{R}^T \right)\Delta t,
 *  \qquad \Delta\mathbf{R} = \left(\mathbf{I}-\tfrac{\Delta t}{2}\mathbf{\Omega}\right)^{-1}\left(\mathbf{I}+\tfrac{\Delta t}{2}\mathbf{\Omega}\right)
 * \f]
 *
 * @param[in] D rate of deformation \f$ \mathbf{D} \f$
 * @param[in] Omega corotational spin \f$ \mathbf{\Omega} \f$ (Jaumann / Green-Naghdi / logarithmic, per the calling rate)
 * @param[in] DTime time difference \f$ \Delta t = t_1 - t_0 \f$
 * @return the logarithmic strain increment \f$ \Delta \boldsymbol{\epsilon} \f$
*/
arma::mat Delta_log_strain(const arma::mat &D, const arma::mat &Omega, const double &DTime);

/**
 * @brief Naive log_F (convected) logarithmic-strain increment.
 *
 * Same midpoint form as Delta_log_strain, but the frame increment is the non-orthogonal
 * \f$ DF = (I-\tfrac{\Delta t}{2}L)^{-1}(I+\tfrac{\Delta t}{2}L) \f$, so the rotated term is the
 * push-forward \f$ DF\,D\,DF^{-1} \f$ — inverse, NOT transpose (F is not orthogonal).
 *
 * @param[in] D rate of deformation
 * @param[in] L velocity gradient
 * @param[in] DTime time difference \f$ \Delta t \f$
 * @return the naive log_F strain increment
*/
arma::mat Delta_log_strain_F(const arma::mat &D, const arma::mat &L, const double &DTime);

/**
 * @brief Corate-dispatched logarithmic-strain increment (the rate-form box accumulates
 *        \f$ \ln V = \tfrac12\ln(\mathbf{F}\mathbf{F}^T) \f$). Picks the integrator matching @p corate_type:
 *  - **2** XBM/logarithmic (\f$ \mathbf{A}=\mathbf{I} \f$): exact closed form
 *    \f$ \ln V_1 - \mathbf{DR}\,\ln V_0\,\mathbf{DR}^T \f$ (\f$ \epsilon=\ln V_1 \f$ to machine precision).
 *  - **3** log_R: \f$ \mathbf{D}_e=\mathbf{A}^{R}\!:\!\mathbf{D} \f$, the R-corotational (Green-Naghdi)
 *    rate of \f$ \ln V \f$ integrated over the orthogonal frame; \f$ \mathbf{A}^{R} \f$ is PD and
 *    \f$ \epsilon\to\ln V \f$, so the F-reconstruction stays well posed.
 *  - **5** log_F: \f$ \mathbf{D}_e=\mathbf{A}^{F}\!:\!\mathbf{D} \f$, the convected (Oldroyd) rate of
 *    \f$ \ln V \f$ integrated over the F-frame (\f$ \mathbf{DF} \f$ built from the velocity gradient
 *    L, carried in @p Omega). \f$ \mathbf{A}^{F} \f$ now recovers \f$ \ln V \f$ like \f$ \mathbf{A}^{R} \f$
 *    (the earlier \f$ -\tfrac12\ln(b_i b_j) \f$ indefinite term was removed). A genuine rate, used for
 *    ALL control_types so inelastic UMATs integrate from a real \f$ \mathbf{D}_e \f$.
 *  - **0/1/4** Jaumann / Green-Naghdi / Truesdell (\f$ \mathbf{A}=\mathbf{I} \f$): \f$ \mathbf{D}_e=\mathbf{D} \f$.
 *
 * Only affects rate-form/hypoelastic UMATs that accumulate \f$ \epsilon \f$; hyperelastic boxes read
 * stress/tangent off \f$ \mathbf{F}_1 \f$ and are unchanged.
 * @return the spatial logarithmic-strain increment for the chosen corate
*/
arma::mat Delta_log_strain_corate(const arma::mat &F0, const arma::mat &F1, const arma::mat &DR, const arma::mat &D, const arma::mat &Omega, const double &DTime, const int &corate_type);

/**
 * @brief Corate spin dispatch: for the chosen objective rate, set the frame increment @p DR and the
 *        rate of deformation @p D / spin (or velocity gradient L) @p Omega from @p F0, @p F1.
 *        Single source of truth for the solver's control_type ladders (predictor + Newton-Raphson),
 *        so every control_type (1..4, NLGEOM) dispatches the corate identically.
 *  - 0 Jaumann, 1 Green-Naghdi, 2 logarithmic/XBM, 3 log_R (DR=R-rotation), 4 Truesdell (DR=DF),
 *    5 log_F (DR=DF; @p Omega receives the velocity gradient L for the convected A^F:D rate).
 */
void corate_kinematics(const int &corate_type, arma::mat &DR, arma::mat &D, arma::mat &Omega, const arma::mat &F0, const arma::mat &F1, const double &DTime);

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
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 *
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * @param[in] DtauDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DtauDe_GreenNaghdiDD_2_DSDE(const arma::mat &DtauDe, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$ from the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin
 *
 * This function takes in the tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 *
 * It returns the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * @param[in] DsigmaDe (6x6 arma::mat) tangent modulus \f$ L^t \f$ that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and the approximation to logarithmic strain \f$ \mathbf{e} \f$ integrated using the Green-Naghdi spin
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
*/
arma::mat DsigmaDe_GreenNaghdiDD_2_DSDE(const arma::mat &DsigmaDe, const arma::mat &F, const arma::mat &sigma);

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

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 *
 * It returns the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
 *
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
*/
arma::mat DSDE_2_Dtau_GreenNaghdiDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Zaremba-Jaumann-Noll spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$,
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
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * This function takes in the tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$,
 * the transformation gradient \f$ \mathbf{F} \f$ and the Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 *
 * It returns the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
 *
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$ that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the Green-Naghdi spin
*/
arma::mat DSDE_2_Dsigma_GreenNaghdiDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Computes the tangent modulus that links the Kirchoff stress tensor \f$ \mathbf{\tau} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * Convenience function that computes get_BBBB(F) internally.
 *
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] tau (3x3 arma::mat) Kirchoff stress tensor \f$ \mathbf{\tau} \f$.
 * @return (6x6 arma::mat) the tangent modulus integrated using the logarithmic spin
*/
arma::mat DSDE_2_Dtau_logarithmicDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Computes the tangent modulus that links the Cauchy stress tensor \f$ \mathbf{\sigma} \f$ and rate of deformation \f$ \mathbf{D} \f$ integrated using the logarithmic spin from the tangent modulus that links the Piola-Kirchoff II stress \f$ \mathbf{S} \f$ to the Green-Lagrange stress \f$ \mathbf{E} \f$
 *
 * Convenience function that computes get_BBBB(F) internally.
 *
 * @param[in] DSDE (6x6 arma::mat) tangent modulus \f$ \frac{\partial \mathbf{S}}{\partial \mathbf{E}} \f$
 * @param[in] F (3x3 arma::mat) transformation gradient \f$ \mathbf{F} \f$
 * @param[in] sigma (3x3 arma::mat) Cauchy stress tensor \f$ \mathbf{\sigma} \f$.
 * @return (6x6 arma::mat) the tangent modulus integrated using the logarithmic spin
*/
arma::mat DSDE_2_Dsigma_logarithmicDD(const arma::mat &DSDE, const arma::mat &F, const arma::mat &sigma);

/**
 * @brief Corate-dispatched material<->box tangent maps. The box convention is
 * \f$ \mathbf{L}_t=\partial\hat{\boldsymbol\tau}/\partial\mathbf{D}_e \f$, the Kirchhoff corotational
 * tangent IN the solver's @p corate_type rate. Each picks the matching transport so the round-trip
 * \f$ \mathrm{d}\mathbf{S}/\mathrm{d}\mathbf{E}\leftrightarrow\mathbf{L}_t \f$ is exact per corate:
 * 0 Jaumann (spin W) | 1 Green-Naghdi | 2 logarithmic/XBM | 3 log_R (R = GN orthogonal frame) |
 * 5 log_F (convected/Oldroyd-Lie, pure F pull-back, B = I). @c DSDE_2_DtauDe_corate maps
 * \f$ \mathrm{d}\mathbf{S}/\mathrm{d}\mathbf{E}\to\mathbf{L}_t \f$; @c DtauDe_corate_2_DSDE is its inverse.
*/
arma::mat DSDE_2_DtauDe_corate(const arma::mat &DSDE, const int &corate_type, const arma::mat &F, const arma::mat &tau);
arma::mat DtauDe_corate_2_DSDE(const arma::mat &Lt, const int &corate_type, const arma::mat &F, const arma::mat &tau);

/**
 * @brief Assemble the canonical box tangent
 * \f$ \mathbf{L}_t=\partial\hat{\boldsymbol\tau}/\partial\mathbf{D}_e \f$ (Kirchhoff, no-J, XBM/log rate)
 * that every finite UMAT must emit -- the single source of truth for the box-tangent convention.
 * @c box_DtauDe_from_dSdE builds it from the material tangent \f$ \mathrm{d}\mathbf{S}/\mathrm{d}\mathbf{E} \f$;
 * @c box_DtauDe_from_spatial from the Cauchy (Oldroyd/Lie) spatial elasticity tensor
 * \f$ \partial\boldsymbol\sigma/\partial\mathbf{D} \f$.
*/
arma::mat box_DtauDe_from_dSdE(const arma::mat &dSdE, const arma::mat &F, const arma::vec &sigma);
arma::mat box_DtauDe_from_spatial(const arma::mat &Lt_spatial, const arma::mat &F, const arma::vec &sigma);

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
