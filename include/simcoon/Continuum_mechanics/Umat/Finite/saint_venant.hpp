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

///@file saint_venant.hpp
///@brief User subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/**
 * @file saint_venant.hpp
 * @brief Finite strain constitutive model.
 */

/** @addtogroup umat_finite
 *  @{
 */


/**
 * @brief Saint-Venant-Kirchhoff hyperelastic UMAT.
 *
 * \f$ \mathbf{S} = \mathbf{L} : \mathbf{E} \f$ with \f$ \mathbf{E} \f$ the Green-Lagrange
 * strain and \f$ \mathbf{L} \f$ the isotropic stiffness of \f$ (E, \nu) \f$. The kernel is
 * Kirchhoff-native: @p sigma carries \f$ \boldsymbol{\tau} = \mathbf{F}\mathbf{S}\mathbf{F}^T \f$,
 * and the tangent is \f$ \mathbf{L} \f$ (which is \f$ \partial \mathbf{S} / \partial \mathbf{E} \f$)
 * mapped once to the box of @p corate_type by DSDE_2_DtauDe_corate. The law is elastic, so
 * \f$ W_m = W_{m,r} \f$ and \f$ W_{m,ir} = W_{m,d} = 0 \f$.
 *
 * **props** (3): \f$ E, \nu, \alpha \f$. The CTE \f$ \alpha \f$ is read but not used: no
 * thermal strain is applied. **statev** (1): @c statev(0) stores the initial temperature.
 *
 * @param[in] umat_name the 5-letter UMAT name
 * @param[in] etot,Detot logarithmic strain at the start of the increment and its increment
 * @param[in] F0,F1 deformation gradient at the start and the end of the increment
 * @param[in,out] sigma Kirchhoff stress \f$ \boldsymbol{\tau} \f$, 6-Voigt (the UMAT interface name)
 * @param[out] Lt box tangent \f$ \partial \hat{\boldsymbol{\tau}} / \partial \mathbf{D}_e \f$ in @p corate_type, 6x6
 * @param[out] L the ground-state stiffness
 * @param[in] DR the increment of rigid-body rotation
 * @param[in] nprops,props the material parameters (see above)
 * @param[in] nstatev,statev the state variables (one, the initial temperature)
 * @param[in] T,DT temperature and its increment
 * @param[in] Time,DTime time and its increment
 * @param[in,out] Wm,Wm_r,Wm_ir,Wm_d the cumulative work terms, per reference volume
 * @param[in] ndi,nshr the number of direct and shear components
 * @param[in] start true on the first call of a block
 * @param[out] tnew_dt the suggested time-step scaling
 * @param[in] corate_type the solver's objective rate, which @p Lt is expressed in
 * @param[in] tangent_mode unused: a hyperelastic law returns its exact tangent
 */
void umat_saint_venant(const std::string &umat_name, const arma::vec &etot, const arma::vec &Detot, const arma::mat &F0, const arma::mat &F1, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode = tangent_default);
                            

/** @} */ // end of umat_finite group

} //namespace simcoon
