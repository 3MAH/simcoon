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

///@file hypoelastic_orthotropic.hpp
///@brief User subroutine for ortothropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/**
 * @file hypoelastic_orthotropic.hpp
 * @brief Finite strain constitutive model.
 */

/** @addtogroup umat_finite
 *  @{
 */


/**
 * @brief Hypoelastic orthotropic finite-strain UMAT (rate form).
 *
 * Integrates the corotational KIRCHHOFF rate incrementally,
 * \f[
    \boldsymbol{\tau}_{n+1} = \boldsymbol{\tau}_n + \mathbf{L} : \left( \Delta\boldsymbol{\varepsilon}
    - \boldsymbol{\alpha} \Delta T \right),
 * \f]
 * where \f$ \boldsymbol{\tau}_n \f$ has already been transported by the solver's corate.
 * It is Kirchhoff-native like every other simcoon kernel: @p sigma carries
 * \f$ \boldsymbol{\tau} \f$, and the Cauchy stress \f$ \boldsymbol{\sigma} =
 * \boldsymbol{\tau}/J \f$ is formed only at the output boundaries.
 *
 * It is the RATE counterpart of ELORT, which evaluates the same orthotropic
 * \f$ \mathbf{L} \f$ in total form, \f$ \boldsymbol{\tau} = \mathbf{L} :
 * (\boldsymbol{\varepsilon} - \boldsymbol{\alpha} \Delta T) \f$. Along a path without rotation
 * increments the two are identical. Under rotation (shear, off-axis loading) they differ in
 * two ways:
 * - where the stress transport does not commute with an anisotropic \f$ \mathbf{L} \f$ -- the
 *   difference this kernel exists to show, of order one for a strongly orthotropic stiffness;
 * - by the time discretisation of the rate form, which sums transported increments where
 *   ELORT evaluates the total strain. That part is first order in the increment size and
 *   present even for an isotropic \f$ \mathbf{L} \f$ (about 1 % at \f$ \Delta\gamma = 0.01 \f$
 *   in simple shear); it vanishes as the increments are refined.
 *
 * **props** (12): \f$ E_x, E_y, E_z, \nu_{xy}, \nu_{xz}, \nu_{yz}, G_{xy}, G_{xz}, G_{yz},
 * \alpha_x, \alpha_y, \alpha_z \f$ ("EnuG" convention, material frame).
 *
 * **statev** (1): @c statev(0) stores the initial temperature. The law is elastic, so
 * \f$ W_m = W_{m,r} \f$ (per reference volume, on \f$ \boldsymbol{\tau} \f$) and
 * \f$ W_{m,ir} = W_{m,d} = 0 \f$.
 *
 * @p Lt is \f$ \mathbf{L} \f$ itself: for a Kirchhoff rate it is already the box
 * \f$ \partial \hat{\boldsymbol{\tau}} / \partial \mathbf{D}_e \f$, and it is in-rate whatever
 * the corate, since the increment arrives corotated. @p corate_type, @p F0 and @p F1 are
 * therefore unused.
*/
void umat_hypoelasticity_ortho(const std::string &umat_name, const arma::vec &etot, const arma::vec &Detot, const arma::mat &F0, const arma::mat &F1, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode = tangent_default);
                              

/** @} */ // end of umat_finite group

} //namespace simcoon
