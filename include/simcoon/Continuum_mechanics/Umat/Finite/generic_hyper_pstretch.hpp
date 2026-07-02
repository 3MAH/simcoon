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

///@file generic_hyper_pstretch.hpp
///@brief User subroutine for generic hyperelastic materials using principal stretches
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file generic_hyper_pstretch.hpp
 * @brief Finite strain constitutive model.
 */

/** @addtogroup umat_finite
 *  @{
 */


/**
 * @brief Generic hyperelastic UMAT for potentials expressed in isochoric principal stretches \f$ \bar{\lambda}_a \f$.
 *
 * The potential is selected by umat_name. Currently available:
 * - OGDEN: \f$ W = \sum_{i=1}^N \frac{2 \mu_i}{\alpha_i^2} \left( \bar{\lambda}_1^{\alpha_i} + \bar{\lambda}_2^{\alpha_i} + \bar{\lambda}_3^{\alpha_i} - 3 \right) + \kappa \left( J \, \textrm{ln} J - J + 1 \right) \f$
 *   with props = { N, \f$ \kappa \f$, \f$ \mu_1 \f$, \f$ \alpha_1 \f$, ..., \f$ \mu_N \f$, \f$ \alpha_N \f$ } (nprops = 2 + 2N).
 *   The ground-state shear modulus is \f$ \mu = \sum_i \mu_i \f$; N=1, \f$ \alpha_1 = 2 \f$ recovers the compressible neo-Hookean potential (NEOHC).
 *
 * The Cauchy stress and the spatial tangent are assembled from the isochoric principal-stretch
 * machinery (sigma_iso_hyper_pstretch / L_iso_hyper_pstretch) plus the volumetric part
 * (sigma_vol_hyper / L_vol_hyper); the returned Lt follows the canonical box convention
 * \f$ \partial \hat{\tau} / \partial D_e \f$ (see generic_hyper_invariants).
 *
 * statev(0) stores the initial temperature; nstatev = 1.
 */
void umat_generic_hyper_pstretch(const std::string &umat_name, const arma::vec &etot, const arma::vec &Detot, const arma::mat &F0, const arma::mat &F1, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm_0, double &Wm_1, double &Wm_2, double &Wm_3, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = 0);
                            

/** @} */ // end of umat_finite group

} //namespace simcoon
