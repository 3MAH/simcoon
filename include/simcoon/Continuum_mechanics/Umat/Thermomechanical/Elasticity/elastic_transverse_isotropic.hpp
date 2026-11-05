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

///@file elastic_transverse_isotropic.hpp
///@brief User subroutine for transversely isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file elastic_transverse_isotropic.hpp
 * @brief Thermomechanical elasticity model.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */

    
/**
 * @brief Thermomechanical transversely isotropic linear elasticity (ELIST, thermomechanical).
 *
 * The stress is \f$ \boldsymbol{\sigma} = \mathbf{L} : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init}) \right) \f$,
 * the tangents \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} = \mathbf{L} \f$ and
 * \f$ \partial\boldsymbol{\sigma}/\partial T = -\mathbf{L} : \boldsymbol{\alpha} \f$; the reversible entropy is
 * \f$ \eta_r = \rho c_p \ln(T/T_{init}) + \boldsymbol{\alpha} : \boldsymbol{\sigma} \f$ (no dissipation),
 * with \f$ \mathbf{L} \f$ = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis) and
 * \f$ \boldsymbol{\alpha} = \alpha_T\,\mathbf{I} + (\alpha_L - \alpha_T)\,\mathbf{e}_{axis} \otimes \mathbf{e}_{axis} \f$.
 *
 * Voigt order [11, 22, 33, 12, 13, 23], engineering shear strains; MPa, K.
 *
 * **Material properties** (10):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | \f$ \rho \f$ | density |
 * | props[1] | \f$ c_p \f$ | specific heat capacity |
 * | props[2] | axis | longitudinal direction: 1, 2 or 3 |
 * | props[3] | \f$ E_L \f$ | longitudinal Young modulus |
 * | props[4] | \f$ E_T \f$ | transverse Young modulus |
 * | props[5] | \f$ \nu_{TL} \f$ | Poisson ratio, loading along the longitudinal axis |
 * | props[6] | \f$ \nu_{TT} \f$ | Poisson ratio in the transverse plane |
 * | props[7] | \f$ G_{LT} \f$ | longitudinal shear modulus |
 * | props[8] | \f$ \alpha_L \f$ | longitudinal coefficient of thermal expansion |
 * | props[9] | \f$ \alpha_T \f$ | transverse coefficient of thermal expansion |
 *
 * **State variables** (1):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | initial (reference) temperature |
 */
void umat_elasticity_trans_iso_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);
    

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
