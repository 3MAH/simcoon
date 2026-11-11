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

///@file elastic_orthotropic.hpp
///@brief User subroutine for ortothropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file elastic_orthotropic.hpp
 * @brief Thermomechanical elasticity model.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */


/**
 * @brief Thermomechanical orthotropic linear elasticity in the material axes (ELORT, thermomechanical).
 *
 * The stress is \f$ \boldsymbol{\sigma} = \mathbf{L} : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init}) \right) \f$,
 * the tangents \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} = \mathbf{L} \f$ and
 * \f$ \partial\boldsymbol{\sigma}/\partial T = -\mathbf{L} : \boldsymbol{\alpha} \f$; the reversible entropy is
 * \f$ \eta_r = \rho c_p \ln(T/T_{init}) + \boldsymbol{\alpha} : \boldsymbol{\sigma} \f$ (no dissipation),
 * with \f$ \mathbf{L} \f$ = L_ortho(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, "EnuG") and
 * \f$ \boldsymbol{\alpha} = (\alpha_1, \alpha_2, \alpha_3, 0, 0, 0) \f$.
 *
 * Voigt order [11, 22, 33, 12, 13, 23], engineering shear strains; MPa, K.
 *
 * **Material properties** (14):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | \f$ \rho \f$ | density |
 * | props[1] | \f$ c_p \f$ | specific heat capacity |
 * | props[2..4] | \f$ E_1, E_2, E_3 \f$ | Young moduli |
 * | props[5..7] | \f$ \nu_{12}, \nu_{13}, \nu_{23} \f$ | Poisson ratios |
 * | props[8..10] | \f$ G_{12}, G_{13}, G_{23} \f$ | shear moduli |
 * | props[11..13] | \f$ \alpha_1, \alpha_2, \alpha_3 \f$ | coefficients of thermal expansion |
 *
 * **State variables** (1):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | initial (reference) temperature |
 */
void umat_elasticity_ortho_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::vec &, arma::vec &, double &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);


/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
