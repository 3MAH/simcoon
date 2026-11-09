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

///@file elastic_isotropic.hpp
///@brief User thermomechanical subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file elastic_isotropic.hpp
 * @brief Thermomechanical elasticity model.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */

    
/**
 * @brief Thermomechanical isotropic linear elasticity (ELISO, thermomechanical).
 *
 * The stress is \f$ \boldsymbol{\sigma} = \mathbf{L} : \left( \boldsymbol{\varepsilon} - \boldsymbol{\alpha}\,(T - T_{init}) \right) \f$,
 * the tangents \f$ \partial\boldsymbol{\sigma}/\partial\boldsymbol{\varepsilon} = \mathbf{L} \f$ and
 * \f$ \partial\boldsymbol{\sigma}/\partial T = -\mathbf{L} : \boldsymbol{\alpha} \f$; the reversible entropy is
 * \f$ \eta_r = \rho c_p \ln(T/T_{init}) + \boldsymbol{\alpha} : \boldsymbol{\sigma} \f$ (no dissipation),
 * with \f$ \mathbf{L} = \mathbf{L}_{iso}(E, \nu) \f$ and \f$ \boldsymbol{\alpha} = \alpha\,\mathbf{I} \f$.
 *
 * Voigt order [11, 22, 33, 12, 13, 23], engineering shear strains; MPa, K.
 *
 * **Material properties** (5):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | props[0] | \f$ \rho \f$ | density |
 * | props[1] | \f$ c_p \f$ | specific heat capacity |
 * | props[2] | \f$ E \f$ | Young modulus |
 * | props[3] | \f$ \nu \f$ | Poisson ratio |
 * | props[4] | \f$ \alpha \f$ | coefficient of thermal expansion |
 *
 * **State variables** (1):
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | statev[0] | \f$ T_{init} \f$ | initial (reference) temperature |
 */
void umat_elasticity_iso_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &, const int & = 0);
    

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
