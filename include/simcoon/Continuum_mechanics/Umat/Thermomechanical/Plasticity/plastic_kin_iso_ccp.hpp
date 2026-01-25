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

///@file plastic_isotropic_ccp.hpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Isotropic hardening with a power-law hardenig is considered
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file plastic_kin_iso_ccp.hpp
 * @brief Thermomechanical plasticity model.
 */

/** @addtogroup umat_thermomechanical
 *  @{
 */


///@brief The thermomechanical elastic-plastic UMAT with kinematic + isotropic hardening requires 9 constants:

///@brief props[0] : density
///@brief props[1] : specific heat capacity
///@brief props[2] : Young modulus
///@brief props[3] : Poisson ratio
///@brief props[4] : CTE
///@brief props[5] : J2 equivalent yield stress limit : sigmaY
///@brief props[6] : hardening parameter k
///@brief props[7] : exponent m
///@brief props[8] : linear kinematical hardening h

///@brief The elastic-plastic UMAT with kinematic + isotropic hardening requires 14 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)
///@brief statev[8] : Backstress 11: X(0,0)
///@brief statev[9] : Backstress 11: X(1,1)
///@brief statev[10] : Backstress 11: X(2,2)
///@brief statev[11] : Backstress 11: X(0,1)
///@brief statev[12] : Backstress 11: X(0,2)
///@brief statev[13] : Backstress 11: X(1,2)

void umat_plasticity_kin_iso_CCP_T(const arma::vec &, const arma::vec &, arma::vec &, double &, arma::mat &, arma::mat &, arma::mat &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);
    

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
