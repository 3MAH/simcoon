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

///@file unified_T.hpp
///@brief Unified model from:
///@brief Constitutive model of SMART LEM3 group - D. Chatziathanasiou, Y. Chemisky, G. Chatzigeorgiou, F. meraghni
///@brief Implemented in 1D-2D-3D

#pragma once
#include <armadillo>

namespace simcoon{

///@brief The unified SMA UMAT for transformation requires 25 constants and 17 statev:
///@brief The mechanical transformation UMAT for SMAs has the following statev and props

///@brief props[0] : flagT: 0 transformation temperatures linearly extrapolated; 1 : smooth temperatures
///@brief props[1] : EA: Young's modulus of Austenite
///@brief props[2] : EM: Young's modulus of Martensite
///@brief props[3] : nuA : Poisson's ratio of Austenite
///@brief props[4] : nuM : Poisson's ratio of Martensite
///@brief props[5] : alphaA_iso : CTE of Austenite
///@brief props[6] : alphaM_iso : CTE of Martensite
///@brief props[7] : Hmin : Minimal transformation strain magnitude
///@brief props[8] : Hmax : Maximal transformation strain magnitude
///@brief props[9] : k1 : Exponential evolution of transformation strain magnitude
///@brief props[10] : sigmacrit : Critical stress for change of transformation strain magnitude
///@brief props[11]: C_A : Slope of martesnite -> austenite parameter
///@brief props[12]: C_M : Slope of austenite -> martensite parameter
///@brief props[13]: Ms0 : Martensite start at zero stress
///@brief props[14]: Mf0 : Martensite finish at zero stress
///@brief props[15]: As0 : Austenite start at zero stress
///@brief props[16]: Af0 : Austenite finish at zero stress
///@brief props[17]: n1 : Martensite start smooth exponent
///@brief props[18]: n2 : Martensite finish smooth exponent
///@brief props[19]: n3 : Austenite start smooth exponent
///@brief props[20]: n4 : Austenite finish smooth exponent

///@brief props[21]: c_lambda : penalty function exponent start point
///@brief props[22]: p0_lambda : penalty function exponent limit penalty value
///@brief props[23]: n_lambda : penalty function power law exponent
///@brief props[24]: alpha_lambda : penalty function power law parameter

///@brief The elastic-plastic UMAT with isotropic hardening requires 14 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : xi : MVF (Martensitic volume fraction)
///@brief statev[2] : Transformation strain 11: ET(0,0)
///@brief statev[3] : Transformation strain 22: ET(1,1)
///@brief statev[4] : Transformation strain 33: ET(2,2)
///@brief statev[5] : Transformation strain 12: ET(0,1) (*2)
///@brief statev[6] : Transformation strain 13: ET(0,2) (*2)
///@brief statev[7] : Transformation strain 23: ET(1,2) (*2)

///@brief statev[8] : xiF : forward MVF
///@brief statev[9] : xiR : reverse MVF
///@brief statev[10] : rhoDs0 difference in entropy for the phases (M - A)
///@brief statev[11] : rhoDs0 difference in internal energy for the phases (M - A)
///@brief statev[12] : parameter for the stress dependance of transformation limits
///@brief statev[13] : a1 : forward hardening parameter
///@brief statev[14] : a2 : reverse hardening parameter
///@brief statev[15] : a3 : Equilibrium hardening parameter
///@brief statev[16] : Y0t : Initial transformation critical value

void umat_sma_unified_T(const arma::vec &, const arma::vec &, arma::vec &, arma::mat &, const arma::mat &, const int &, const arma::vec &, const int &, arma::vec &, const double &, const double &,const double &,const double &, double &, double &, double &, double &, const int &, const int &, const bool &, double &);
    
} //namespace simcoon
