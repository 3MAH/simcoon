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
///@brief Unified SMA model supporting multiple elastic symmetries and transformation criteria
///@brief Based on constitutive model of D. Chatziathanasiou Ph.D Thesis
///@brief Implemented by Y. Chemisky and D. Chatziathanasiou
///@brief Implemented in 1D-2D-3D

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

///@brief The unified SMA UMAT for transformation supports 4 variants via umat_name:
///@brief   SMADI: Drucker criteria, Isotropic elasticity (28 props)
///@brief   SMADC: Drucker criteria, Cubic elasticity (30 props)
///@brief   SMAAI: Anisotropic Drucker criteria, Isotropic elasticity (35 props)
///@brief   SMAAC: Anisotropic Drucker criteria, Cubic elasticity (37 props)
///@brief
///@brief For SMADI (isotropic elasticity, Drucker criteria) - 28 props:
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
///@brief props[11]: C_A : Slope of martensite -> austenite parameter
///@brief props[12]: C_M : Slope of austenite -> martensite parameter
///@brief props[13]: Ms0 : Martensite start at zero stress
///@brief props[14]: Mf0 : Martensite finish at zero stress
///@brief props[15]: As0 : Austenite start at zero stress
///@brief props[16]: Af0 : Austenite finish at zero stress
///@brief props[17]: n1 : Martensite start smooth exponent
///@brief props[18]: n2 : Martensite finish smooth exponent
///@brief props[19]: n3 : Austenite start smooth exponent
///@brief props[20]: n4 : Austenite finish smooth exponent
///@brief props[21]: sigmacaliber : Stress at which the slopes CA and CM are identified
///@brief props[22]: prager_b : Tension-compression asymmetry parameter
///@brief props[23]: prager_n : Tension-compression asymmetry exponent
///@brief props[24]: c_lambda : penalty function exponent start point
///@brief props[25]: p0_lambda : penalty function exponent limit penalty value
///@brief props[26]: n_lambda : penalty function power law exponent
///@brief props[27]: alpha_lambda : penalty function power law parameter
///@brief
///@brief For SMADC (cubic elasticity, Drucker criteria) - 30 props:
///@brief props[0] : flagT
///@brief props[1] : EA : Young's modulus of Austenite (in [100] direction)
///@brief props[2] : EM : Young's modulus of Martensite (in [100] direction)
///@brief props[3] : nuA : Poisson's ratio of Austenite
///@brief props[4] : nuM : Poisson's ratio of Martensite
///@brief props[5] : GA : Shear modulus of Austenite
///@brief props[6] : GM : Shear modulus of Martensite
///@brief props[7-29] : same as props[5-27] for SMADI
///@brief
///@brief For SMAAI (isotropic elasticity, anisotropic Drucker criteria) - 35 props:
///@brief Same as SMADI plus:
///@brief props[28]: F_dfa : F parameter of DFA criteria
///@brief props[29]: G_dfa : G parameter of DFA criteria
///@brief props[30]: H_dfa : H parameter of DFA criteria
///@brief props[31]: L_dfa : L parameter of DFA criteria
///@brief props[32]: M_dfa : M parameter of DFA criteria
///@brief props[33]: N_dfa : N parameter of DFA criteria
///@brief props[34]: K_dfa : K parameter of DFA criteria
///@brief
///@brief For SMAAC (cubic elasticity, anisotropic Drucker criteria) - 37 props:
///@brief Same as SMADC plus 7 DFA parameters at the end
///@brief
///@brief State variables (17 statev):
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
///@brief statev[11] : rhoDE0 difference in internal energy for the phases (M - A)
///@brief statev[12] : D parameter for the stress dependance of transformation limits
///@brief statev[13] : a1 : forward hardening parameter
///@brief statev[14] : a2 : reverse hardening parameter
///@brief statev[15] : a3 : Equilibrium hardening parameter
///@brief statev[16] : Y0t : Initial transformation critical value

void umat_sma_unified_T(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

} //namespace simcoon
