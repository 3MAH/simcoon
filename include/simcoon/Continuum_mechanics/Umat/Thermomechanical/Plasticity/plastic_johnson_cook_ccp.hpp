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

///@file plastic_johnson_cook_ccp.hpp
///@brief Thermomechanical UMAT for Johnson-Cook plasticity using the Convex Cutting Plane algorithm
///@brief Rate-dependent and temperature-dependent hardening with exact thermodynamic coupling
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @brief Thermomechanical Johnson-Cook plasticity with exact heat production (no Taylor-Quinney)
 *
 * @details Heat production r is computed from the Gibbs free energy framework,
 * providing exact decomposition of dissipated vs stored energy. The key difference
 * from EPICP_T is that dPhi/dT != 0 due to JC thermal softening (1 - T*^m),
 * which strengthens the thermomechanical coupling.
 *
 * **Material Properties (props, 13 constants):**
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | 0 | rho | Density |
 * | 1 | c_p | Specific heat capacity |
 * | 2 | E | Young's modulus |
 * | 3 | nu | Poisson's ratio |
 * | 4 | alpha | CTE |
 * | 5 | A | JC initial yield stress |
 * | 6 | B | JC hardening coefficient |
 * | 7 | n_JC | JC hardening exponent |
 * | 8 | C | JC strain rate sensitivity |
 * | 9 | edot0 | JC reference strain rate |
 * | 10 | m_JC | JC thermal softening exponent |
 * | 11 | T_ref | JC reference temperature |
 * | 12 | T_melt | JC melting temperature |
 */
void umat_plasticity_johnson_cook_CCP_T(const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drdE, arma::mat &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

} //namespace simcoon
