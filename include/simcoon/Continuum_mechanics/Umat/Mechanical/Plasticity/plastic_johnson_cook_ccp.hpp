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
///@brief User subroutine for elastic-plastic materials with Johnson-Cook hardening in 1D-2D-3D case
///@brief This subroutine uses a convex cutting plane algorithm
///@brief Johnson-Cook rate-dependent and temperature-dependent hardening is considered
///@version 1.0

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file plastic_johnson_cook_ccp.hpp
 * @brief Elastic-plastic material model with Johnson-Cook hardening using the Convex Cutting Plane algorithm
 * @version 1.0
 */

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Elastic-plastic constitutive model with Johnson-Cook hardening solved by the Convex Cutting Plane (CCP) algorithm
 *
 * @details This function implements an elastic-plastic material model for small and finite strain analysis.
 * The model features:
 * - J2 (von Mises) plasticity with associative flow rule
 * - Johnson-Cook hardening with strain rate sensitivity and thermal softening:
 *   \f$ \sigma_Y = (A + B \cdot p^n)(1 + C \cdot \ln(\dot{\varepsilon}/\dot{\varepsilon}_0))(1 - T^{*m}) \f$
 *   where \f$ T^* = (T - T_{ref})/(T_{melt} - T_{ref}) \f$
 * - Fully implicit strain rate treatment within the CCP loop
 * - Consistent tangent modulus for implicit FE analysis
 *
 * **Material Properties (props):**
 * | Index | Symbol            | Description                          |
 * |-------|-------------------|--------------------------------------|
 * | 0     | E                 | Young's modulus                      |
 * | 1     | nu                | Poisson's ratio                      |
 * | 2     | alpha             | Coefficient of thermal expansion     |
 * | 3     | A                 | JC initial yield stress              |
 * | 4     | B                 | JC hardening coefficient             |
 * | 5     | n_JC              | JC hardening exponent                |
 * | 6     | C                 | JC strain rate sensitivity           |
 * | 7     | edot0             | JC reference strain rate             |
 * | 8     | m_JC              | JC thermal softening exponent        |
 * | 9     | T_ref             | JC reference temperature             |
 * | 10    | T_melt            | JC melting temperature               |
 *
 * **State Variables (statev):**
 * | Index | Symbol  | Description                              |
 * |-------|---------|------------------------------------------|
 * | 0     | T_init  | Initial temperature                      |
 * | 1     | p       | Accumulated plastic strain               |
 * | 2-7   | EP      | Plastic strain tensor (Voigt notation)   |
 * | 8     | edot_p  | Plastic strain rate (output/diagnostic)  |
 *
 * @param umat_name Name of the constitutive law (EPJCK)
 * @param Etot Total strain at t_n (6x1 Voigt)
 * @param DEtot Strain increment (6x1 Voigt)
 * @param sigma Cauchy stress [input/output] (6x1 Voigt)
 * @param Lt Consistent tangent modulus [output] (6x6)
 * @param L Elastic stiffness tensor (6x6)
 * @param DR Rotation increment (3x3)
 * @param nprops Number of material properties (11)
 * @param props Material properties vector
 * @param nstatev Number of state variables (9)
 * @param statev State variables vector [input/output]
 * @param T Temperature at t_n
 * @param DT Temperature increment
 * @param Time Time at t_n
 * @param DTime Time increment
 * @param Wm Total mechanical work [input/output]
 * @param Wm_r Recoverable mechanical work [input/output]
 * @param Wm_ir Irrecoverable mechanical work [input/output]
 * @param Wm_d Dissipated work [input/output]
 * @param ndi Number of direct stress components
 * @param nshr Number of shear stress components
 * @param start Flag for the first increment
 * @param tnew_dt Suggested time step ratio [output]
 */
void umat_plasticity_johnson_cook_CCP(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */

} //namespace simcoon
