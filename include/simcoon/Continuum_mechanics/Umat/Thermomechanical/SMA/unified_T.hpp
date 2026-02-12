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

/**
 * @file unified_T.hpp
 * @brief Unified thermomechanical SMA model supporting multiple elastic symmetries and transformation criteria
 * @author Y. Chemisky, D. Chatziathanasiou
 * @version 1.0
 *
 * Based on the constitutive model of D. Chatziathanasiou Ph.D Thesis.
 * Thermomechanical extension of the unified SMA model with coupled heat equation.
 * Implemented in 1D-2D-3D.
 */

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Unified thermomechanical SMA model with coupled phase transformation and heat equation
 *
 * @details This function implements the fully coupled thermomechanical version of the unified
 * phenomenological SMA model. Compared to the mechanical-only version (umat_sma_unified_T),
 * this adds thermal coupling through the heat equation, requiring density and specific heat
 * capacity as additional material parameters, and producing thermal tangent operators as output.
 *
 * **Elastic Symmetry and Criteria Selection:**
 *
 * The elastic behavior and transformation criteria are determined by the umat_name parameter:
 *
 * | umat_name | Elasticity | Criteria | Number of props |
 * |-----------|------------|----------|-----------------|
 * | SMADI | Isotropic | Drucker | 31 |
 * | SMADC | Cubic | Drucker | 33 |
 * | SMAAI | Isotropic | Anisotropic Drucker (DFA) | 38 |
 * | SMAAC | Cubic | Anisotropic Drucker (DFA) | 40 |
 *
 * **Material Parameters (props):**
 *
 * **SMADI (Isotropic elasticity, Drucker criteria) — 31 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ \rho \f$ | Density | kg/m^3 |
 * | 1 | \f$ c_{pA} \f$ | Specific heat capacity of Austenite | J/(kg.K) |
 * | 2 | \f$ c_{pM} \f$ | Specific heat capacity of Martensite | J/(kg.K) |
 * | 3 | flagT | Temperature extrapolation (0=linear, 1=smooth) | - |
 * | 4 | \f$ E_A \f$ | Young's modulus of Austenite | MPa |
 * | 5 | \f$ E_M \f$ | Young's modulus of Martensite | MPa |
 * | 6 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 7 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 8 | \f$ \alpha_A \f$ | CTE of Austenite | 1/K |
 * | 9 | \f$ \alpha_M \f$ | CTE of Martensite | 1/K |
 * | 10 | \f$ H_{min} \f$ | Minimal transformation strain magnitude | - |
 * | 11 | \f$ H_{max} \f$ | Maximal transformation strain magnitude | - |
 * | 12 | \f$ k_1 \f$ | Exponential evolution parameter for \f$ H^{cur} \f$ | - |
 * | 13 | \f$ \sigma_{crit} \f$ | Critical stress for \f$ H^{cur} \f$ evolution | MPa |
 * | 14 | \f$ C_A \f$ | Clausius-Clapeyron slope (M \f$ \rightarrow \f$ A) | MPa/K |
 * | 15 | \f$ C_M \f$ | Clausius-Clapeyron slope (A \f$ \rightarrow \f$ M) | MPa/K |
 * | 16 | \f$ M_{s0} \f$ | Martensite start temperature at zero stress | K |
 * | 17 | \f$ M_{f0} \f$ | Martensite finish temperature at zero stress | K |
 * | 18 | \f$ A_{s0} \f$ | Austenite start temperature at zero stress | K |
 * | 19 | \f$ A_{f0} \f$ | Austenite finish temperature at zero stress | K |
 * | 20 | \f$ n_1 \f$ | Martensite start smooth exponent | - |
 * | 21 | \f$ n_2 \f$ | Martensite finish smooth exponent | - |
 * | 22 | \f$ n_3 \f$ | Austenite start smooth exponent | - |
 * | 23 | \f$ n_4 \f$ | Austenite finish smooth exponent | - |
 * | 24 | \f$ \sigma_{caliber} \f$ | Calibration stress for \f$ C_A \f$ / \f$ C_M \f$ | MPa |
 * | 25 | \f$ b \f$ | Tension-compression asymmetry parameter (Prager) | - |
 * | 26 | \f$ n \f$ | Tension-compression asymmetry exponent (Prager) | - |
 * | 27 | \f$ c_\lambda \f$ | Penalty function exponent start point | - |
 * | 28 | \f$ p_{0\lambda} \f$ | Penalty function limit value | - |
 * | 29 | \f$ n_\lambda \f$ | Penalty function power law exponent | - |
 * | 30 | \f$ \alpha_\lambda \f$ | Penalty function power law parameter | - |
 *
 * **SMADC (Cubic elasticity, Drucker criteria) — 33 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ \rho \f$ | Density | kg/m^3 |
 * | 1 | \f$ c_{pA} \f$ | Specific heat capacity of Austenite | J/(kg.K) |
 * | 2 | \f$ c_{pM} \f$ | Specific heat capacity of Martensite | J/(kg.K) |
 * | 3 | flagT | Temperature extrapolation (0=linear, 1=smooth) | - |
 * | 4 | \f$ E_A \f$ | Young's modulus of Austenite ([100] direction) | MPa |
 * | 5 | \f$ E_M \f$ | Young's modulus of Martensite ([100] direction) | MPa |
 * | 6 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 7 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 8 | \f$ G_A \f$ | Shear modulus of Austenite | MPa |
 * | 9 | \f$ G_M \f$ | Shear modulus of Martensite | MPa |
 * | 10–32 | | Same as SMADI props[8–30] (common SMA parameters) | |
 *
 * **SMAAI (Isotropic elasticity, anisotropic Drucker criteria) — 38 props:**
 *
 * Same as SMADI (props[0–30]) followed by 7 DFA parameters:
 *
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | 31 | \f$ F_{dfa} \f$ | F parameter of DFA criteria |
 * | 32 | \f$ G_{dfa} \f$ | G parameter of DFA criteria |
 * | 33 | \f$ H_{dfa} \f$ | H parameter of DFA criteria |
 * | 34 | \f$ L_{dfa} \f$ | L parameter of DFA criteria |
 * | 35 | \f$ M_{dfa} \f$ | M parameter of DFA criteria |
 * | 36 | \f$ N_{dfa} \f$ | N parameter of DFA criteria |
 * | 37 | \f$ K_{dfa} \f$ | K parameter of DFA criteria |
 *
 * **SMAAC (Cubic elasticity, anisotropic Drucker criteria) — 40 props:**
 *
 * Same as SMADC (props[0–32]) followed by 7 DFA parameters (props[33–39]).
 *
 * **State Variables (statev) — 17 variables:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ T_{init} \f$ | Initial/reference temperature | K |
 * | 1 | \f$ \xi \f$ | Martensitic volume fraction | - |
 * | 2 | \f$ \varepsilon^{tr}_{11} \f$ | Transformation strain component 11 | - |
 * | 3 | \f$ \varepsilon^{tr}_{22} \f$ | Transformation strain component 22 | - |
 * | 4 | \f$ \varepsilon^{tr}_{33} \f$ | Transformation strain component 33 | - |
 * | 5 | \f$ \gamma^{tr}_{12} \f$ | Transformation strain component 12 (engineering) | - |
 * | 6 | \f$ \gamma^{tr}_{13} \f$ | Transformation strain component 13 (engineering) | - |
 * | 7 | \f$ \gamma^{tr}_{23} \f$ | Transformation strain component 23 (engineering) | - |
 * | 8 | \f$ \xi_F \f$ | Forward martensitic volume fraction | - |
 * | 9 | \f$ \xi_R \f$ | Reverse martensitic volume fraction | - |
 * | 10 | \f$ \rho \Delta s_0 \f$ | Entropy difference between phases (M - A) | MPa/K |
 * | 11 | \f$ \rho \Delta E_0 \f$ | Internal energy difference between phases (M - A) | MPa |
 * | 12 | \f$ D \f$ | Stress dependence parameter for transformation limits | - |
 * | 13 | \f$ a_1 \f$ | Forward hardening parameter | MPa |
 * | 14 | \f$ a_2 \f$ | Reverse hardening parameter | MPa |
 * | 15 | \f$ a_3 \f$ | Equilibrium hardening parameter | MPa |
 * | 16 | \f$ Y_{0t} \f$ | Initial transformation critical value | MPa |
 *
 * @param umat_name Model variant (SMADI, SMADC, SMAAI, SMAAC)
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6x1)
 * @param DEtot Strain increment tensor (Voigt notation: 6x1)
 * @param sigma Cauchy stress tensor (Voigt notation: 6x1) [output]
 * @param r Heat source rate [output]
 * @param dSdE Mechanical tangent modulus \f$ \partial \sigma / \partial \varepsilon \f$ (6x6) [output]
 * @param dSdT Thermal stress tangent \f$ \partial \sigma / \partial T \f$ (6x1) [output]
 * @param drdE Thermal-mechanical tangent \f$ \partial r / \partial \varepsilon \f$ (1x6) [output]
 * @param drdT Thermal tangent \f$ \partial r / \partial T \f$ (scalar) [output]
 * @param DR Rotation increment matrix (3x3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (layout depends on umat_name)
 * @param nstatev Number of state variables (17)
 * @param statev State variables vector [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable (elastic) mechanical work [output]
 * @param Wm_ir Irrecoverable mechanical work stored in transformation [output]
 * @param Wm_d Dissipated mechanical work (hysteresis) [output]
 * @param Wt Total thermal work [output]
 * @param Wt_r Recoverable thermal work [output]
 * @param Wt_ir Irrecoverable thermal work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step ratio for adaptive time stepping [output]
 *
 * @note Elastic convention: isotropic uses "Enu", cubic uses "EnuG" for L_iso / L_cubic
 * @note The flagT parameter controls temperature extrapolation: 0 for linear, 1 for smooth
 * @note Legacy aliases SMAUT and SMANI map to SMADI and SMAAI respectively
 *
 * @see umat_sma_unified_T() for the mechanical-only version (no thermal coupling)
 * @see L_iso() for isotropic stiffness tensor (SMADI, SMAAI)
 * @see L_cubic() for cubic stiffness tensor (SMADC, SMAAC)
 *
 * **References:**
 * - Chatziathanasiou, D. (2016). Ph.D. Thesis — Phenomenological constitutive modeling
 *   of shape memory alloys.
 * - Chemisky, Y., Chatziathanasiou, D., Kumar, P., & Lagoudas, D. C. (2014).
 *   "A constitutive model for cyclic actuation of high-temperature shape memory alloys."
 *   *Mechanics of Materials*, 68, 120-136.
 */
void umat_sma_unified_T_T(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drdE, arma::mat &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
