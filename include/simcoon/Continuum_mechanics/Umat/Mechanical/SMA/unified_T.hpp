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
 * @brief Unified SMA model supporting multiple elastic symmetries and transformation criteria
 * @author Y. Chemisky, D. Chatziathanasiou
 * @version 1.0
 *
 * Based on the constitutive model of D. Chatziathanasiou Ph.D Thesis.
 * Implemented in 1D-2D-3D.
 */

#pragma once
#include <string>
#include <armadillo>

namespace simcoon{

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Unified phenomenological SMA model with phase transformation
 *
 * @details This function implements the unified phenomenological model for shape memory alloys
 * based on Chatziathanasiou's framework. The model describes the macroscopic thermomechanical
 * behavior through a single scalar internal variable (martensitic volume fraction \f$ \xi \f$)
 * and a tensorial transformation strain \f$ \boldsymbol{\varepsilon}^{tr} \f$.
 *
 * **Elastic Symmetry and Criteria Selection:**
 *
 * The elastic behavior and transformation criteria are determined by the umat_name parameter:
 *
 * | umat_name | Elasticity | Criteria | Number of props |
 * |-----------|------------|----------|-----------------|
 * | SMADI | Isotropic | Drucker | 28 |
 * | SMADC | Cubic | Drucker | 30 |
 * | SMAAI | Isotropic | Anisotropic Drucker (DFA) | 35 |
 * | SMAAC | Cubic | Anisotropic Drucker (DFA) | 37 |
 *
 * **Effective Properties:**
 *
 * The effective compliance is obtained via a Reuss mixing rule:
 * \f[
 * \mathbf{M}_{eff} = \xi \mathbf{M}_M + (1 - \xi) \mathbf{M}_A
 * \f]
 * and the effective CTE is similarly interpolated.
 *
 * **Current Maximum Transformation Strain:**
 *
 * The maximum transformation strain magnitude \f$ H^{cur} \f$ evolves with stress:
 * \f[
 * H^{cur}(\bar{\sigma}) = H_{min} + (H_{max} - H_{min})(1 - e^{-k_1 \bar{\sigma} / \sigma_{crit}})
 * \f]
 *
 * **Material Parameters (props):**
 *
 * **SMADI (Isotropic elasticity, Drucker criteria) — 28 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | flagT | Temperature extrapolation (0=linear, 1=smooth) | - |
 * | 1 | \f$ E_A \f$ | Young's modulus of Austenite | MPa |
 * | 2 | \f$ E_M \f$ | Young's modulus of Martensite | MPa |
 * | 3 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 4 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 5 | \f$ \alpha_A \f$ | CTE of Austenite | 1/K |
 * | 6 | \f$ \alpha_M \f$ | CTE of Martensite | 1/K |
 * | 7 | \f$ H_{min} \f$ | Minimal transformation strain magnitude | - |
 * | 8 | \f$ H_{max} \f$ | Maximal transformation strain magnitude | - |
 * | 9 | \f$ k_1 \f$ | Exponential evolution parameter for \f$ H^{cur} \f$ | - |
 * | 10 | \f$ \sigma_{crit} \f$ | Critical stress for \f$ H^{cur} \f$ evolution | MPa |
 * | 11 | \f$ C_A \f$ | Clausius-Clapeyron slope (M \f$ \rightarrow \f$ A) | MPa/K |
 * | 12 | \f$ C_M \f$ | Clausius-Clapeyron slope (A \f$ \rightarrow \f$ M) | MPa/K |
 * | 13 | \f$ M_{s0} \f$ | Martensite start temperature at zero stress | K |
 * | 14 | \f$ M_{f0} \f$ | Martensite finish temperature at zero stress | K |
 * | 15 | \f$ A_{s0} \f$ | Austenite start temperature at zero stress | K |
 * | 16 | \f$ A_{f0} \f$ | Austenite finish temperature at zero stress | K |
 * | 17 | \f$ n_1 \f$ | Martensite start smooth exponent | - |
 * | 18 | \f$ n_2 \f$ | Martensite finish smooth exponent | - |
 * | 19 | \f$ n_3 \f$ | Austenite start smooth exponent | - |
 * | 20 | \f$ n_4 \f$ | Austenite finish smooth exponent | - |
 * | 21 | \f$ \sigma_{caliber} \f$ | Calibration stress for \f$ C_A \f$ / \f$ C_M \f$ | MPa |
 * | 22 | \f$ b \f$ | Tension-compression asymmetry parameter (Prager) | - |
 * | 23 | \f$ n \f$ | Tension-compression asymmetry exponent (Prager) | - |
 * | 24 | \f$ c_\lambda \f$ | Penalty function exponent start point | - |
 * | 25 | \f$ p_{0\lambda} \f$ | Penalty function limit value | - |
 * | 26 | \f$ n_\lambda \f$ | Penalty function power law exponent | - |
 * | 27 | \f$ \alpha_\lambda \f$ | Penalty function power law parameter | - |
 *
 * **SMADC (Cubic elasticity, Drucker criteria) — 30 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | flagT | Temperature extrapolation (0=linear, 1=smooth) | - |
 * | 1 | \f$ E_A \f$ | Young's modulus of Austenite ([100] direction) | MPa |
 * | 2 | \f$ E_M \f$ | Young's modulus of Martensite ([100] direction) | MPa |
 * | 3 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 4 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 5 | \f$ G_A \f$ | Shear modulus of Austenite | MPa |
 * | 6 | \f$ G_M \f$ | Shear modulus of Martensite | MPa |
 * | 7–29 | | Same as SMADI props[5–27] (common SMA parameters) | |
 *
 * **SMAAI (Isotropic elasticity, anisotropic Drucker criteria) — 35 props:**
 *
 * Same as SMADI (props[0–27]) followed by 7 DFA parameters:
 *
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | 28 | \f$ F_{dfa} \f$ | F parameter of DFA criteria |
 * | 29 | \f$ G_{dfa} \f$ | G parameter of DFA criteria |
 * | 30 | \f$ H_{dfa} \f$ | H parameter of DFA criteria |
 * | 31 | \f$ L_{dfa} \f$ | L parameter of DFA criteria |
 * | 32 | \f$ M_{dfa} \f$ | M parameter of DFA criteria |
 * | 33 | \f$ N_{dfa} \f$ | N parameter of DFA criteria |
 * | 34 | \f$ K_{dfa} \f$ | K parameter of DFA criteria |
 *
 * **SMAAC (Cubic elasticity, anisotropic Drucker criteria) — 37 props:**
 *
 * Same as SMADC (props[0–29]) followed by 7 DFA parameters (props[30–36]).
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
 * @param Lt Consistent tangent modulus (6x6) [output]
 * @param L Elastic stiffness tensor (6x6) [output]
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
 * @param Wm_r Recoverable (elastic) work [output]
 * @param Wm_ir Irrecoverable work stored in transformation [output]
 * @param Wm_d Dissipated work (hysteresis) [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step ratio for adaptive time stepping [output]
 *
 * @note Elastic convention: isotropic uses "Enu", cubic uses "EnuG" for L_iso / L_cubic
 * @note The flagT parameter controls temperature extrapolation: 0 for linear, 1 for smooth
 * @note Legacy aliases SMAUT and SMANI map to SMADI and SMAAI respectively
 *
 * @see L_iso() for isotropic stiffness tensor (SMADI, SMAAI)
 * @see L_cubic() for cubic stiffness tensor (SMADC, SMAAC)
 * @see umat_sma_mono() for the micromechanical monocrystal SMA model
 *
 * **References:**
 * - Chatziathanasiou, D. (2016). Ph.D. Thesis — Phenomenological constitutive modeling
 *   of shape memory alloys.
 * - Chemisky, Y., Chatziathanasiou, D., Kumar, P., & Lagoudas, D. C. (2014).
 *   "A constitutive model for cyclic actuation of high-temperature shape memory alloys."
 *   *Mechanics of Materials*, 68, 120-136.
 */
void umat_sma_unified_T(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt);

/** @} */ // end of umat_mechanical group

} //namespace simcoon
