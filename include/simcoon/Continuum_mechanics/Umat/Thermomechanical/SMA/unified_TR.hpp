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
 * @file unified_TR.hpp
 * @brief Thermomechanical face of the unified SMA model with martensite reorientation (SMR*)
 */

#pragma once
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon{

/** @addtogroup umat_thermomechanical
 *  @{
 */

/**
 * @brief Unified SMA model with transformation and reorientation, coupled to the heat equation
 *        (SMRDI, SMRDC, SMRAI, SMRAC).
 *
 * @details Fully coupled thermomechanical version of umat_sma_unified_TR: transformation (forward and
 * reverse) and martensite reorientation, with the heat equation. The mechanical model (criteria,
 * flows, back-strain \f$ \mathbf{v}^{re} \f$, both integrators, active set, closest-point
 * construction) is the one of umat_sma_unified_TR; this kernel adds the density and the specific
 * heats and returns the thermal tangent operators.
 *
 * **Elastic Symmetry and Criteria Selection:**
 *
 * | umat_name | Elasticity | Criteria | Number of props |
 * |-----------|------------|----------|-----------------|
 * | SMRDI | Isotropic | Drucker | 38 |
 * | SMRDC | Cubic | Drucker | 40 |
 * | SMRAI | Isotropic | Anisotropic Drucker (DFA) | 45 |
 * | SMRAC | Cubic | Anisotropic Drucker (DFA) | 47 |
 *
 * **Material Parameters (props):**
 *
 * **SMRDI (Isotropic elasticity, Drucker criteria) — 38 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ \rho \f$ | Density | kg/m^3 |
 * | 1 | \f$ c_{pA} \f$ | Specific heat capacity of Austenite | J/(kg.K) |
 * | 2 | \f$ c_{pM} \f$ | Specific heat capacity of Martensite | J/(kg.K) |
 * | 3 | flagT | 0: the given transformation temperatures are those of the measured curve, converted to the model's with n1..n4; 1: used as given | - |
 * | 4 | \f$ E_A \f$ | Young's modulus of Austenite | MPa |
 * | 5 | \f$ E_M \f$ | Young's modulus of Martensite | MPa |
 * | 6 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 7 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 8 | \f$ \alpha_A \f$ | CTE of Austenite | 1/K |
 * | 9 | \f$ \alpha_M \f$ | CTE of Martensite | 1/K |
 * | 10 | \f$ H_{min} \f$ | Minimal transformation strain magnitude | - |
 * | 11 | \f$ H_{max} \f$ | Maximal transformation strain magnitude | - |
 * | 12 | \f$ k_1 \f$ | Exponential evolution parameter for \f$ H^{cur} \f$ | 1/MPa |
 * | 13 | \f$ \sigma_{crit} \f$ | Von Mises stress below which \f$ H^{cur} = H_{min} \f$ | MPa |
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
 * | 31 | \f$ Y^{Reo} \f$ | Stress limit for onset of reorientation | MPa |
 * | 32 | \f$ H^{Reo} \f$ | Reorientation kinematic hardening coefficient | MPa |
 * | 33 | \f$ E_{T}^{Reo,max} \f$ | Maximum reorientation back-strain magnitude | - |
 * | 34 | \f$ c_{\lambda Reo} \f$ | Reorientation penalty exponent start point | - |
 * | 35 | \f$ p_{0\lambda Reo} \f$ | Reorientation penalty limit value | - |
 * | 36 | \f$ n_{\lambda Reo} \f$ | Reorientation penalty power exponent | - |
 * | 37 | \f$ \alpha_{\lambda Reo} \f$ | Reorientation penalty power parameter | - |
 *
 * **SMRDC (Cubic elasticity, Drucker criteria) — 40 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0–3 | | Same as SMRDI props[0–3] (\f$ \rho, c_{pA}, c_{pM} \f$, flagT) | |
 * | 4 | \f$ E_A \f$ | Young's modulus of Austenite ([100] direction) | MPa |
 * | 5 | \f$ E_M \f$ | Young's modulus of Martensite ([100] direction) | MPa |
 * | 6 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 7 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 8 | \f$ G_A \f$ | Shear modulus of Austenite | MPa |
 * | 9 | \f$ G_M \f$ | Shear modulus of Martensite | MPa |
 * | 10–32 | | Same as SMRDI props[8–30] (common SMA parameters) | |
 * | 33–39 | | Same as SMRDI props[31–37] (reorientation parameters) | |
 *
 * **SMRAI (Isotropic elasticity, anisotropic Drucker criteria) — 45 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0–30 | | Same as SMRDI props[0–30] | |
 * | 31 | \f$ F_{dfa} \f$ | F parameter of DFA criteria | - |
 * | 32 | \f$ G_{dfa} \f$ | G parameter of DFA criteria | - |
 * | 33 | \f$ H_{dfa} \f$ | H parameter of DFA criteria | - |
 * | 34 | \f$ L_{dfa} \f$ | L parameter of DFA criteria | - |
 * | 35 | \f$ M_{dfa} \f$ | M parameter of DFA criteria | - |
 * | 36 | \f$ N_{dfa} \f$ | N parameter of DFA criteria | - |
 * | 37 | \f$ K_{dfa} \f$ | K parameter of DFA criteria | - |
 * | 38–44 | | Same as SMRDI props[31–37] (reorientation parameters) | |
 *
 * **SMRAC (Cubic elasticity, anisotropic Drucker criteria) — 47 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0–32 | | Same as SMRDC props[0–32] | |
 * | 33–39 | | DFA parameters \f$ F, G, H, L, M, N, K \f$ (as SMRAI props[31–37]) | - |
 * | 40–46 | | Same as SMRDI props[31–37] (reorientation parameters) | |
 *
 * **State Variables (statev) — 30 variables:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ T_{init} \f$ | Initial/reference temperature | K |
 * | 1 | \f$ \xi \f$ | Martensitic volume fraction | - |
 * | 2..7 | \f$ \boldsymbol{\varepsilon}^{tr} \f$ | Transformation strain (Voigt, engineering shears) | - |
 * | 8 | \f$ \xi^F \f$ | Cumulated forward transformation, non-decreasing (\f$ \xi = \xi^F - \xi^R \f$) | - |
 * | 9 | \f$ \xi^R \f$ | Cumulated reverse transformation, non-decreasing | - |
 * | 10 | \f$ \rho \Delta s_0 \f$ | Entropy difference between phases (M - A) | MPa/K |
 * | 11 | \f$ \rho \Delta E_0 \f$ | Internal energy difference between phases (M - A) | MPa |
 * | 12 | \f$ D \f$ | Stress dependence parameter for transformation limits | - |
 * | 13 | \f$ a_1 \f$ | Forward hardening parameter | MPa |
 * | 14 | \f$ a_2 \f$ | Reverse hardening parameter | MPa |
 * | 15 | \f$ a_3 \f$ | Equilibrium hardening parameter | MPa |
 * | 16 | \f$ Y_{0t} \f$ | Initial transformation critical value | MPa |
 * | 17 | \f$ p^{TR} \f$ | Cumulative reorientation multiplier \f$ \sum \Delta s^{Reo} \f$ | - |
 * | 18..23 | \f$ \mathbf{v}^{re} \f$ | Reorientation back-strain (Voigt, engineering shears) | - |
 * | 24..29 | \f$ \mathbf{E}^{Reo} \f$ | Macroscopic reorientation strain accumulator (Voigt, engineering shears) | - |
 *
 * @param umat_name Model variant (SMRDI, SMRDC, SMRAI, SMRAC)
 * @param Etot Total strain tensor at beginning of increment (Voigt notation: 6x1)
 * @param DEtot Strain increment tensor (Voigt notation: 6x1)
 * @param sigma Stress tensor (Voigt notation: 6x1), see the stress-measure note of umat_sma_unified_TR [input/output]
 * @param r Heat of the increment per reference volume, storage included (heat_source_increment()) [output]
 * @param dSdE Mechanical tangent modulus \f$ \partial \sigma / \partial \varepsilon \f$ (6x6) [output]
 * @param dSdT Thermal stress tangent \f$ \partial \sigma / \partial T \f$ (6) [output]
 * @param drdE Thermal-mechanical tangent \f$ \partial r / \partial \varepsilon \f$ (6) [output]
 * @param drdT Thermal tangent \f$ \partial r / \partial T \f$ (scalar) [output]
 * @param DR Rotation increment matrix (3x3) for objective integration
 * @param nprops Number of material properties
 * @param props Material properties vector (layout depends on umat_name)
 * @param nstatev Number of state variables (30)
 * @param statev State variables vector [input/output]
 * @param T Temperature at beginning of increment
 * @param DT Temperature increment
 * @param Time Time at beginning of increment
 * @param DTime Time increment
 * @param Wm Total mechanical work [output]
 * @param Wm_r Recoverable work: elastic, transformation forces and back-strain force [output]
 * @param Wm_ir Irrecoverable work (always 0) [output]
 * @param Wm_d Dissipated work (hysteresis) [output]
 * @param Wt Total thermal work [output]
 * @param Wt_r Recoverable thermal work [output]
 * @param Wt_ir Irrecoverable thermal work [output]
 * @param ndi Number of direct stress components (typically 3)
 * @param nshr Number of shear stress components (typically 3)
 * @param start Flag indicating first increment (true) or continuation (false)
 * @param tnew_dt Suggested new time step ratio for adaptive time stepping [output]
 * @param tangent_mode kernel tangent code (parameter.hpp): 0 none, 1 continuum, 2 algorithmic,
 *        3 to 5 the closest-point integrator with the algorithmic, continuum or elastic operator
 * @param ncp NCPMethod code of the local step (0 Fischer-Burmeister, the default; see sma_ncp_method)
 *
 * @note Elastic convention: isotropic uses "Enu", cubic uses "EnuG" for L_iso / L_cubic
 * @note flagT = 0 converts the given Ms, Mf, As, Af (read on the measured curve) into the model's
 *       temperatures through the smoothing exponents n1..n4; flagT = 1 uses them as given
 * @note Convergence measure, reverse row of the local Jacobian, active set, reorientation flows and
 * closest-point construction: as @ref umat_sma_unified_TR (see its notes).
 * @note **Tangent modes.** As @ref umat_sma_unified_T_T: under the cutting plane `tangent_none`,
 * `tangent_algorithmic` and `tangent_continuum` all give the continuum tangent; the closest point
 * (codes 3 to 5) takes the operator of its code on the converged solve, the elastic choice taking the
 * continuum sensitivities; with the algorithmic operator the thermal face is the exact derivative of
 * the update.
 * @note **Thermal face.** The one of @ref umat_sma_unified_T_T with the reorientation row as a third
 * mechanism (\f$ N = 3 \f$). With \f$ \boldsymbol{\alpha}_q \f$ the thermal expansion of the elastic
 * prediction and \f$ \tilde{\mathbf{L}} \f$, \f$ \tilde{\boldsymbol{\kappa}} \f$ those of the operator,
 * \f[
 * \mathbf{P}_\theta = \hat{\mathbf{B}}^{-1}\left(\frac{\partial \Phi_l}{\partial \theta}
 * - \frac{\partial \Phi_l}{\partial \boldsymbol{\sigma}} : \tilde{\mathbf{L}}\,\boldsymbol{\alpha}_q\right)_{l = F, R, Reo}, \qquad
 * \frac{\partial \boldsymbol{\sigma}}{\partial T} = -\tilde{\mathbf{L}}\,\boldsymbol{\alpha}_q - \sum_j \tilde{\boldsymbol{\kappa}}_j\,P_{\theta,j},
 * \f]
 * where the reorientation criterion does not depend on the temperature
 * (\f$ \partial\Phi^{Reo}/\partial\theta = 0 \f$). The back-strain channel carries the force
 * \f$ \mathbf{A}_{Reo} = -(1 + \lambda_1^{Reo})\,\mathbf{X} \f$, and the dissipation of the increment is
 * \f[
 * \Delta\gamma = \tfrac12(\boldsymbol{\sigma}_n + \boldsymbol{\sigma}):\Delta\boldsymbol{\varepsilon}^{T}
 * + \sum_{j = F, R} \tfrac12 (A_{j,n} + A_j)\,\Delta\xi_j
 * + \tfrac12(\mathbf{A}_{Reo,n} + \mathbf{A}_{Reo}):\Delta\mathbf{E}^{Reo},
 * \f]
 * \f$ \Delta\boldsymbol{\varepsilon}^{T} \f$ the increment of the transformation strain, reorientation
 * included: the force terms are stored (Wm_r, as on the mechanical face) and the rest heats. The
 * reorientation does not enter the entropy,
 * \f$ \eta = c_0(\xi)\ln(T/T_{init}) + \boldsymbol{\alpha}(\xi_n):\boldsymbol{\sigma} + \rho\Delta s_0\,\xi \f$;
 * heat_source_increment() forms \f$ r \f$, \f$ \partial r/\partial\boldsymbol{\varepsilon} \f$,
 * \f$ \partial r/\partial T \f$.
 *
 * @see umat_sma_unified_TR() for the mechanical-only version (no thermal coupling)
 * @see umat_sma_unified_T_T() for the transformation-only thermomechanical model
 * @see L_iso() for isotropic stiffness tensor (SMRDI, SMRAI)
 * @see L_cubic() for cubic stiffness tensor (SMRDC, SMRAC)
 *
 * **Complementarity function**: @p ncp drives the local step of the closest point; the cutting
 * plane takes Fischer-Burmeister and refuses @c min / @c active_set (sma_ncp_method).
 *
 * **References:**
 * - Chatziathanasiou, D. (2016). Ph.D. Thesis.
 * - Chatziathanasiou, Chemisky, Meraghni, Echchorfi, Patoor (2015).
 *   *Smart Materials and Structures*.
 */
void umat_sma_unified_TR_T(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, double &r, arma::mat &dSdE, arma::vec &dSdT, arma::vec &drdE, double &drdT, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default, const int &ncp = 0);

/** @} */ // end of umat_thermomechanical group

} //namespace simcoon
