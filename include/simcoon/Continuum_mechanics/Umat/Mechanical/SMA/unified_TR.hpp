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
 * @brief Unified SMA model with transformation AND martensite reorientation
 * @author Y. Chemisky, D. Chatziathanasiou
 *
 * Adds a third leading mechanism (reorientation) to umat_sma_unified_T.
 * Implemented in 1D-2D-3D. Based on Chatziathanasiou Ph.D Thesis (2016).
 */

#pragma once
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon {

/** @addtogroup umat_mechanical
 *  @{
 */

/**
 * @brief Unified phenomenological SMA model with phase transformation and reorientation
 *
 * @details Extends umat_sma_unified_T with a third leading mechanism: martensite
 * reorientation. The macroscopic state is the martensitic volume fraction
 * \f$ \xi = \xi^F - \xi^R \f$, where the cumulated forward and reverse transformations
 * \f$ \xi^F, \xi^R \f$ (non-decreasing: their increments are the non-negative multipliers of the
 * two mechanisms) are kept apart so that both can be active in the same increment, the
 * transformation strain \f$ \boldsymbol{\varepsilon}^{tr} \f$, the reorientation strain
 * accumulator \f$ \mathbf{E}^{Reo} \f$ and its cumulative multiplier \f$ p^{TR} \f$, and a
 * single back-strain \f$ \mathbf{v}^{re} \f$, which grows with the reorientation increments only.
 *
 * **Variants:**
 *
 * | umat_name | Elasticity | Criteria | nprops |
 * |-----------|------------|----------|--------|
 * | SMRDI | Isotropic | Drucker | 35 |
 * | SMRDC | Cubic | Drucker | 37 |
 * | SMRAI | Isotropic | Anisotropic Drucker (DFA) | 42 |
 * | SMRAC | Cubic | Anisotropic Drucker (DFA) | 44 |
 *
 * **Effective back-stress and saturation-coupled reorientation surface (Form B):**
 *
 * \f[
 * \mathbf{X} = H^{Reo} \, (\mathbf{v}^{re} \odot \mathbf{I}_{r05})
 * \f]
 * \f[
 * \Phi^{Reo} = \mathrm{Prager}\!\bigl( \boldsymbol{\sigma} - (1+\lambda_{1}^{Reo}) \, \mathbf{X} \bigr) - Y^{Reo}
 * \f]
 *
 * where \f$ \lambda_{1}^{Reo} = \mathrm{lagrange}_{\mathrm{pow}_1}(\|\mathbf{v}^{re}\|/E_{T}^{Reo,max}, \ldots) \f$
 * is the saturation penalty on the intrinsic (per-unit-martensite) back-strain
 * \f$ \mathbf{v}^{re} \f$.
 *
 * @note **Stress measure.** The stress returned by this model (the `stress` argument,
 * written \f$ \boldsymbol{\sigma} \f$ in the relations above) is the Cauchy stress under
 * infinitesimal strain; under finite strain the update runs in a corotational frame, so it
 * is the rotated Kirchhoff stress
 * \f$ \hat{\boldsymbol{\tau}} = \boldsymbol{Q}^{T}\boldsymbol{\tau}\,\boldsymbol{Q} \f$ on the
 * frame fixed by the chosen objective rate (\f$ \boldsymbol{Q} = \boldsymbol{R} \f$ for
 * Green--Naghdi and \f$ \log_R \f$, the logarithmic frame for the XBM/log rate,
 * \f$ \boldsymbol{F} \f$ for \f$ \log_F \f$).
 *
 * **Material Parameters (props):**
 *
 * Each variant takes the props of its SMA* counterpart of umat_sma_unified_T (same order, same
 * meaning) followed by the 7 reorientation parameters.
 *
 * **SMRDI (Isotropic elasticity, Drucker criteria) — 35 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | flagT | 0: the given transformation temperatures are those of the measured curve, converted to the model's with n1..n4; 1: used as given | - |
 * | 1 | \f$ E_A \f$ | Young's modulus of Austenite | MPa |
 * | 2 | \f$ E_M \f$ | Young's modulus of Martensite | MPa |
 * | 3 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 4 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 5 | \f$ \alpha_A \f$ | CTE of Austenite | 1/K |
 * | 6 | \f$ \alpha_M \f$ | CTE of Martensite | 1/K |
 * | 7 | \f$ H_{min} \f$ | Minimal transformation strain magnitude | - |
 * | 8 | \f$ H_{max} \f$ | Maximal transformation strain magnitude | - |
 * | 9 | \f$ k_1 \f$ | Exponential evolution parameter for \f$ H^{cur} \f$ | 1/MPa |
 * | 10 | \f$ \sigma_{crit} \f$ | Von Mises stress below which \f$ H^{cur} = H_{min} \f$ | MPa |
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
 * | 28 | \f$ Y^{Reo} \f$ | Stress limit for onset of reorientation | MPa |
 * | 29 | \f$ H^{Reo} \f$ | Reorientation kinematic hardening coefficient | MPa |
 * | 30 | \f$ E_{T}^{Reo,max} \f$ | Maximum reorientation back-strain magnitude | - |
 * | 31 | \f$ c_{\lambda Reo} \f$ | Reorientation penalty exponent start point | - |
 * | 32 | \f$ p_{0\lambda Reo} \f$ | Reorientation penalty limit value | - |
 * | 33 | \f$ n_{\lambda Reo} \f$ | Reorientation penalty power exponent | - |
 * | 34 | \f$ \alpha_{\lambda Reo} \f$ | Reorientation penalty power parameter | - |
 *
 * **SMRDC (Cubic elasticity, Drucker criteria) — 37 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | flagT | as SMRDI | - |
 * | 1 | \f$ E_A \f$ | Young's modulus of Austenite ([100] direction) | MPa |
 * | 2 | \f$ E_M \f$ | Young's modulus of Martensite ([100] direction) | MPa |
 * | 3 | \f$ \nu_A \f$ | Poisson's ratio of Austenite | - |
 * | 4 | \f$ \nu_M \f$ | Poisson's ratio of Martensite | - |
 * | 5 | \f$ G_A \f$ | Shear modulus of Austenite | MPa |
 * | 6 | \f$ G_M \f$ | Shear modulus of Martensite | MPa |
 * | 7–29 | | Same as SMRDI props[5–27] (common SMA parameters) | |
 * | 30–36 | | Same as SMRDI props[28–34] (reorientation parameters) | |
 *
 * **SMRAI (Isotropic elasticity, anisotropic Drucker criteria) — 42 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0–27 | | Same as SMRDI props[0–27] | |
 * | 28 | \f$ F_{dfa} \f$ | F parameter of DFA criteria | - |
 * | 29 | \f$ G_{dfa} \f$ | G parameter of DFA criteria | - |
 * | 30 | \f$ H_{dfa} \f$ | H parameter of DFA criteria | - |
 * | 31 | \f$ L_{dfa} \f$ | L parameter of DFA criteria | - |
 * | 32 | \f$ M_{dfa} \f$ | M parameter of DFA criteria | - |
 * | 33 | \f$ N_{dfa} \f$ | N parameter of DFA criteria | - |
 * | 34 | \f$ K_{dfa} \f$ | K parameter of DFA criteria | - |
 * | 35–41 | | Same as SMRDI props[28–34] (reorientation parameters) | |
 *
 * **SMRAC (Cubic elasticity, anisotropic Drucker criteria) — 44 props:**
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0–29 | | Same as SMRDC props[0–29] | |
 * | 30–36 | | DFA parameters \f$ F, G, H, L, M, N, K \f$ (as SMRAI props[28–34]) | - |
 * | 37–43 | | Same as SMRDI props[28–34] (reorientation parameters) | |
 *
 * **State variables (statev) — 30 entries:**
 *
 * Per Chatziathanasiou §4.3, the internal variables of the model are
 * \f$ \{\boldsymbol\sigma, T, \boldsymbol\varepsilon^F, \boldsymbol\varepsilon^R, \boldsymbol\varepsilon^{re}, \mathbf{v}^{re}, \xi^F, \xi^R\} \f$.
 * There is only one back-strain in the entire model — \f$ \mathbf{v}^{re} \f$ —
 * associated with the reorientation channel. Forward and reverse transformation
 * use isotropic hardening only and do not carry back-strains.
 *
 * | Index | Symbol | Description | Units |
 * |-------|--------|-------------|-------|
 * | 0 | \f$ T_{init} \f$ | Initial/reference temperature | K |
 * | 1 | \f$ \xi \f$ | Martensitic volume fraction | - |
 * | 2..7 | \f$ \boldsymbol{\varepsilon}^{tr} \f$ | Transformation strain (Voigt, engineering shears) | - |
 * | 8 | \f$ \xi^F \f$ | Cumulated forward transformation, non-decreasing (\f$ \xi = \xi^F - \xi^R \f$) | - |
 * | 9 | \f$ \xi^R \f$ | Cumulated reverse transformation, non-decreasing, kept apart from \f$ \xi^F \f$ so that both can be active in the same increment | - |
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
 * @param stress Stress tensor (Voigt notation: 6x1), see the stress-measure note [input/output]
 * @param Lt Consistent tangent modulus (6x6) [output]
 * @param L Elastic stiffness tensor (6x6) [output]
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
 *
 * @note `SMRDI`/`SMRDC`/`SMRAI`/`SMRAC` route through @ref umat_sma_unified_TR.
 * @note Internal vectorial state variables are rotated as **strain** by `DR`.
 * @note **Active set.** The local system couples the three mechanisms in one
 * Fischer-Burmeister solve. The reorientation row has purely kinematic hardening,
 * \f$ B_{22} \simeq -(1+\lambda_1^{Reo})\,H^{Reo}\,\boldsymbol{\eta}^{Reo}:\boldsymbol{\eta}^{Reo} \f$,
 * and at the apex of the Drucker cone (\f$ \boldsymbol{\sigma}_{eff} = 0 \f$: a thermal step at
 * zero stress; or \f$ \xi = 0 \f$ with \f$ H^{Reo} = 0 \f$) the direction
 * \f$ \boldsymbol{\eta}^{Reo} \f$ is undefined and both the stress-flow and the back-strain
 * terms of the row vanish, so \f$ B_{22} = 0 \f$ and the row
 * \f$ -|B_{22}|\,\Delta s^{Reo} = 0 \f$ is singular. The transformation rows keep the
 * isotropic term \f$ \partial \Phi / \partial \xi \f$ on their diagonal and never degenerate.
 * A mechanism with \f$ \Phi_j < 0 \f$ and \f$ B_{jj} = 0 \f$ is therefore dropped from the
 * iteration, \f$ \Delta s_j = 0 \f$, and the reduced system is solved: the row of an
 * inactive mechanism with a finite diagonal decouples and gives \f$ \Delta s_j = 0 \f$
 * exactly, so the reduced solve returns the same increments as the full one wherever the
 * full one is solvable.
 * @note Convergence measure: same magnitude scales \f$ Y^{crit}_j \f$ as
 * @ref umat_sma_unified_T (see its note); \f$ Y^{crit}_{Reo} = Y^{Reo} \f$.
 * @note When \f$ \Delta s^{Reo} = 0 \f$, the result agrees with @ref umat_sma_unified_T to
 * round-off under either integrator (same elastic mixing, same forward/reverse equations, the
 * same cutting-plane scheme; measured within 2e-12 MPa on a superelastic loop).
 * @note **Reorientation flows** (thesis §4, eq. 503). Reorientation only re-orients the
 * martensite present at the start of the increment, \f$ \xi_{start} \f$ (new martensite forms
 * already aligned with \f$ \boldsymbol{\sigma} \f$), so with
 * \f$ \bar{\boldsymbol{\Lambda}}^{re} = \partial\,\mathrm{Prager}/\partial\boldsymbol{\sigma}\,(\boldsymbol{\sigma}_{eff}) \f$:
 * \f[
 * \boldsymbol{\Lambda}_{\varepsilon}^{re} = \xi_{start}\,\bar{\boldsymbol{\Lambda}}^{re}, \qquad
 * \dot{\mathbf{v}}^{re} = \dot{p}^{re}\,\bar{\boldsymbol{\Lambda}}^{re}, \qquad
 * \frac{\partial \Phi^{Reo}}{\partial \boldsymbol{\sigma}} = \bar{\boldsymbol{\Lambda}}^{re}.
 * \f]
 * The criterion and the back-strain flow are per unit martensite (not \f$ \xi \f$-scaled), which
 * keeps \f$ B_{22} \f$ finite at \f$ \xi = 0 \f$ through the back-strain term of \f$ K_{22} \f$.
 * The saturation argument is the intrinsic ratio \f$ \|\mathbf{v}^{re}\| / E_T^{Reo,max} \f$:
 * dividing by \f$ \xi \f$ would blow the penalty up as \f$ \xi \to 0 \f$. Hence
 * \f$ K_{2F} = K_{2R} = 0 \f$ and
 * \f[
 * K_{22} = \frac{\partial \Phi^{Reo}}{\partial \mathbf{v}^{re}} : \bar{\boldsymbol{\Lambda}}^{re}, \qquad
 * \frac{\partial \Phi^{Reo}}{\partial \mathbf{v}^{re}} = -(1+\lambda_1^{Reo})\,H^{Reo}\,(\boldsymbol{\eta}^{Reo}\odot\mathbf{I}_{r05})
 * - (\boldsymbol{\eta}^{Reo}:\mathbf{X})\,\frac{\partial \lambda_1^{Reo}}{\partial \mathbf{v}^{re}}.
 * \f]
 * The reverse row uses the same ETMean chain rule as @ref umat_sma_unified_T, with
 * \f$ \boldsymbol{\Lambda}_{\bar{\varepsilon}}^{re} = \bar{\boldsymbol{\Lambda}}^{re} \f$
 * (\f$ \xi \f$ unchanged by reorientation).
 * @note **Closest-point mode** (`tangent_closest_point`): same frozen-mixture construction as
 * @ref umat_sma_unified_T, with \f$ \boldsymbol{\Lambda}^{re} = \xi_{start}\,\mathbf{r} \f$,
 * \f$ \mathbf{r} = \bar{\boldsymbol{\Lambda}}^{re}(\boldsymbol{\sigma}_{eff}) \f$. The back-strain is
 * solved by Newton on its backward-Euler equation
 * \f$ \mathbf{R}(\mathbf{v}) = \mathbf{v} - \mathbf{v}_n - \Delta s^{re}\,\mathbf{r}(\boldsymbol{\sigma}_{eff}(\mathbf{v})) = \mathbf{0} \f$,
 * Jacobian \f$ \mathbf{J} = \mathbf{I} - \Delta s^{re}\,\mathbf{G}\,\mathbf{A} \f$ with
 * \f$ \mathbf{G} = \partial\mathbf{r}/\partial\boldsymbol{\sigma}_{eff} \f$ and
 * \f$ \mathbf{A} = \partial\boldsymbol{\sigma}_{eff}/\partial\mathbf{v} \f$ (a fixed-point sweep only
 * contracts while \f$ \Delta s^{re} H^{Reo}(1+\lambda_1^{Reo})/\bar{\sigma}_{eff} < 1 \f$); no root
 * within the budget rejects the trial. Then
 * \f[
 * \frac{\partial\mathbf{v}}{\partial\boldsymbol{\sigma}} = \mathbf{J}^{-1}\Delta s^{re}\,\mathbf{G}, \quad
 * \frac{\partial\mathbf{v}}{\partial\Delta s^{re}} = \mathbf{J}^{-1}\mathbf{r}, \quad
 * \frac{\partial\boldsymbol{\sigma}_{eff}}{\partial\boldsymbol{\sigma}} = \mathbf{I} + \mathbf{A}\,\frac{\partial\mathbf{v}}{\partial\boldsymbol{\sigma}}, \quad
 * \frac{\partial\boldsymbol{\sigma}_{eff}}{\partial\Delta s^{re}} = \mathbf{A}\,\frac{\partial\mathbf{v}}{\partial\Delta s^{re}},
 * \f]
 * and \f$ \bar{\boldsymbol{\varepsilon}}^{T}(\xi_n + \Delta\xi_F) = \mathrm{dev}(\boldsymbol{\varepsilon}^{T}_n)
 * + \Delta\xi_F\,\mathrm{dev}(\boldsymbol{\Lambda}^F) + \Delta s^{re}\,\mathrm{dev}(\boldsymbol{\Lambda}^{re}) \f$.
 * @note **Energy.** \f$ W_m^{ir} = 0 \f$. The back-strain channel carries the force
 * \f$ \mathbf{A}_{Reo} = -(1+\lambda_1^{Reo})\,\mathbf{X} \f$ (force/variable pair, as the
 * kinematic hardening of EPCHA), and with \f$ A_{\xi F}, A_{\xi R} \f$ the transformation forces
 * (bound penalties of \f$ \xi \f$ included) the dissipation of the increment is
 * \f[
 * \Delta\gamma = \tfrac12(\boldsymbol{\sigma}_n + \boldsymbol{\sigma}):(\Delta\boldsymbol{\varepsilon}^{F} - \Delta\boldsymbol{\varepsilon}^{R} + \Delta\mathbf{E}^{Reo})
 * + \tfrac12(A_{\xi F,n} + A_{\xi F})\,\Delta\xi^F + \tfrac12(A_{\xi R,n} + A_{\xi R})\,\Delta\xi^R
 * + \tfrac12(\mathbf{A}_{Reo,n} + \mathbf{A}_{Reo}):\Delta\mathbf{E}^{Reo},
 * \f]
 * \f$ W_m^{d} \mathrel{+}= \Delta\gamma \f$ and \f$ W_m^{r} = W_m - W_m^{d} \f$: the force terms are
 * stored and recovered at unloading.
 * @note `tangent_algorithmic` is clamped to the continuum tangent (see @ref umat_sma_unified_T).
 *
 * @see umat_sma_unified_T for the transformation-only base model.
 * @see assemble_continuum_tangent for the shared 3-mechanism tangent operator.
 *
 * **Complementarity function**: @p ncp (NCPMethod code, 0 Fischer-Burmeister by default) drives
 * the local step of the closest point (ReturnMappingControl::ncp, UMAT catalogue), whose merit is
 * then the one of ncp_step; the cutting plane takes Fischer-Burmeister and refuses
 * @c min / @c active_set (sma_ncp_method).
 *
 * **References:**
 * - Chatziathanasiou, D. (2016). Ph.D. Thesis.
 * - Chatziathanasiou, Chemisky, Meraghni, Echchorfi, Patoor (2015).
 *   *Smart Materials and Structures*.
 */
void umat_sma_unified_TR(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &stress, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default, const int &ncp = 0);

/** @} */

} //namespace simcoon
