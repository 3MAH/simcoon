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
 * **Property layout for SMRDI (28 + 7 = 35 props):**
 *
 * | Index | Symbol | Description |
 * |-------|--------|-------------|
 * | 0..27 | -- | Same as SMADI (flagT, EA, EM, nuA, nuM, alphaA, alphaM, Hmin, Hmax, k1, sigmacrit, C_A, C_M, Ms0, Mf0, As0, Af0, n1..n4, sigmacaliber, prager_b, prager_n, c_lambda, p0_lambda, n_lambda, alpha_lambda) |
 * | 28 | \f$ Y^{Reo} \f$ | Stress limit for onset of reorientation |
 * | 29 | \f$ H^{Reo} \f$ | Reorientation kinematic hardening coefficient |
 * | 30 | \f$ E_{T}^{Reo,max} \f$ | Maximum reorientation back-strain magnitude |
 * | 31 | \f$ c_{\lambda Reo} \f$ | Reorientation penalty exponent start point |
 * | 32 | \f$ p_{0\lambda Reo} \f$ | Reorientation penalty limit value |
 * | 33 | \f$ n_{\lambda Reo} \f$ | Reorientation penalty power exponent |
 * | 34 | \f$ \alpha_{\lambda Reo} \f$ | Reorientation penalty power parameter |
 *
 * For SMRDC, SMRAI, SMRAC the same 7 reorientation parameters are appended
 * after the SMADC / SMAAI / SMAAC property block.
 *
 * **State variables (statev) — 30 entries:**
 *
 * Per Chatziathanasiou §4.3, the internal variables of the model are
 * \f$ \{\boldsymbol\sigma, T, \boldsymbol\varepsilon^F, \boldsymbol\varepsilon^R, \boldsymbol\varepsilon^{re}, \mathbf{v}^{re}, \xi^F, \xi^R\} \f$.
 * There is only one back-strain in the entire model — \f$ \mathbf{v}^{re} \f$ —
 * associated with the reorientation channel. Forward and reverse transformation
 * use isotropic hardening only and do not carry back-strains.
 *
 * | Index | Description |
 * |-------|-------------|
 * | 0 | \f$ T_{init} \f$ : initial temperature |
 * | 1 | \f$ \xi \f$ : martensitic volume fraction |
 * | 2..7 | \f$ \boldsymbol{\varepsilon}^{tr} \f$ : transformation strain |
 * | 8 | \f$ \xi^F \f$ : cumulated forward transformation, non-decreasing (\f$ \xi = \xi^F - \xi^R \f$) |
 * | 9 | \f$ \xi^R \f$ : cumulated reverse transformation, non-decreasing, kept apart from \f$ \xi^F \f$ so that both can be active in the same increment |
 * | 10..16 | calibrated constants, as in umat_sma_unified_T: \f$ \rho\Delta s_0, \rho\Delta E_0, D, a_1, a_2, a_3, Y_0^t \f$ |
 * | 17 | \f$ p^{TR} \f$ : cumulative reorientation multiplier \f$ \sum \Delta s^{Reo} \f$ |
 * | 18..23 | \f$ \mathbf{v}^{re} \f$ : reorientation back-strain (Voigt, strain convention with engineering shears) |
 * | 24..29 | \f$ \mathbf{E}^{Reo} \f$ : macroscopic reorientation strain accumulator |
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
 * @note When \f$ \Delta s^{Reo} = 0 \f$, the result is bit-identical to
 * @ref umat_sma_unified_T (same elastic mixing, same forward/reverse equations).
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
 * @note **Energy.** \f$ W_m^{ir} = 0 \f$: the back-strain is a kinematic history variable that
 * shifts the reorientation criterion and stores no free energy (no \f$ \mathbf{X}:\Delta\mathbf{v} \f$
 * term); the reorientation strain increment enters \f$ W_m^{r} \f$ and \f$ W_m^{d} \f$ like the
 * transformation strain.
 * @note `tangent_algorithmic` is clamped to the continuum tangent (see @ref umat_sma_unified_T).
 *
 * @see umat_sma_unified_T for the transformation-only base model.
 * @see assemble_continuum_tangent for the shared 3-mechanism tangent operator.
 *
 * **References:**
 * - Chatziathanasiou, D. (2016). Ph.D. Thesis.
 * - Chatziathanasiou, Chemisky, Meraghni, Echchorfi, Patoor (2015).
 *   *Smart Materials and Structures*.
 */
void umat_sma_unified_TR(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &stress, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode = tangent_default);

/** @} */

} //namespace simcoon
