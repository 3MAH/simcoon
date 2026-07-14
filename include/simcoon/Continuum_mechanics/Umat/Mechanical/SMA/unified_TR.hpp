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
 * \f$ \xi \f$, the transformation strain \f$ \boldsymbol{\varepsilon}^{tr} \f$, the
 * reorientation strain accumulator \f$ \mathbf{E}^{Reo} \f$, and three
 * back-strain channels \f$ \mathbf{a}_F, \mathbf{a}_R, \mathbf{a}_{reo} \f$ that
 * grow with the forward / reverse / reorientation increments respectively.
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
 * \mathbf{X} = H^{Reo} \, (\mathbf{a} \odot \mathbf{I}_{r05}), \qquad
 * \mathbf{a} = \mathbf{a}_F + \mathbf{a}_R + \mathbf{a}_{reo}
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
 * | 0..16 | Identical to umat_sma_unified_T (T_init, xi, ET(6), xiF, xiR, calibrated rho-Ds0, rho-DE0, D, a1, a2, a3, Y0t) |
 * | 17 | \f$ p^{TR} \f$ : cumulative reorientation multiplier \f$ \sum \Delta s^{Reo} \f$ |
 * | 18..23 | \f$ \mathbf{v}^{re} \f$ : reorientation back-strain (Voigt, strain convention with engineering shears) |
 * | 24..29 | \f$ \mathbf{E}^{Reo} \f$ : macroscopic reorientation strain accumulator |
 *
 * @note `SMRDI`/`SMRDC`/`SMRAI`/`SMRAC` route through @ref umat_sma_unified_TR.
 * @note Internal vectorial state variables are rotated as **strain** by `DR`.
 * @note When \f$ \Delta s^{Reo} = 0 \f$, the result is bit-identical to
 * @ref umat_sma_unified_T (same elastic mixing, same forward/reverse equations).
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
