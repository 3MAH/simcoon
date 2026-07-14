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

///@file tangent_assembly.hpp
///@brief Generic tangent-modulus assembly for dissipative UMATs.
///@version 1.0

#pragma once

#include <armadillo>
#include <vector>
#include <functional>

namespace simcoon {

/**
 * @brief Result of the continuum-tangent assembly.
 *
 * Bundles the assembled tangent operator and the per-mechanism strain
 * sensitivities, both of which some UMATs (notably the thermomechanical
 * mirrors and the SMA family) need downstream.
 */
struct ContinuumTangent {
    arma::mat Lt;                       ///< Tangent operator \f$ \mathbf{L}_t \f$ (\f$ 6\times 6 \f$).
    std::vector<arma::vec> P_epsilon;   ///< Per-mechanism strain sensitivity \f$ \mathbf{P}_\varepsilon^l \f$ (each 6).
    arma::mat invBhat;                  ///< Active-set-masked inverse \f$ \hat{B}^{-1} \f$
                                        ///< (\f$ N_{\mathrm{mech}}\times N_{\mathrm{mech}} \f$), exposed for
                                        ///< thermomechanical UMATs that build \f$ \mathbf{P}_\theta \f$ downstream.
};

/**
 * @brief Assemble the continuum elasto-(visco)plastic tangent following the
 * leading-mechanism framework of simcoon-documentation §7.4.
 *
 * For a UMAT with \f$ N_{\mathrm{mech}} \f$ leading mechanisms, given
 * - the local Jacobian
 *   \f$ \hat{B}^{lj} = \frac{\partial \Phi^l}{\partial \boldsymbol{\sigma}} : \boldsymbol{\kappa}^j - K^{lj} \f$,
 *   already built by the UMAT during the local Newton loop (@p Bhat),
 * - the flux directions \f$ \boldsymbol{\kappa}^j = \mathbf{L}\,\boldsymbol{\Lambda}_\varepsilon^j \f$ (@p kappa_j),
 * - the criterion gradients \f$ \partial \Phi^l / \partial \boldsymbol{\sigma} \f$ (@p dPhidsigma_l),
 * - the increments \f$ \Delta s^j \f$ over the time step (@p Ds_j), which drive the active set,
 *
 * this routine computes
 * \f[
 *   \mathbf{P}_\varepsilon^l = \sum_m \bigl(\hat{B}^{-1}\bigr)_{ml}\,
 *                              \mathbf{L}\,\frac{\partial \Phi^m}{\partial \boldsymbol{\sigma}},
 *   \qquad
 *   \mathbf{L}_t = \mathbf{L} - \sum_l \boldsymbol{\kappa}^l \otimes \mathbf{P}_\varepsilon^l .
 * \f]
 *
 * with the standard active-set masking: a mechanism \f$ j \f$ is treated as inactive when
 * \f$ \Delta s^j \le \f$ @c simcoon::iota; the corresponding row/column of \f$ \hat{B} \f$ is
 * replaced by the identity before inversion and masked back to zero, so an inactive
 * mechanism contributes nothing.
 *
 * This is the shared implementation of the tangent-assembly idiom: the iterative UMATs in
 * @c Mechanical/{Plasticity,SMA,Viscoelasticity,Damage}/ and their Thermomechanical mirrors
 * call it instead of re-deriving the operator inline.
 */
ContinuumTangent assemble_continuum_tangent(
    const arma::mat& Bhat,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L);

/**
 * @brief Single-mechanism convenience overload.
 *
 * Forwards to the N-mechanism form with `{kappa}`, `{dPhidsigma}`, and the
 * scalar `Bhat`/`Ds` lifted to \f$1 \times 1\f$ / size-1 quantities. Removes the
 * `{ vec }` brace-init dance from the 9 single-mechanism UMAT call sites.
 */
ContinuumTangent assemble_continuum_tangent(
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L);

/**
 * @brief Assemble the Simo–Hughes algorithmic (consistent) tangent following the same
 * leading-mechanism framework as @ref assemble_continuum_tangent, but accounting for the
 * variation of the strain-flow direction \f$ \boldsymbol{\Lambda}_\varepsilon^j \f$ with
 * stress during the increment.
 *
 * Given the converged-state derivative
 * \f$ \partial \boldsymbol{\Lambda}_\varepsilon^j / \partial \boldsymbol{\sigma} \f$
 * (a \f$ 6\times 6 \f$ matrix per mechanism, @p dLambda_dsigma_l), the discrete return map
 * \f[
 *   \boldsymbol{\sigma}_{n+1} = \mathbf{L}\,\Bigl( \boldsymbol{\varepsilon}_{n+1}
 *       - \boldsymbol{\varepsilon}^p_n
 *       - \sum_j \Delta s^j\, \boldsymbol{\Lambda}_\varepsilon^j(\boldsymbol{\sigma}_{n+1}, \ldots) \Bigr)
 * \f]
 * linearised at convergence yields the modified operator
 * \f[
 *   \mathbf{M} = \mathbf{I}_6 + \mathbf{L}\,\sum_j \Delta s^j\,
 *                \frac{\partial \boldsymbol{\Lambda}_\varepsilon^j}{\partial \boldsymbol{\sigma}},
 *   \quad
 *   \tilde{\mathbf{L}} = \mathbf{M}^{-1}\mathbf{L},
 *   \quad
 *   \tilde{\boldsymbol{\kappa}}^j = \mathbf{M}^{-1}\boldsymbol{\kappa}^j,
 * \f]
 * \f[
 *   \tilde{B}^{lj} = \frac{\partial \Phi^l}{\partial \boldsymbol{\sigma}} : \tilde{\boldsymbol{\kappa}}^j - K^{lj}
 *                  = \hat{B}^{lj} + \frac{\partial \Phi^l}{\partial \boldsymbol{\sigma}}
 *                    : (\mathbf{M}^{-1} - \mathbf{I})\,\boldsymbol{\kappa}^j ,
 * \f]
 * after which the assembly is structurally identical to the continuum helper, with
 * \f$ \tilde{\mathbf{L}} \f$ in place of \f$ \mathbf{L} \f$ and
 * \f$ \tilde{\boldsymbol{\kappa}} \f$ in place of \f$ \boldsymbol{\kappa} \f$:
 * \f[
 *   \mathbf{P}_\varepsilon^l = \sum_m \bigl(\tilde{B}^{-1}\bigr)_{ml}\,
 *                              \tilde{\mathbf{L}}\,\frac{\partial \Phi^m}{\partial \boldsymbol{\sigma}},
 *   \qquad
 *   \mathbf{L}_t = \tilde{\mathbf{L}} - \sum_l \tilde{\boldsymbol{\kappa}}^l \otimes \mathbf{P}_\varepsilon^l .
 * \f]
 * Active-set masking and \f$ \hat{B} \f$ regularisation reuse the same
 * \f$ \Delta s^j > \f$ @c simcoon::iota convention as the continuum helper.
 *
 * @note **Parity.** When every @p dLambda_dsigma_l is zero (the frozen-direction assumption
 * underlying the continuum tangent), \f$ \mathbf{M} = \mathbf{I} \f$,
 * \f$ \tilde{\mathbf{L}} = \mathbf{L} \f$,
 * \f$ \tilde{\boldsymbol{\kappa}} = \boldsymbol{\kappa} \f$, \f$ \tilde{B} = \hat{B} \f$, and
 * this routine returns results *bit-identical* to @ref assemble_continuum_tangent. That is
 * what lets the @c tangent_mode dispatch work in UMATs not yet extended with
 * \f$ \partial \boldsymbol{\Lambda}/\partial \boldsymbol{\sigma} \f$ data — the algorithmic
 * call is a strict superset of the continuum call.
 *
 * @note **State coupling.** This first version linearises
 * \f$ \boldsymbol{\Lambda}_\varepsilon^j \f$ w.r.t. \f$ \boldsymbol{\sigma} \f$ only. The fully
 * consistent tangent additionally requires
 * \f$ \partial \boldsymbol{\Lambda}_\varepsilon^j / \partial \mathbf{V} \f$ and
 * \f$ \boldsymbol{\Lambda}_V^{j,q} \f$; that extension will land as an overload once UMATs
 * expose the corresponding data (doc §7.4).
 *
 * @param Bhat_continuum   **Continuum** local Jacobian
 *   \f$ \hat{B}^{lj} = \partial \Phi^l/\partial \boldsymbol{\sigma} : \boldsymbol{\kappa}^j - K^{lj} \f$,
 *   i.e. the same quantity passed to @ref assemble_continuum_tangent. Do **not** pass an
 *   already-corrected algorithmic \f$ \tilde{B} \f$ — the routine builds \f$ \tilde{B} \f$
 *   internally and would double-correct, silently producing wrong tangents.
 * @param kappa_j          Per-mechanism flux \f$ \boldsymbol{\kappa}^j = \mathbf{L}\,\boldsymbol{\Lambda}_\varepsilon^j \f$ (each 6).
 * @param dPhidsigma_l     Per-mechanism criterion gradient \f$ \partial \Phi^l/\partial \boldsymbol{\sigma} \f$ (each 6).
 * @param Ds_j             Per-mechanism increment \f$ \Delta s^j \f$ over the step.
 * @param L                Elastic stiffness \f$ \mathbf{L} \f$ (\f$ 6\times 6 \f$).
 * @param dLambda_dsigma_l Per-mechanism \f$ \partial \boldsymbol{\Lambda}_\varepsilon^j/\partial \boldsymbol{\sigma} \f$
 *   (each \f$ 6\times 6 \f$). Pass a zero matrix where the direction is independent of
 *   \f$ \boldsymbol{\sigma} \f$ (the algorithmic correction degenerates to the continuum tangent there).
 */
ContinuumTangent assemble_algorithmic_tangent(
    const arma::mat& Bhat_continuum,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const std::vector<arma::mat>& dLambda_dsigma_l);

/**
 * @brief Single-mechanism convenience overload for the algorithmic tangent.
 */
ContinuumTangent assemble_algorithmic_tangent(
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L,
    const arma::mat& dLambda_dsigma);

/**
 * @brief Lazy per-mechanism Hessian provider for compute_tangent_operator().
 *
 * Invoked ONLY when the algorithmic operator is requested — expensive Hessians
 * (e.g. the SMA finite-difference $ \partial oldsymbol{\Lambda}/\partial
 * oldsymbol{\sigma} $) must not be evaluated in the other modes.
 */
using HessianProvider = std::function<std::vector<arma::mat>()>;

// NOTE: only the iterative (plasticity/SMA/thermomechanical) kernels dispatch
// on tangent_mode. The finite-strain hyperelastic, viscoelastic, damage and
// elastic kernels always return their exact tangent and accept-but-ignore the
// mode (documented in docs/simulation/umat_catalog.rst).

/**
 * @brief Mode dispatch over the shared tangent assemblies — the single entry
 * point UMATs call after their return mapping, keeping only the physical
 * inputs at the call site.
 *
 * - @c tangent_none: no assembly; returns { Lt = L, per-mechanism zero
 *   P_epsilon, zero invBhat } so downstream thermomechanical algebra
 *   (P_theta from P_epsilon/invBhat) degrades consistently to the elastic
 *   operator (explicit integration).
 * - @c tangent_continuum: assemble_continuum_tangent().
 * - @c tangent_algorithmic: assemble_algorithmic_tangent() with
 *   @p dLambda_dsigma(); falls back to the continuum operator when no
 *   provider is given (flow independent of stress).
 * - @c tangent_closest_point: reserved — throws std::invalid_argument.
 *
 * @param tangent_mode One of the tangent_* constants (parameter.hpp)
 * @param Bhat Local Jacobian (N x N), as in assemble_continuum_tangent
 * @param kappa_j Flux directions (N x 6-vector)
 * @param dPhidsigma_l Criterion gradients (N x 6-vector)
 * @param Ds_j Multiplier increments (active-set mask, N)
 * @param L Elastic stiffness (6 x 6)
 * @param dLambda_dsigma Lazy provider of the per-mechanism flow Hessians
 */
ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    const arma::mat& Bhat,
    const std::vector<arma::vec>& kappa_j,
    const std::vector<arma::vec>& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const HessianProvider& dLambda_dsigma = nullptr);

/**
 * @brief Single-mechanism convenience overload of compute_tangent_operator().
 * @param dLambda_dsigma Lazy provider of the single flow Hessian (6 x 6)
 */
ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L,
    const std::function<arma::mat()>& dLambda_dsigma = nullptr);

} // namespace simcoon
