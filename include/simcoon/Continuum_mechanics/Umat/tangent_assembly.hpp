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
    arma::mat P_epsilon;                ///< Strain sensitivities \f$ \mathbf{P}_\varepsilon^l \f$, one column per mechanism (\f$ 6\times N_{\mathrm{mech}} \f$).
    arma::mat invBhat;                  ///< Active-set-masked inverse \f$ \hat{B}^{-1} \f$
                                        ///< (\f$ N_{\mathrm{mech}}\times N_{\mathrm{mech}} \f$), exposed for
                                        ///< thermomechanical UMATs that build \f$ \mathbf{P}_\theta \f$ downstream.
    arma::mat::fixed<6, 6> L_tilde;     ///< Elastic operator the sensitivities are taken on: \f$ \mathbf{L} \f$ for the
                                        ///< continuum assembly, \f$ \tilde{\mathbf{L}} = \mathbf{M}^{-1}\mathbf{L} \f$ for
                                        ///< the algorithmic one, so a trial stress rate \f$ \mathbf{q} \f$ of another
                                        ///< origin (thermal) enters the rows as \f$ \tilde{\mathbf{L}}\mathbf{L}^{-1}\mathbf{q} \f$
                                        ///< and its flux as \f$ \tilde{\boldsymbol{\kappa}} \f$ below.
    arma::mat kappa_tilde;              ///< \f$ \tilde{\boldsymbol{\kappa}} = \mathbf{M}^{-1}\boldsymbol{\kappa} \f$ of the
                                        ///< algorithmic assembly (\f$ 6\times N_{\mathrm{mech}} \f$); empty when it is
                                        ///< \f$ \boldsymbol{\kappa} \f$ itself (continuum assembly, no correction).
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
 * this routine eliminates the multipliers from the linearised consistency conditions
 * \f$ \sum_j \hat{B}^{lj}\, \mathrm{d}s^j = \frac{\partial \Phi^l}{\partial \boldsymbol{\sigma}} : \mathbf{L} : \mathrm{d}\boldsymbol{\varepsilon} \f$
 * (criterion on the row, mechanism on the column) and computes
 * \f[
 *   \mathbf{P}_\varepsilon^l = \sum_m \bigl(\hat{B}^{-1}\bigr)_{lm}\,
 *                              \mathbf{L}\,\frac{\partial \Phi^m}{\partial \boldsymbol{\sigma}},
 *   \qquad
 *   \mathbf{L}_t = \mathbf{L} - \sum_l \boldsymbol{\kappa}^l \otimes \mathbf{P}_\varepsilon^l ,
 * \f]
 * so that \f$ \mathrm{d}s^l = \mathbf{P}_\varepsilon^l : \mathrm{d}\boldsymbol{\varepsilon} \f$.
 * \f$ \hat{B} \f$ is not symmetric as soon as two coupled mechanisms are active together
 * (reverse transformation and reorientation in the SMR* laws): the index order matters.
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
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L);

/**
 * @brief Single-mechanism convenience overload.
 *
 * The N-mechanism form with the scalar `Bhat`/`Ds` lifted to \f$1 \times 1\f$ / size-1
 * quantities, without the \f$ 1 \times 1 \f$ inverse.
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
 * expose the corresponding data (doc §7.4). Its consequence: the operator is the exact
 * derivative of the update only when the flow direction does not rotate during the step
 * (von Mises with isotropic hardening); with kinematic hardening or a non-radial criterion
 * it is an approximation (see the scope note in docs/simulation/solver.rst).
 *
 * @param Bhat_continuum   **Continuum** local Jacobian
 *   \f$ \hat{B}^{lj} = \partial \Phi^l/\partial \boldsymbol{\sigma} : \boldsymbol{\kappa}^j - K^{lj} \f$,
 *   i.e. the same quantity passed to @ref assemble_continuum_tangent. Do **not** pass an
 *   already-corrected algorithmic \f$ \tilde{B} \f$ — the routine builds \f$ \tilde{B} \f$
 *   internally and would double-correct, silently producing wrong tangents.
 * @param kappa_j          Fluxes \f$ \boldsymbol{\kappa}^j = \mathbf{L}\,\boldsymbol{\Lambda}_\varepsilon^j \f$, one column per mechanism (\f$ 6\times N \f$).
 * @param dPhidsigma_l     Criterion gradients \f$ \partial \Phi^l/\partial \boldsymbol{\sigma} \f$, one column per criterion (\f$ 6\times N \f$).
 * @param Ds_j             Per-mechanism increment \f$ \Delta s^j \f$ over the step.
 * @param L                Elastic stiffness \f$ \mathbf{L} \f$ (\f$ 6\times 6 \f$).
 * @param dLambda_dsigma_l \f$ \partial \boldsymbol{\Lambda}_\varepsilon^j/\partial \boldsymbol{\sigma} \f$ of every mechanism side by
 *   side, columns \f$ 6j..6j+5 \f$ for mechanism \f$ j \f$ (\f$ 6\times 6N \f$). Pass a zero block where the direction
 *   is independent of \f$ \boldsymbol{\sigma} \f$ (the algorithmic correction degenerates to the continuum tangent there).
 */
ContinuumTangent assemble_algorithmic_tangent(
    const arma::mat& Bhat_continuum,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
    const arma::vec& Ds_j,
    const arma::mat& L,
    const arma::mat& dLambda_dsigma_l);

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
 * @brief Lazy Hessian provider for compute_tangent_operator(): the \f$ 6\times 6N \f$ matrix of the
 * per-mechanism Hessians side by side (assemble_algorithmic_tangent()).
 *
 * Invoked ONLY when the algorithmic operator is requested — expensive Hessians
 * (e.g. the SMA finite-difference \f$ \partial \boldsymbol{\Lambda}/\partial
 * \boldsymbol{\sigma} \f$) must not be evaluated in the other modes.
 */
using HessianProvider = std::function<arma::mat()>;

// NOTE: only the iterative (plasticity/SMA/thermomechanical) kernels dispatch
// on tangent_mode. The finite-strain hyperelastic, linear viscoelastic, damage
// and elastic kernels always return their exact tangent and accept-but-ignore
// the mode except for tangent_none (documented in docs/simulation/umat_catalog.rst).

/**
 * @brief Mode dispatch over the shared tangent assemblies — the single entry
 * point UMATs call after their return mapping, keeping only the physical
 * inputs at the call site.
 *
 * - @c tangent_none: no assembly; returns { Lt = L, zero P_epsilon
 *   (\f$ 6\times N \f$), zero invBhat } so downstream thermomechanical algebra
 *   (P_theta from P_epsilon/invBhat) degrades consistently to the elastic
 *   operator (explicit integration).
 * - @c tangent_continuum: assemble_continuum_tangent().
 * - @c tangent_algorithmic: assemble_algorithmic_tangent() with
 *   @p dLambda_dsigma(); falls back to the continuum operator when no
 *   provider is given (flow independent of stress).
 * - @c tangent_closest_point: for a kernel whose local loop is the cutting-plane one, the
 *   algorithmic operator (same as @c tangent_algorithmic — the closest-point exact operator
 *   is only available from a closest_point_return_mapping() solve, see return_mapping.hpp,
 *   and such kernels take it from cpp_consistent_tangent() instead of calling this).
 *
 * @param tangent_mode One of the tangent_* constants (parameter.hpp)
 * @param Bhat Local Jacobian (N x N), as in assemble_continuum_tangent
 * @param kappa_j Flux directions (6 x N, one column per mechanism; a vec for N = 1)
 * @param dPhidsigma_l Criterion gradients (6 x N)
 * @param Ds_j Multiplier increments (active-set mask, N)
 * @param L Elastic stiffness (6 x 6)
 * @param dLambda_dsigma Lazy provider of the flow Hessians (6 x 6N)
 */
ContinuumTangent compute_tangent_operator(
    int tangent_mode,
    const arma::mat& Bhat,
    const arma::mat& kappa_j,
    const arma::mat& dPhidsigma_l,
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

/**
 * @brief Heat released over a thermomechanical increment and its exact derivatives.
 *
 * The thermomechanical UMATs return the heat of the increment, per unit reference volume,
 * \f[ R = \Delta t \, r = \Delta\gamma_{loc} - \theta_{n+1} \left( \eta_{n+1} - \eta_n \right), \f]
 * with \f$ \Delta\gamma_{loc} \f$ the dissipation the kernel adds to \f$ W_m^d \f$ and \f$ \eta \f$ its
 * entropy, a state function. \f$ R \f$ includes the thermal storage, \f$ R = 0 \f$ is the adiabatic
 * condition (also at \f$ \Delta t = 0 \f$, where it gives the instantaneous adiabatic temperature
 * jump) and a closed elastic cycle under \f$ R = 0 \f$ returns to its initial temperature. The
 * derivatives follow from those of \f$ \Delta\gamma_{loc} \f$ and \f$ \eta_{n+1} \f$:
 * \f[ \frac{\partial R}{\partial \boldsymbol{\varepsilon}} = \frac{\partial \Delta\gamma_{loc}}{\partial \boldsymbol{\varepsilon}}
 *     - \theta_{n+1} \frac{\partial \eta_{n+1}}{\partial \boldsymbol{\varepsilon}}, \qquad
 *     \frac{\partial R}{\partial \theta} = \frac{\partial \Delta\gamma_{loc}}{\partial \theta} - \Delta\eta
 *     - \theta_{n+1} \frac{\partial \eta_{n+1}}{\partial \theta}. \f]
 * For the Gibbs energies of the thermomechanical kernels,
 * \f$ \eta = c_0 \ln(\theta/\theta_0) + \boldsymbol{\alpha} : \boldsymbol{\sigma} + \ldots \f$, so that
 * \f$ \theta \, \partial\eta/\partial\theta \f$ at fixed strain is the heat capacity at constant strain
 * \f$ c_0 - \theta \, \boldsymbol{\alpha} : \mathbf{L} : \boldsymbol{\alpha} \f$ (with \f$ c_0 = \rho c_p \f$
 * the heat capacity at constant stress).
 *
 * @param Dgamma dissipation of the increment, \f$ \Delta\gamma_{loc} \f$
 * @param dDgamma_dE \f$ \partial\Delta\gamma_{loc}/\partial\boldsymbol{\varepsilon} \f$ (6)
 * @param dDgamma_dT \f$ \partial\Delta\gamma_{loc}/\partial\theta \f$
 * @param Deta entropy increment \f$ \eta_{n+1} - \eta_n \f$
 * @param deta_dE \f$ \partial\eta_{n+1}/\partial\boldsymbol{\varepsilon} \f$ (6)
 * @param deta_dT \f$ \partial\eta_{n+1}/\partial\theta \f$
 * @param T1 temperature at the end of the increment, \f$ \theta_{n+1} \f$
 * @param[out] r heat of the increment \f$ R \f$ (MPa, i.e. MJ/m\f$^3\f$)
 * @param[out] drdE \f$ \partial R/\partial\boldsymbol{\varepsilon} \f$ (6)
 * @param[out] drdT \f$ \partial R/\partial\theta \f$
 */
void heat_source_increment(const double &Dgamma, const arma::vec &dDgamma_dE, const double &dDgamma_dT,
                           const double &Deta, const arma::vec &deta_dE, const double &deta_dT,
                           const double &T1, double &r, arma::vec &drdE, double &drdT);

} // namespace simcoon
