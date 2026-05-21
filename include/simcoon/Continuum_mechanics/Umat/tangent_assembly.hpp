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

namespace simcoon {

/**
 * @brief Result of the continuum-tangent assembly.
 *
 * Bundles the assembled tangent operator and the per-mechanism strain
 * sensitivities, both of which some UMATs (notably the thermomechanical
 * mirrors and the SMA family) need downstream.
 */
struct ContinuumTangent {
    arma::mat Lt;                       ///< Tangent operator (6x6).
    std::vector<arma::vec> P_epsilon;   ///< Per-mechanism strain sensitivity P_eps^l (each 6).
    arma::mat invBhat;                  ///< Active-set-masked inverse of Bhat (Nmech x Nmech).
                                        ///< Exposed for thermomechanical UMATs that need it to
                                        ///< build P_theta downstream.
};

/**
 * @brief Assemble the continuum elasto-(visco)plastic tangent following the
 * leading-mechanism framework of simcoon-documentation §7.4.
 *
 * For an UMAT with Nmech leading mechanisms, given:
 *   - the local Jacobian @c Bhat^{lj} = ∂Φ^l/∂σ:κ^j − K^{lj} already built
 *     by the UMAT during the local Newton loop,
 *   - the flux directions @c kappa_j[j] = L·Λ_ε^j,
 *   - the criterion gradients @c dPhidsigma_l[l] = ∂Φ^l/∂σ,
 *   - the increments @c Ds_j[j] = Δs^j over the time step (drives the
 *     active set),
 *
 * this routine computes:
 *
 *     P_ε^l = Σ_m (B̂⁻¹)_{m,l} · L · ∂Φ^m/∂σ
 *     Lt    = L − Σ_l κ^l ⊗ P_ε^l
 *
 * with the standard active-set masking: a mechanism j is treated as
 * inactive when Δs^j ≤ simcoon::iota; the corresponding row/column of
 * B̂ is replaced by the identity before inversion and the result is
 * masked back to zero so an inactive mechanism contributes nothing.
 *
 * This routine reproduces — identically — the tangent-assembly idiom
 * that is currently duplicated across every iterative UMAT in
 * src/Continuum_mechanics/Umat/Mechanical/{Plasticity,SMA,Viscoelasticity,Combined}/
 * and their Thermomechanical mirrors.
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
 * scalar `Bhat`/`Ds` lifted to 1×1 / size-1 quantities. Removes the
 * `{ vec }` brace-init dance from the 9 single-mechanism UMAT call sites.
 */
ContinuumTangent assemble_continuum_tangent(
    double Bhat_scalar,
    const arma::vec& kappa,
    const arma::vec& dPhidsigma,
    double Ds,
    const arma::mat& L);

/**
 * @brief Assemble the Simo–Hughes algorithmic (consistent) tangent
 * following the same leading-mechanism framework as
 * @ref assemble_continuum_tangent but accounting for the variation of the
 * strain-flow direction Λ_ε^j with stress during the increment.
 *
 * Given the converged-state derivative @c dLambda_dsigma_l[j] = ∂Λ_ε^j/∂σ
 * (a 6×6 matrix, per mechanism), the discrete return map
 *
 *     σ_{n+1} = L · (ε_{n+1} − ε^p_n − Σ_j Δs^j Λ_ε^j(σ_{n+1}, …))
 *
 * linearised at convergence yields the modified operator
 *
 *     M = I_6 + L · Σ_j Δs^j (∂Λ_ε^j/∂σ)
 *     L̃ = M⁻¹ L
 *     κ̃^j = M⁻¹ κ^j
 *     B̃^{lj} = ∂Φ^l/∂σ · κ̃^j − K^{lj} = Bhat^{lj} + ∂Φ^l/∂σ · (M⁻¹ − I) · κ^j
 *
 * after which the assembly is structurally identical to the continuum
 * helper, with L̃ in place of L and κ̃ in place of κ:
 *
 *     P_ε^l = Σ_m (B̃⁻¹)_{m,l} · L̃ · ∂Φ^m/∂σ
 *     Lt    = L̃ − Σ_l κ̃^l ⊗ P_ε^l
 *
 * Active-set masking and Bhat regularisation reuse the same Ds_j > iota
 * convention as the continuum helper.
 *
 * **Parity property.** When every @c dLambda_dsigma_l[j] is zero (the
 * frozen-direction assumption underlying the continuum tangent), `M = I`,
 * `L̃ = L`, `κ̃ = κ`, `B̃ = Bhat`, and this routine returns *bit-identical*
 * results to @ref assemble_continuum_tangent. This is what enables the
 * `tangent_mode` dispatch in UMATs that have not yet been extended with
 * ∂Λ/∂σ data — the algorithmic call is a strict superset of the continuum
 * call.
 *
 * **State coupling.** This first version assumes Λ_ε^j has been linearised
 * w.r.t. σ only. The fully consistent tangent additionally requires
 * `∂Λ_ε^j/∂V` and `Λ_V^{j,q}`; that extension will land as an overload
 * once UMATs expose the corresponding data (doc §7.4).
 *
 * @param Bhat_continuum     **Continuum** local Jacobian
 *                           B̂^{lj} = ∂Φ^l/∂σ · κ^j − K^{lj}, i.e. the same
 *                           quantity passed to @ref assemble_continuum_tangent.
 *                           Do **NOT** pass an already-corrected algorithmic B̃ —
 *                           the routine builds B̃ internally and would
 *                           double-correct, silently producing wrong tangents.
 * @param kappa_j            Per-mechanism flux κ^j = L · Λ_ε^j (each 6).
 * @param dPhidsigma_l       Per-mechanism criterion gradient ∂Φ^l/∂σ (each 6).
 * @param Ds_j               Per-mechanism increment Δs^j over the step.
 * @param L                  Elastic stiffness (6×6).
 * @param dLambda_dsigma_l   Per-mechanism ∂Λ_ε^j/∂σ (each 6×6). Pass a
 *                           zero matrix where the direction is independent
 *                           of σ (the algorithmic correction degenerates
 *                           to the continuum tangent in that block).
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

} // namespace simcoon
