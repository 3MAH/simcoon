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
 * @file viscoelastic_mechanism.hpp
 * @brief Generalized Maxwell viscoelasticity, porting the Prony_Nfast kernel.
 *
 * Each Prony branch i is a Maxwell element with its own stiffness L_i and
 * viscosity tensor H_i (bulk/shear split): driving stress L_i (eps - EV_i)
 * produces a strain rate invH_i . L_i (eps - EV_i). The lead scalar v_i is
 * the accumulated flow magnitude ("path length") and EV_i is the branch
 * viscous strain tensor. The constraint is an equality:
 *
 *   Phi_i = || invH_i . L_i . (eps + Δeps - EV_i) ||_strain - Δv_i / Δt = 0
 *
 * Solved via the same FB complementarity residual used by plasticity; for
 * pure viscoelasticity (no activation threshold) the FB "inactive" branch is
 * unreachable except at the trivial rest state, so FB degrades cleanly to
 * the equality case.
 *
 * Props per branch (4 values): E_i, nu_i, etaB_i, etaS_i.
 * Statev per branch (7 values): v_i (scalar), EV_i (6 Voigt).
 *
 * @see StrainMechanism, ModularUMAT, Prony_Nfast.cpp (original kernel)
 * @version 2.0
 */

#pragma once

#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

class ViscoelasticMechanism final : public StrainMechanism {
private:
    int N_prony_;
    std::vector<double> E_i_;       ///< Branch Young's moduli
    std::vector<double> nu_i_;      ///< Branch Poisson ratios
    std::vector<double> etaB_i_;    ///< Branch bulk viscosities
    std::vector<double> etaS_i_;    ///< Branch shear viscosities

    // Cached per-branch tensors (set in configure; constant over the step)
    std::vector<arma::mat> L_i_;        ///< Branch stiffnesses
    std::vector<arma::mat> H_i_;        ///< Branch viscosity tensors
    std::vector<arma::mat> invH_i_;     ///< Cached inv(H_i)
    std::vector<arma::mat> M0_L_i_;     ///< Cached M_0 · L_i (for inelastic_strain)

    arma::mat M_0_;                     ///< Reference compliance = inv(L_0), set by orchestrator

    // Fully-prefixed IVC keys, cached at register_variables (avoids string
    // concatenation and key() calls inside the FB hot loop).
    std::vector<std::string> ev_key_;   ///< key("EV_" + i) per branch
    std::vector<std::string> v_key_;    ///< key("v_" + i) per branch

    // Per-iteration caches (CCP: frozen flow direction from previous iteration)
    mutable std::vector<arma::vec> flow_i_;       ///< Strain-rate vector per branch
    mutable std::vector<arma::vec> Lambda_i_;     ///< eta_norm_strain(flow_i)
    mutable std::vector<arma::vec> kappa_i_;      ///< L_i . Lambda_i (used in tangent)
    mutable std::vector<arma::vec> dPhi_i_dv_;    ///< invH_i . (eta_norm_strain(flow_i) % Ir05())
    mutable arma::vec K_diag_;                    ///< Cached diagonal K(i,i) = -dPhi_i_dv . kappa_i - 1/Δt

public:
    explicit ViscoelasticMechanism(int N_prony);

    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;

    [[nodiscard]] int num_constraints() const override { return N_prony_; }
    [[nodiscard]] MechanismType type() const override { return MechanismType::VISCOELASTICITY; }

    void compute_constraints(
        const arma::vec& sigma,
        const arma::vec& E_total,
        const arma::mat& L,
        double DTime,
        const InternalVariableCollection& ivc,
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const override;

    void compute_flow_directions(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        std::map<std::string, arma::vec>& Lambda_map
    ) const override;

    void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const InternalVariableCollection& ivc,
        arma::mat& B,
        int row_offset
    ) const override;

    // Viscoelastic Phi is written in strain form (Prony_Nfast) and does not
    // depend on sigma, so dPhi_dsigma returns empty → no cross-mechanism
    // coupling inbound from stress. The mechanism still contributes kappa to
    // outbound stress perturbations seen by other mechanisms.
    [[nodiscard]] const std::vector<arma::vec>& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref,
        const InternalVariableCollection& ivc) const override;

    [[nodiscard]] arma::vec inelastic_strain(const InternalVariableCollection& ivc) const override;

    void update(
        const arma::vec& ds,
        int offset,
        InternalVariableCollection& ivc
    ) override;

    void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        const InternalVariableCollection& ivc,
        arma::mat& Lt
    ) const override;

    void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const override;

    // Accessors
    [[nodiscard]] int num_prony_terms() const noexcept { return N_prony_; }
    [[nodiscard]] double E(int i)    const { return E_i_[i]; }
    [[nodiscard]] double nu(int i)   const { return nu_i_[i]; }
    [[nodiscard]] double etaB(int i) const { return etaB_i_[i]; }
    [[nodiscard]] double etaS(int i) const { return etaS_i_[i]; }
    [[nodiscard]] const arma::mat& L_branch(int i) const { return L_i_[i]; }

    /**
     * @brief Set reference compliance used to convert branch-stored viscous
     * strain into the mechanism's contribution to the total inelastic strain.
     * Called by ModularUMAT after the elasticity module is configured.
     */
    void set_reference_stiffness(const arma::mat& L_0);
};

} // namespace simcoon
