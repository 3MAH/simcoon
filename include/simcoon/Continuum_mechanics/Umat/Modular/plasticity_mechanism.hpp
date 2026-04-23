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
 * @file plasticity_mechanism.hpp
 * @brief Plasticity strain mechanism for modular UMAT.
 *
 * Combines yield criterion, isotropic hardening, and kinematic hardening
 * into a complete plasticity model.
 *
 * @see StrainMechanism, YieldCriterion, IsotropicHardening, KinematicHardening, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/yield_criterion.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>

namespace simcoon {

/**
 * @brief Plasticity strain mechanism
 *
 * This class implements rate-independent plasticity by combining:
 * - A yield criterion (von Mises, Hill, etc.)
 * - Isotropic hardening (power-law, Voce, etc.)
 * - Kinematic hardening (Prager, Chaboche, etc.)
 *
 * Internal variables:
 * - "p": accumulated plastic strain (scalar)
 * - "EP": plastic strain tensor (6 Voigt)
 * - "X_i": backstress tensors from kinematic hardening
 */
class PlasticityMechanism final : public StrainMechanism {
private:
    YieldType yield_type_;             ///< Yield criterion type (for deferred configuration)
    std::unique_ptr<YieldCriterion> yield_;
    std::unique_ptr<IsotropicHardening> iso_hard_;
    std::unique_ptr<KinematicHardening> kin_hard_;
    double sigma_Y_;  ///< Initial yield stress

    // IVC keys, cached at register_variables.
    std::string p_key_;
    std::string EP_key_;

    // Per-iteration caches populated by compute_constraints. dPhi_dsigma() and
    // kappa() return const-refs into these single-element buffers to avoid
    // re-constructing std::vector<arma::vec> on every FB iteration.
    mutable arma::vec flow_dir_;                ///< ∂Φ/∂σ (associated flow direction)
    mutable arma::vec kappa_;                   ///< L_ref · flow_dir_
    mutable std::vector<arma::vec> dPhi_dsigma_cache_{arma::zeros(6)};
    mutable std::vector<arma::vec> kappa_cache_{arma::zeros(6)};
    mutable double H_total_{0.0};               ///< Hardening modulus (iso + kin)

public:
    /**
     * @brief Constructor
     * @param yield_type Type of yield criterion
     * @param iso_type Type of isotropic hardening
     * @param kin_type Type of kinematic hardening
     * @param N_iso Number of isotropic hardening terms (for COMBINED_VOCE)
     * @param N_kin Number of kinematic hardening terms (for CHABOCHE)
     */
    PlasticityMechanism(
        YieldType yield_type = YieldType::VON_MISES,
        IsoHardType iso_type = IsoHardType::POWER_LAW,
        KinHardType kin_type = KinHardType::NONE,
        int N_iso = 1,
        int N_kin = 1
    );

    // Default move, no copy (due to unique_ptr)
    PlasticityMechanism(const PlasticityMechanism&) = delete;
    PlasticityMechanism& operator=(const PlasticityMechanism&) = delete;
    PlasticityMechanism(PlasticityMechanism&&) = default;
    PlasticityMechanism& operator=(PlasticityMechanism&&) = default;
    ~PlasticityMechanism() override = default;

    // ========== Configuration ==========

    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    void set_ivc_prefix(const std::string& prefix) override;

    // ========== Mechanism Properties ==========

    [[nodiscard]] MechanismType type() const override { return MechanismType::PLASTICITY; }
    [[nodiscard]] int num_constraints() const override { return 1; }  // Single yield surface

    /**
     * @brief Get the yield criterion
     * @return Reference to yield criterion
     */
    [[nodiscard]] const YieldCriterion& yield_criterion() const noexcept { return *yield_; }

    /**
     * @brief Get the isotropic hardening
     * @return Reference to isotropic hardening
     */
    [[nodiscard]] const IsotropicHardening& isotropic_hardening() const noexcept { return *iso_hard_; }

    /**
     * @brief Get the kinematic hardening
     * @return Reference to kinematic hardening
     */
    [[nodiscard]] const KinematicHardening& kinematic_hardening() const noexcept { return *kin_hard_; }

    /**
     * @brief Get initial yield stress
     * @return sigma_Y
     */
    [[nodiscard]] double initial_yield_stress() const noexcept { return sigma_Y_; }

    // ========== Yield-side queries ==========
    //
    // These helpers expose the stress-side yield machinery (criteria.cpp
    // Eq_stress dispatch) with the backstress shift applied internally, so
    // callers do not need to reach into kin_hard_->total_backstress(ivc).

    /**
     * @brief Equivalent stress of (sigma - X) under the configured yield criterion.
     *
     * Dispatches through YieldCriterion to criteria.cpp::Eq_stress
     * (Mises / Tresca / Drucker / Hill / DFA / Ani).
     */
    [[nodiscard]] double equivalent_stress(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc) const;

    /// tensor2-typed overload.
    [[nodiscard]] double equivalent_stress(
        const tensor2& sigma,
        const InternalVariableCollection& ivc) const;

    /**
     * @brief Yield function value Phi = sigma_eq(sigma - X) - R(p) - sigma_Y.
     *
     * Phi <= 0 inside the elastic domain; the return-mapping drives Phi to 0
     * when the trial state lies outside.
     */
    [[nodiscard]] double yield_function(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc) const;

    [[nodiscard]] double yield_function(
        const tensor2& sigma,
        const InternalVariableCollection& ivc) const;

    /// Flow direction dPhi/dsigma at (sigma - X), strain-typed in the tensor2 overload.
    [[nodiscard]] arma::vec flow_direction(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc) const;

    [[nodiscard]] tensor2 flow_direction(
        const tensor2& sigma,
        const InternalVariableCollection& ivc) const;

    // ========== Constitutive Computations ==========

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

    [[nodiscard]] const std::vector<arma::vec>& dPhi_dsigma(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc) const override;

    [[nodiscard]] const std::vector<arma::vec>& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref,
        const InternalVariableCollection& ivc) const override;

    arma::vec inelastic_strain(const InternalVariableCollection& ivc) const override;

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
};

} // namespace simcoon
