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

    // Cached quantities from last constraint computation
    mutable arma::vec flow_dir_;      ///< Flow direction
    mutable arma::vec kappa_;         ///< L * flow_direction
    mutable double dPhi_dp_;          ///< dPhi/dp
    mutable double H_total_;          ///< Total hardening modulus

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

    // ========== Constitutive Computations ==========

    void compute_constraints(
        const arma::vec& sigma,
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
