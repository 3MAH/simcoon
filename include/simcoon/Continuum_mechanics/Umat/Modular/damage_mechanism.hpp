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
 * @file damage_mechanism.hpp
 * @brief Scalar damage mechanism for continuum damage mechanics (CDM)
 *
 * This mechanism implements isotropic scalar damage following the
 * effective stress concept:
 *   sigma_eff = sigma / (1 - D)
 *
 * Damage evolution is driven by an energy-based criterion:
 *   Y = (1/2) * sigma : S : sigma  (damage driving force)
 *
 * Damage evolution follows:
 *   D = f(Y_max)
 *
 * where Y_max is the maximum damage driving force reached in history.
 *
 * Supported damage evolution laws:
 * - Linear: D = (Y - Y_0) / (Y_c - Y_0) for Y > Y_0
 *   Props: Y_0, Y_c
 * - Exponential: D = 1 - exp(-A * (Y - Y_0) / (Y_c - Y_0)) for Y > Y_0
 *   Props: Y_0, Y_c, A
 * - Power law: D = ((Y - Y_0) / (Y_c - Y_0))^n for Y > Y_0
 *   Props: Y_0, Y_c, n
 * - Weibull: D = 1 - exp(-((Y - Y_0) / A)^n) for Y > Y_0
 *   Props: Y_0, Y_c, A (scale), n (shape)
 *
 * @see StrainMechanism, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

/**
 * @brief Type of damage evolution law
 */
enum class DamageType {
    LINEAR = 0,         ///< Linear damage evolution
    EXPONENTIAL = 1,    ///< Exponential damage evolution
    POWER_LAW = 2,      ///< Power-law damage evolution
    WEIBULL = 3         ///< Weibull-based damage
};

/**
 * @brief Scalar damage mechanism
 *
 * Implements isotropic damage following continuum damage mechanics.
 * The damage variable D reduces the effective stiffness:
 *   L_damaged = (1 - D) * L
 */
class DamageMechanism final : public StrainMechanism {
private:
    DamageType damage_type_;    ///< Type of damage evolution law

    // Damage parameters
    double Y_0_;                ///< Damage threshold (no damage below this)
    double Y_c_;                ///< Critical damage driving force
    double D_c_;                ///< Critical damage value (default: 0.99)
    double A_;                  ///< Damage evolution parameter
    double n_;                  ///< Power law exponent

    // Cached values
    mutable double Y_current_;  ///< Current damage driving force
    mutable double D_current_;  ///< Current damage
    mutable double dD_dY_;      ///< Derivative of damage w.r.t. driving force
    mutable arma::mat M_cached_; ///< Cached compliance (set from L on first use)
    mutable bool M_cached_valid_; ///< Whether M_cached_ is valid

public:
    /**
     * @brief Constructor
     * @param type Type of damage evolution law (default: LINEAR)
     */
    explicit DamageMechanism(DamageType type = DamageType::LINEAR);

    // StrainMechanism interface
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;

    [[nodiscard]] int num_constraints() const override { return 1; }
    [[nodiscard]] MechanismType type() const override { return MechanismType::DAMAGE; }

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

    // Damage-specific methods

    /**
     * @brief Compute damage driving force from stress
     * @param sigma Stress (6-component Voigt)
     * @param S Compliance tensor (6x6)
     * @return Damage driving force Y
     */
    double compute_driving_force(const arma::vec& sigma, const arma::mat& S) const;

    /**
     * @brief Compute damage from driving force
     * @param Y Damage driving force
     * @param Y_max Maximum driving force in history
     * @return Damage value D
     */
    double compute_damage(double Y, double Y_max) const;

    /**
     * @brief Get the current damage value
     * @param ivc Internal variable collection
     * @return Current damage D
     */
    double get_damage(const InternalVariableCollection& ivc) const;

    /**
     * @brief Get damaged stiffness tensor
     * @param L Undamaged stiffness
     * @param D Damage value
     * @return Damaged stiffness (1-D)*L
     */
    static arma::mat damaged_stiffness(const arma::mat& L, double D);

    // Accessors
    [[nodiscard]] DamageType damage_type() const noexcept { return damage_type_; }
    [[nodiscard]] double Y_0() const noexcept { return Y_0_; }
    [[nodiscard]] double Y_c() const noexcept { return Y_c_; }
    [[nodiscard]] double D_c() const noexcept { return D_c_; }
};

} // namespace simcoon
