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
 * @file strain_mechanism.hpp
 * @brief Base class for strain mechanisms in modular UMAT.
 *
 * Strain mechanisms are pluggable components that contribute to the
 * total strain decomposition: Etot = Eel + E_plastic + E_viscous + ...
 *
 * Each mechanism:
 * - Registers its internal variables
 * - Contributes constraint functions (yield/evolution equations)
 * - Contributes to the Jacobian matrix
 * - Updates its internal variables during return mapping
 *
 * @see PlasticityMechanism, ViscoelasticMechanism, DamageMechanism, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <string>
#include <map>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>

namespace simcoon {

/**
 * @brief Types of strain mechanisms
 */
enum class MechanismType {
    PLASTICITY = 0,       ///< Rate-independent plasticity
    VISCOELASTICITY = 1,  ///< Viscoelasticity (Prony series)
    DAMAGE = 2            ///< Damage mechanics
};

/**
 * @brief Base class for strain mechanisms
 *
 * This abstract class defines the interface for strain mechanisms that
 * can be composed in a modular UMAT. Each mechanism contributes to the
 * total strain decomposition and provides constraint functions for the
 * return mapping algorithm.
 */
class StrainMechanism {
protected:
    std::string name_;  ///< Mechanism name for debugging
    bool active_;       ///< Whether this mechanism is currently active

public:
    /**
     * @brief Constructor
     * @param name Mechanism identifier
     */
    explicit StrainMechanism(const std::string& name) : name_(name), active_(true) {}

    virtual ~StrainMechanism() = default;

    // ========== Configuration ==========

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Register internal variables for this mechanism
     * @param ivc Internal variable collection
     */
    virtual void register_variables(InternalVariableCollection& ivc) = 0;

    // ========== Mechanism Properties ==========

    /**
     * @brief Get mechanism name
     * @return Mechanism identifier
     */
    [[nodiscard]] const std::string& name() const noexcept { return name_; }

    /**
     * @brief Get mechanism type
     * @return Mechanism type enum
     */
    [[nodiscard]] virtual MechanismType type() const = 0;

    /**
     * @brief Check if mechanism is active
     * @return True if mechanism contributes to constitutive response
     */
    [[nodiscard]] bool is_active() const noexcept { return active_; }

    /**
     * @brief Set mechanism active state
     * @param active Whether mechanism should be active
     */
    void set_active(bool active) noexcept { active_ = active; }

    /**
     * @brief Get number of constraint equations
     * @return Number of active constraints
     */
    [[nodiscard]] virtual int num_constraints() const = 0;

    // ========== Constitutive Computations ==========

    /**
     * @brief Compute constraint functions (yield/evolution equations)
     * @param sigma Current stress tensor (6 Voigt)
     * @param L Elastic stiffness tensor (6x6)
     * @param DTime Time increment
     * @param ivc Internal variable collection
     * @param Phi Output: constraint function values
     * @param Y_crit Output: critical values for convergence
     */
    virtual void compute_constraints(
        const arma::vec& sigma,
        const arma::mat& L,
        double DTime,
        const InternalVariableCollection& ivc,
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const = 0;

    /**
     * @brief Compute flow directions for internal variables
     * @param sigma Current stress tensor (6 Voigt)
     * @param ivc Internal variable collection
     * @param Lambda_map Output: map of variable name to flow direction
     */
    virtual void compute_flow_directions(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        std::map<std::string, arma::vec>& Lambda_map
    ) const = 0;

    /**
     * @brief Compute contribution to Jacobian matrix B
     * @param sigma Current stress tensor (6 Voigt)
     * @param L Elastic stiffness tensor (6x6)
     * @param ivc Internal variable collection
     * @param B Jacobian matrix to update
     * @param row_offset Starting row for this mechanism's contributions
     */
    virtual void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const InternalVariableCollection& ivc,
        arma::mat& B,
        int row_offset
    ) const = 0;

    /**
     * @brief Get inelastic strain from this mechanism
     * @param ivc Internal variable collection
     * @return Inelastic strain tensor (6 Voigt)
     */
    virtual arma::vec inelastic_strain(const InternalVariableCollection& ivc) const = 0;

    /**
     * @brief Update internal variables given multiplier increments
     * @param ds Multiplier increment vector
     * @param offset Starting index in ds for this mechanism
     * @param ivc Internal variable collection to update
     */
    virtual void update(
        const arma::vec& ds,
        int offset,
        InternalVariableCollection& ivc
    ) = 0;

    /**
     * @brief Compute contribution to consistent tangent
     * @param sigma Current stress tensor
     * @param L Elastic stiffness tensor
     * @param Ds Total multiplier increments
     * @param offset Starting index in Ds for this mechanism
     * @param ivc Internal variable collection
     * @param Lt Tangent modulus to update
     */
    virtual void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        const InternalVariableCollection& ivc,
        arma::mat& Lt
    ) const = 0;

    // ========== Work Quantities ==========

    /**
     * @brief Compute work decomposition
     * @param sigma_start Stress at start of increment
     * @param sigma Current stress
     * @param ivc Internal variable collection
     * @param Wm_r Output: recoverable work increment
     * @param Wm_ir Output: irrecoverable (stored) work increment
     * @param Wm_d Output: dissipated work increment
     */
    virtual void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        const InternalVariableCollection& ivc,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const = 0;
};

} // namespace simcoon
