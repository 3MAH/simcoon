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
#include <vector>
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
    std::string name_;        ///< Mechanism name for debugging
    std::string ivc_prefix_;  ///< Prefix for IVC variable keys (empty = bare names)
    bool active_;             ///< Whether this mechanism is currently active

    /// Compose the IVC key for a base variable name; empty prefix yields the
    /// bare name, preserving legacy behavior when a mechanism is used standalone.
    [[nodiscard]] std::string key(const std::string& base) const {
        return ivc_prefix_.empty() ? base : ivc_prefix_ + base;
    }

public:
    /**
     * @brief Constructor
     * @param name Mechanism identifier
     */
    explicit StrainMechanism(const std::string& name) : name_(name), ivc_prefix_(), active_(true) {}

    virtual ~StrainMechanism() = default;

    /// Set a prefix applied to every IVC key registered or queried by this
    /// mechanism. ModularUMAT uses this to disambiguate same-type mechanisms.
    virtual void set_ivc_prefix(const std::string& prefix) { ivc_prefix_ = prefix; }
    [[nodiscard]] const std::string& ivc_prefix() const noexcept { return ivc_prefix_; }

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
     * @param E_total Total mechanical strain at current iterate (Etot + DEtot, 6 Voigt)
     * @param L Elastic stiffness tensor (6x6)
     * @param DTime Time increment
     * @param ivc Internal variable collection
     * @param Phi Output: constraint function values
     * @param Y_crit Output: critical values for convergence
     *
     * E_total is required by rate-dependent mechanisms (viscoelasticity,
     * viscoplasticity) to compute the branch driving force; rate-independent
     * mechanisms (plasticity, damage) ignore it.
     */
    virtual void compute_constraints(
        const arma::vec& sigma,
        const arma::vec& E_total,
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
     *
     * Fills its own diagonal block (and within-block off-diagonals, if any).
     * Cross-mechanism off-diagonals are pre-filled by the orchestrator using
     * dPhi_dsigma() and kappa() before this method is called.
     */
    virtual void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const InternalVariableCollection& ivc,
        arma::mat& B,
        int row_offset
    ) const = 0;

    // ========== Cross-mechanism Jacobian assembly ==========

    /**
     * @brief Return dPhi^l/dsigma for each constraint l this mechanism owns.
     *
     * Used by the orchestrator to assemble the off-diagonal Jacobian entries
     *   B_{lj} = -dot(dPhi^l/dsigma, kappa^j)
     *
     * Returned by const-ref into mechanism-internal storage populated during
     * compute_constraints — callers must not retain the reference past the
     * next call to compute_constraints on this mechanism.
     *
     * Default: empty — meaning Phi does not depend on sigma (as in the
     * Prony_Nfast-style viscoelastic mechanism, whose Phi is written in
     * terms of strain and branch-internal EV_i only).
     */
    virtual const std::vector<arma::vec>& dPhi_dsigma(
        const arma::vec& sigma,
        const InternalVariableCollection& ivc) const {
        (void)sigma;
        (void)ivc;
        static const std::vector<arma::vec> empty;
        return empty;
    }

    /**
     * @brief Return kappa^j = L_ref * Lambda_eps^j (+ stiffness/CTE
     * sensitivity corrections) for each lead variable j this mechanism owns.
     *
     * kappa^j is the "total stress-influence tensor" of lead variable s^j,
     * i.e., delta_sigma += -kappa^j * delta_s^j during the return mapping.
     * Length must equal num_constraints(). Same const-ref lifetime contract
     * as dPhi_dsigma — reference is valid until the next compute_constraints.
     */
    virtual const std::vector<arma::vec>& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref,
        const InternalVariableCollection& ivc) const = 0;

    // ========== Weak-coupling hooks (see theory chapter 7) ==========

    /**
     * @brief K^{lj}: how this mechanism's criterion l depends directly on
     * another mechanism j's internal variables via the chain
     *   K^{lj} = sum_i (∂Φ^l/∂V_i^j) · Λ_i^j
     *
     * Non-zero only for genuine cross-mechanism couplings (e.g. damage-
     * weakened plastic hardening where R(p) depends on D). Default: 0.
     *
     * Must be called only after compute_constraints has run this FB
     * iteration on *both* this mechanism and `other`, since their internal
     * caches (flow direction, kappa, dD_dY, ...) drive the partial
     * derivatives.
     *
     * @param l_constraint_idx Index within this mechanism's constraints.
     * @param other The mechanism whose lead variable s^j we differentiate with
     *              respect to (may be this — diagonal hardening is handled in
     *              compute_jacobian_contribution, not here).
     * @param j_lead_idx Index within other's lead variables.
     */
    virtual double K_cross(
        int l_constraint_idx,
        const StrainMechanism& other,
        int j_lead_idx,
        const InternalVariableCollection& ivc) const {
        (void)l_constraint_idx;
        (void)other;
        (void)j_lead_idx;
        (void)ivc;
        return 0.0;
    }

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
