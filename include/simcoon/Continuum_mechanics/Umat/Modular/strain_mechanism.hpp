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
 * - Owns its internal variables (a mechanism-local InternalVariableCollection)
 * - Contributes constraint functions (yield/evolution equations)
 * - Contributes to the Jacobian matrix
 * - Updates its internal variables during return mapping
 *
 * Owning the state makes a mechanism self-contained: it can back a composed
 * ModularUMAT or a dedicated single-mechanism UMAT without an external
 * registry, and two mechanisms of the same type can never collide on
 * variable names. The orchestrator serializes state by iterating mechanisms
 * (see compute_offsets / pack / unpack below).
 *
 * @see PlasticityMechanism, ViscoelasticMechanism, DamageMechanism, ModularUMAT
 * @version 2.0
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
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
    std::string name_;                  ///< Mechanism name for debugging
    InternalVariableCollection ivc_;    ///< Mechanism-owned internal variables

public:
    /**
     * @brief Constructor
     * @param name Mechanism identifier
     */
    explicit StrainMechanism(const std::string& name) : name_(name) {}

    virtual ~StrainMechanism() = default;

    // ========== State access & statev serialization ==========
    //
    // Non-virtual, implemented once against the owned collection. The
    // orchestrator serializes the composed state as
    //   statev = [orchestrator scalars | mech 0 | mech 1 | ...]
    // by assigning each mechanism a base offset then delegating pack/unpack.

    /// The mechanism's internal-variable collection (state lives here).
    InternalVariableCollection& variables() noexcept { return ivc_; }
    const InternalVariableCollection& variables() const noexcept { return ivc_; }

    /// Assign absolute statev offsets to the owned variables, starting at
    /// @p base_offset. Must be called after register_variables() and before
    /// pack()/unpack().
    void compute_offsets(unsigned int base_offset) { ivc_.compute_offsets(base_offset); }

    /// Number of statev slots consumed by this mechanism's variables.
    [[nodiscard]] unsigned int statev_size() const noexcept { return ivc_.total_statev_size(); }

    /// Serialize the owned variables into statev at their absolute offsets.
    void pack(arma::vec& statev) const { ivc_.pack_all(statev); }

    /// Restore the owned variables from statev (also resets start values).
    void unpack(const arma::vec& statev) { ivc_.unpack_all(statev); }

    /// Co-rotate every objective owned variable with the rotation increment.
    void rotate(const arma::mat& DR) { ivc_.rotate_all(DR); }

    /// Copy current values to start values for all owned variables.
    void to_start() { ivc_.to_start_all(); }

    // ========== Configuration ==========

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Register this mechanism's internal variables into its owned
     * collection. Call once, after configure().
     */
    virtual void register_variables() = 0;

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
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const = 0;

    /**
     * @brief Compute contribution to Jacobian matrix B
     * @param sigma Current stress tensor (6 Voigt)
     * @param L Elastic stiffness tensor (6x6)
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
        arma::mat& B,
        int row_offset
    ) const = 0;

    // ========== Cross-mechanism Jacobian assembly ==========

    /**
     * @brief Return dPhi^l/dsigma for each constraint l this mechanism owns.
     *
     * Used by the orchestrator to assemble the off-diagonal Jacobian entries
     * \f[ B_{lj} = - \frac{\partial \Phi^l}{\partial \boldsymbol{\sigma}}
     *              \cdot \boldsymbol{\kappa}^j \f]
     * The contraction is evaluated on the engineering Voigt components
     * (arma::dot of .voigt()), the work-conjugate pairing of a strain-typed
     * dPhi with a stress-typed kappa.
     *
     * Typing: dPhi/dsigma is strain-typed (VoigtType::strain — the dEq_stress
     * convention, shear slots carry the doubled gamma components).
     *
     * Returned by const-ref into mechanism-internal storage populated during
     * compute_constraints — callers must not retain the reference past the
     * next call to compute_constraints on this mechanism.
     *
     * Default: empty — meaning Phi does not depend on sigma (as in the
     * Prony_Nfast-style viscoelastic mechanism, whose Phi is written in
     * terms of strain and branch-internal EV_i only).
     */
    virtual const std::vector<tensor2>& dPhi_dsigma(
        const arma::vec& sigma) const {
        (void)sigma;
        static const std::vector<tensor2> empty;
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
     *
     * Typing: stress-typed (VoigtType::stress) for plasticity/viscoelasticity
     * (kappa = L·Lambda). DamageMechanism returns a strain-typed kappa
     * (kappa = ∂M/∂D · σ, an M·σ product) — see its override for why the
     * orchestrator's engineering-Voigt dot is kept for that pairing too.
     */
    virtual const std::vector<tensor2>& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref) const = 0;

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
     * Read the other mechanism's state through other.variables().
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
        int j_lead_idx) const {
        (void)l_constraint_idx;
        (void)other;
        (void)j_lead_idx;
        return 0.0;
    }

    // ========== tangent_mode >= tangent_algorithmic hooks (algorithmic tangent / CPP) ==========

    /**
     * @brief Flow-direction Hessians \f$ \mathrm{d}\boldsymbol{\Lambda}_\varepsilon^j/
     * \mathrm{d}\boldsymbol{\sigma} \f$ per constraint, for tangent_mode >= tangent_algorithmic.
     *
     * Compliance-typed tensor4 (a stress → strain map, so the engineering↔Mandel
     * congruence carries the correct shear factors by construction — the Hessian
     * shear coefficient is 1, not 2/3). Evaluated at the same (shifted) stress as
     * dPhi_dsigma; same const-ref lifetime contract.
     *
     * Default: nullptr — the mechanism opts out and keeps its continuum
     * tangent_contribution() in every tangent_mode (correct for mechanisms whose
     * flow does not depend on stress: Prony viscoelasticity, scalar damage).
     */
    [[nodiscard]] virtual const std::vector<tensor4>* dLambda_dsigma(
        const arma::vec& sigma) const {
        (void)sigma;
        return nullptr;
    }

    /**
     * @brief Backward-Euler state refresh from the IVC start values, for the
     * closest-point (tangent_mode == tangent_closest_point) integrator.
     *
     * Contract (matches ReturnStateHooks::update_state of return_mapping.hpp):
     * rebuild the mechanism's internal variables as
     * \f$ \mathbf{V} = \mathbf{V}_n + \sum_j \Delta s^j\,\boldsymbol{\Lambda}_V^j \f$
     * from the *start* values stored in the owned collection and the TOTAL
     * multiplier increments — NOT an incremental update; each call starts over
     * from \f$ \mathbf{V}_n \f$ so the CPP Newton can re-evaluate the state at
     * every iterate. Closed forms preferred (e.g. Armstrong–Frederick backward
     * Euler).
     *
     * @return false if the mechanism does not support the implicit refresh
     * (default) — the orchestrator then falls back to the CCP integrator.
     */
    virtual bool refresh_state(
        const arma::vec& sigma,
        const arma::vec& Ds_total,
        int offset) {
        (void)sigma;
        (void)Ds_total;
        (void)offset;
        return false;
    }

    /**
     * @brief Multiplicative stiffness-reduction factor applied to the elastic
     * stress prediction.
     *
     * CDM-style mechanisms return (1 - D) so the orchestrator computes the
     * nominal stress sigma = f * L : Eel, consistent with the (1 - D) * L
     * scaling their tangent_contribution applies. Factors from several
     * mechanisms compose multiplicatively. Default: 1 (no reduction).
     */
    [[nodiscard]] virtual double stiffness_reduction() const {
        return 1.0;
    }

    /**
     * @brief Get inelastic strain from this mechanism
     * @return Inelastic strain tensor (6 Voigt)
     */
    virtual arma::vec inelastic_strain() const = 0;

    /**
     * @brief Update internal variables given multiplier increments
     * @param ds Multiplier increment vector
     * @param offset Starting index in ds for this mechanism
     */
    virtual void update(
        const arma::vec& ds,
        int offset
    ) = 0;

    /**
     * @brief Compute contribution to consistent tangent
     * @param sigma Current stress tensor
     * @param L Elastic stiffness tensor
     * @param Ds Total multiplier increments
     * @param offset Starting index in Ds for this mechanism
     * @param Lt Tangent modulus to update
     */
    virtual void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        arma::mat& Lt
    ) const = 0;

    // ========== Work Quantities ==========

    /**
     * @brief Compute work decomposition
     * @param sigma_start Stress at start of increment
     * @param sigma Current stress
     * @param Wm_r Output: recoverable work increment
     * @param Wm_ir Output: irrecoverable (stored) work increment
     * @param Wm_d Output: dissipated work increment
     */
    virtual void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const = 0;
};

} // namespace simcoon
