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
 * @file modular_umat.hpp
 * @brief Modular UMAT orchestrator for composing constitutive models.
 *
 * This class orchestrates the composition of:
 * - Elasticity module (isotropic, cubic, transverse isotropic, orthotropic)
 * - Multiple strain mechanisms (plasticity, viscoelasticity, damage)
 *
 * It provides a unified return mapping algorithm that handles all
 * coupled constraints from the different mechanisms.
 *
 * Usage example (isotropic elasticity + von Mises plasticity with Voce hardening):
 * @code
 *   ModularUMAT mumat;
 *   int offset = 0;
 *   arma::vec props = {210000, 0.3, 1.2e-5, 300, 100, 10};
 *   mumat.set_elasticity(ElasticityType::ISOTROPIC, props, offset);
 *   mumat.add_plasticity(YieldType::VON_MISES, IsoHardType::VOCE,
 *                        KinHardType::NONE, 1, 1, props, offset);
 *   mumat.initialize(nstatev, statev);
 *   mumat.run(...);
 * @endcode
 *
 * @see ElasticityModule, StrainMechanism, PlasticityMechanism,
 *      ViscoelasticMechanism, DamageMechanism
 * @version 1.0
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/elasticity_module.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

// Forward declarations
class PlasticityMechanism;
class ViscoelasticMechanism;
class DamageMechanism;
enum class YieldType;
enum class IsoHardType;
enum class KinHardType;
enum class DamageType;

/**
 * @brief Modular UMAT orchestrator
 *
 * This class provides a composable constitutive model framework where
 * different elasticity types and strain mechanisms can be combined.
 *
 * State ownership: each StrainMechanism owns its internal variables. The
 * orchestrator serializes the composed state as
 *   statev = [T_init | mechanism 0 | mechanism 1 | ...]
 * in composition order, by assigning each mechanism a base offset at
 * initialize() and delegating pack/unpack/rotate/to_start per mechanism.
 */
class ModularUMAT {
private:
    // Modules
    ElasticityModule elasticity_;
    std::vector<std::unique_ptr<StrainMechanism>> mechanisms_;

    // Constraint-row offset per mechanism; computed once in initialize() and
    // reused by every FB iteration (it only changes if mechanisms are
    // added/removed, which is disallowed after initialize()).
    std::vector<int> mech_offset_;

    // State
    double T_init_;         ///< Initial temperature
    arma::vec sigma_start_; ///< Stress at start of increment
    bool initialized_;      ///< Whether initialize() has been called

    // Return-mapping controls (local internal-equilibrium Newton loop)
    int maxiter_;           ///< Maximum iterations for return mapping
    double precision_;      ///< Convergence tolerance

public:
    /**
     * @brief Default constructor
     */
    ModularUMAT();

    // Disable copy, allow move
    ModularUMAT(const ModularUMAT&) = delete;
    ModularUMAT& operator=(const ModularUMAT&) = delete;
    ModularUMAT(ModularUMAT&&) = default;
    ModularUMAT& operator=(ModularUMAT&&) = default;
    ~ModularUMAT() = default;

    // ========== Configuration ==========

    /**
     * @brief Set elasticity module
     * @param type Elasticity type
     * @param props Material properties
     * @param offset Current offset in props (will be updated)
     */
    void set_elasticity(ElasticityType type, const arma::vec& props, int& offset);

    /**
     * @brief Add a plasticity mechanism
     * @param yield_type Type of yield criterion
     * @param iso_type Type of isotropic hardening
     * @param kin_type Type of kinematic hardening
     * @param N_iso Number of isotropic hardening terms
     * @param N_kin Number of kinematic hardening terms
     * @param props Material properties
     * @param offset Current offset in props (will be updated)
     * @return Reference to the added mechanism
     */
    PlasticityMechanism& add_plasticity(
        YieldType yield_type,
        IsoHardType iso_type,
        KinHardType kin_type,
        int N_iso,
        int N_kin,
        const arma::vec& props,
        int& offset
    );

    /**
     * @brief Add a viscoelastic mechanism (Prony series)
     * @param N_prony Number of Prony terms
     * @param props Material properties: per Prony branch a quadruple
     *        (E_i, nu_i, etaB_i, etaS_i) — same layout as the legacy PRONK/ZENNK
     * @param offset Current offset in props (will be updated)
     * @return Reference to the added mechanism
     */
    ViscoelasticMechanism& add_viscoelasticity(
        int N_prony,
        const arma::vec& props,
        int& offset
    );

    /**
     * @brief Add a damage mechanism
     * @param damage_type Type of damage evolution law
     * @param props Material properties
     * @param offset Current offset in props (will be updated)
     * @return Reference to the added mechanism
     */
    DamageMechanism& add_damage(
        DamageType damage_type,
        const arma::vec& props,
        int& offset
    );

    /**
     * @brief Configure from props array
     *
     * Props format:
     * - props[0]: elasticity_type (0=iso, 1=cubic, 2=trans_iso, 3=ortho)
     * - props[1]: elastic-constant convention code (IsoConv/CubicConv/...,
     *   selects the interpretation of the constant slots — see
     *   ElasticityModule::configure)
     * - props[2..N_el]: elasticity parameters
     * - props[N_el+1]: num_mechanisms
     * - For each mechanism:
     *   - props[i]: mechanism_type (0=plasticity, 1=viscoelasticity, 2=damage)
     *   - For plasticity (type 0):
     *     - yield_type, iso_type, kin_type, N_iso, N_kin, sigma_Y, [yield params], [iso params], [kin params]
     *   - For viscoelasticity (type 1):
     *     - N_prony, then N_prony quadruples of (E_i, nu_i, etaB_i, etaS_i)
     *   - For damage (type 2):
     *     - damage_type (0=linear, 1=exp, 2=power, 3=weibull), Y_0, Y_c, [type-specific params]
     *
     * @param props Material properties vector
     * @param offset Starting offset (default: 0)
     */
    void configure_from_props(const arma::vec& props, int offset = 0);

    /**
     * @brief Initialize internal variables
     * @param nstatev Number of state variables
     * @param statev State variable vector
     */
    void initialize(int nstatev, arma::vec& statev);

    // ========== Accessors ==========

    /**
     * @brief Get the elasticity module
     * @return Const reference to elasticity module
     */
    [[nodiscard]] const ElasticityModule& elasticity() const noexcept { return elasticity_; }

    /**
     * @brief Get number of mechanisms
     * @return Number of strain mechanisms
     */
    [[nodiscard]] size_t num_mechanisms() const noexcept { return mechanisms_.size(); }

    /**
     * @brief Get mechanism by index (its internal variables are reachable
     * through StrainMechanism::variables())
     * @param i Index
     * @return Reference to mechanism
     */
    StrainMechanism& mechanism(size_t i) { return *mechanisms_[i]; }
    const StrainMechanism& mechanism(size_t i) const { return *mechanisms_[i]; }

    /// Controls of the local return-mapping Newton loop — the internal
    /// (material-point) equilibrium iteration, not the global solver.
    /// Defaults: 100 iterations, 1e-9.
    void set_return_mapping_params(int maxiter, double precision) {
        maxiter_ = maxiter;
        precision_ = precision;
    }

    /**
     * @brief Get total number of state variables required
     * @return nstatev
     */
    int required_nstatev() const;

    /**
     * @brief Check if initialized
     * @return True if initialize() has been called
     */
    [[nodiscard]] bool is_initialized() const noexcept { return initialized_; }

    // ========== Main UMAT Entry Point ==========

    /**
     * @brief Run the constitutive model update
     *
     * This is the main entry point that performs:
     * 1. Unpack state variables
     * 2. Apply rotation for objectivity
     * 3. Elastic prediction
     * 4. Return mapping (if inelastic)
     * 5. Compute consistent tangent
     * 6. Compute work quantities
     * 7. Pack state variables
     *
     * @param umat_name UMAT identifier
     * @param Etot Total strain at end of increment
     * @param DEtot Strain increment
     * @param sigma Output: stress at end of increment
     * @param Lt Output: consistent tangent modulus
     * @param L Output: elastic stiffness
     * @param DR Rotation increment matrix
     * @param nprops Number of properties
     * @param props Material properties
     * @param nstatev Number of state variables
     * @param statev State variables
     * @param T Temperature at end of increment
     * @param DT Temperature increment
     * @param Time Current time
     * @param DTime Time increment
     * @param Wm Total mechanical work
     * @param Wm_r Recoverable work
     * @param Wm_ir Irrecoverable work
     * @param Wm_d Dissipated work
     * @param ndi Number of direct stress components
     * @param nshr Number of shear stress components
     * @param start True if first increment
     * @param tnew_dt Suggested new time step ratio
     */
    void run(
        const std::string& umat_name,
        const arma::vec& Etot,
        const arma::vec& DEtot,
        arma::vec& sigma,
        arma::mat& Lt,
        arma::mat& L,
        const arma::mat& DR,
        int nprops,
        const arma::vec& props,
        int nstatev,
        arma::vec& statev,
        double T,
        double DT,
        double Time,
        double DTime,
        double& Wm,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d,
        int ndi,
        int nshr,
        bool start,
        double& tnew_dt,
        int tangent_mode = tangent_default
    );

private:
    /**
     * @brief Perform return mapping algorithm
     *
     * Uses Newton iteration with Fischer-Burmeister complementarity
     * to solve the coupled constraint equations from all mechanisms.
     * The algorithm:
     * 1. Computes elastic prediction (trial stress)
     * 2. Evaluates constraint functions (Phi) from all mechanisms
     * 3. Solves for multiplier increments via Fischer-Burmeister
     * 4. Updates internal variables and recomputes stress
     * 5. Repeats until convergence (error < precision_)
     *
     * @param Etot Total strain at start of increment
     * @param DEtot Strain increment
     * @param sigma Output: stress at end of increment
     * @param T_init Reference temperature
     * @param T Current temperature
     * @param DT Temperature increment
     * @param DTime Time increment
     * @param ndi Number of direct stress components
     * @param Ds_total Output: converged multiplier increments
     *
     * @see Fischer_Burmeister_m() in num_solve.hpp
     */
    void return_mapping(
        const arma::vec& Etot,
        const arma::vec& DEtot,
        arma::vec& sigma,
        double T_init,
        double T,
        double DT,
        double DTime,
        int ndi,
        arma::vec& Ds_total
    );

    /**
     * @brief Compute the tangent modulus.
     *
     * tangent_none (0): no assembly — Lt stays the elastic operator
     * (explicit integration).
     *
     * tangent_continuum (1): each mechanism applies its continuum
     * tangent_contribution() in composition order (pre-2.0 mode 0).
     *
     * tangent_algorithmic (2): mechanisms exposing a flow Hessian
     * (dLambda_dsigma() != nullptr, stress-dependent Phi) are assembled
     * together through assemble_algorithmic_tangent() — the Simo–Hughes
     * consistent operator of the coupled sub-system, built from the converged
     * local Jacobian \f$ \hat{B} = -B \f$ and the mechanism caches. The
     * remaining mechanisms (Prony viscoelasticity, scalar damage — flows
     * independent of stress) keep their continuum contribution, applied on
     * top in composition order, exactly as in the continuum mode.
     *
     * @param sigma Current (converged) stress
     * @param Ds_total Total multiplier increments
     * @param Lt Output: tangent modulus
     * @param tangent_mode tangent_* constant (parameter.hpp): 0 = none
     *        (Lt = elastic L), 1 = continuum, 2 = algorithmic (Simo–Hughes),
     *        3 = closest-point (reserved, throws)
     */
    void compute_tangent(
        const arma::vec& sigma,
        const arma::vec& Ds_total,
        arma::mat& Lt,
        int tangent_mode = tangent_default
    );

    /**
     * @brief Assemble the local multiplier Jacobian B (phases 2 and 3 of the
     * return mapping) from the mechanism caches at the current state.
     * Shared by the FB loop and the mode-1 tangent assembly.
     */
    void assemble_jacobian(
        const arma::vec& sigma,
        double DT,
        arma::mat& B
    );

    /**
     * @brief Product of the mechanisms' stiffness-reduction factors (CDM
     * (1-D) scaling) applied to the elastic stress prediction. 1 when no
     * damage-type mechanism is present.
     */
    [[nodiscard]] double stiffness_reduction() const;
};

/**
 * @brief UMAT function for modular constitutive model
 *
 * Standard UMAT interface for the modular constitutive model.
 * Can be registered in umat_smart.cpp for dispatch.
 */
void umat_modular(
    const std::string& umat_name,
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    arma::mat& Lt,
    arma::mat& L,
    const arma::mat& DR,
    const int& nprops,
    const arma::vec& props,
    const int& nstatev,
    arma::vec& statev,
    const double& T,
    const double& DT,
    const double& Time,
    const double& DTime,
    double& Wm,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d,
    const int& ndi,
    const int& nshr,
    const bool& start,
    double& tnew_dt,
    const int& tangent_mode = tangent_default
);

} // namespace simcoon
