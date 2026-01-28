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
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>
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
 */
class ModularUMAT {
private:
    // Modules
    ElasticityModule elasticity_;
    std::vector<std::unique_ptr<StrainMechanism>> mechanisms_;
    InternalVariableCollection ivc_;

    // State
    double T_init_;         ///< Initial temperature
    arma::vec sigma_start_; ///< Stress at start of increment
    bool initialized_;      ///< Whether initialize() has been called

    // Solver parameters
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
     * @param props Material properties (g_i, tau_i pairs for each term)
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
     * - props[1..N_el]: elasticity parameters
     * - props[N_el+1]: num_mechanisms
     * - For each mechanism:
     *   - props[i]: mechanism_type (0=plasticity, 1=viscoelasticity, 2=damage)
     *   - For plasticity (type 0):
     *     - yield_type, iso_type, kin_type, N_iso, N_kin, sigma_Y, [yield params], [iso params], [kin params]
     *   - For viscoelasticity (type 1):
     *     - N_prony, then N_prony pairs of (g_i, tau_i)
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
     * @brief Get mechanism by index
     * @param i Index
     * @return Reference to mechanism
     */
    StrainMechanism& mechanism(size_t i) { return *mechanisms_[i]; }
    const StrainMechanism& mechanism(size_t i) const { return *mechanisms_[i]; }

    /**
     * @brief Get internal variable collection
     * @return Reference to internal variable collection
     */
    InternalVariableCollection& internal_variables() { return ivc_; }
    const InternalVariableCollection& internal_variables() const { return ivc_; }

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

    // ========== Solver Parameters ==========

    /**
     * @brief Set maximum iterations
     * @param maxiter Maximum number of return mapping iterations
     */
    void set_max_iterations(int maxiter) { maxiter_ = maxiter; }

    /**
     * @brief Set convergence precision
     * @param precision Relative tolerance for convergence
     */
    void set_precision(double precision) { precision_ = precision; }

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
        double& tnew_dt
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
     * @brief Compute consistent tangent modulus
     * @param sigma Current stress
     * @param Ds_total Total multiplier increments
     * @param Lt Output: tangent modulus
     */
    void compute_tangent(
        const arma::vec& sigma,
        const arma::vec& Ds_total,
        arma::mat& Lt
    );
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
    double& tnew_dt
);

} // namespace simcoon
