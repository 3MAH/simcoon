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
 * Damage evolution is driven by the thermodynamic force conjugate to D, the undamaged
 * stored energy:
 *   Y = -dpsi/dD = psi_0 = (1/2) * sigma_eff : S : sigma_eff
 * with S the compliance of the undamaged elasticity block (exact for a linear block; for a
 * hyperelastic one the tangent compliance makes it an approximation of psi_0).
 *
 * Coupling (strain equivalence, Lemaitre): the other mechanisms work on the effective stress,
 * the damage scales the result, sigma = (1 - D) sigma_eff. D is not an FB multiplier: it is the
 * fixed point D = f(max(Y_max at the step start, Y(sigma_eff))), iterated until
 * consistency_residual() vanishes; its tangent is the left factor stress_map().
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
#include <vector>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Modular/strain_mechanism.hpp>

namespace simcoon {

/**
 * @brief Type of damage evolution law.
 *
 * In all laws \f$ Y_{eff} = \max(Y, Y_{max}) \f$ (history), \f$ D = 0 \f$
 * below the threshold \f$ Y_0 \f$, and \f$ D \f$ is capped at the critical
 * value \f$ D_c \f$.
 */
enum class DamageType {
    LINEAR = 0,       ///< \f$ D = \frac{Y_{eff} - Y_0}{Y_c - Y_0} \f$
    EXPONENTIAL = 1,  ///< \f$ D = 1 - \exp\!\left(-A \frac{Y_{eff} - Y_0}{Y_c - Y_0}\right) \f$
    POWER_LAW = 2,    ///< \f$ D = \left(\frac{Y_{eff} - Y_0}{Y_c - Y_0}\right)^{n} \f$
    WEIBULL = 3       ///< \f$ D = 1 - \exp\!\left(-\left(\frac{Y_{eff} - Y_0}{A}\right)^{n}\right) \f$
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
    mutable double Y_current_;          ///< Current damage driving force
    mutable double D_current_;          ///< Current damage
    mutable double dD_dY_;              ///< Derivative of damage w.r.t. driving force
    mutable arma::mat M_cached_;        ///< Cached compliance (raw arma form)
    mutable tensor4   M_cached_t_;      ///< Same compliance, typed Tensor4
                                        ///< (compliance type) — built once per
                                        ///< invalidation; reused by
                                        ///< compute_driving_force, dPhi_dsigma,
                                        ///< kappa, tangent_contribution.
    mutable bool M_cached_valid_;       ///< Whether the cached compliances are valid
    mutable arma::mat L_cached_;        ///< Stiffness the cached compliance was inverted from
    arma::mat kappa_zero_ = arma::mat(6, 1, arma::fill::zeros);   ///< zero: no multiplier carried
    arma::vec::fixed<6> no_inelastic_strain_{arma::fill::zeros};   ///< inelastic_strain(): damage adds none
    size_t D_index_ = 0;        ///< ivc_ index of "D", resolved by register_variables
    size_t Y_max_index_ = 0;    ///< ivc_ index of "Y_max", resolved by register_variables

    // last driving force and the effective stress it was evaluated at, for the cached compliance (reset with it)
    mutable arma::vec::fixed<6> Y_sigma_;
    mutable double Y_sigma_value_ = 0.0;
    mutable bool Y_sigma_valid_ = false;

    /// compute_driving_force(sigma_eff, M_cached_), memoised on the last stress
    double driving_force(const arma::vec& sigma_eff) const;

public:
    /**
     * @brief Constructor
     * @param type Type of damage evolution law (default: LINEAR)
     */
    explicit DamageMechanism(DamageType type = DamageType::LINEAR);

    /**
     * @brief Reads the damage parameters from @p props at @p offset (advanced past them).
     *
     * Layout: props[offset] = \f$ Y_0 \f$ (damage threshold), props[offset+1] = \f$ Y_c \f$
     * (critical driving force, \f$ Y_c > Y_0 \f$), then the parameters of the law: none (LINEAR),
     * \f$ A \f$ (EXPONENTIAL), \f$ n \f$ (POWER_LAW), \f$ A, n \f$ (WEIBULL). The law itself comes
     * from the mechanism header, not from these props; \f$ D_c \f$ is fixed at 0.99.
     * @throw std::invalid_argument when \f$ Y_c \le Y_0 \f$
     */
    void configure(const arma::vec& props, int& offset) override;
    void register_variables() override;

    [[nodiscard]] int num_constraints() const override { return 1; }
    [[nodiscard]] MechanismType type() const override { return MechanismType::DAMAGE; }

    void compute_constraints(
        const arma::vec& sigma,
        const arma::vec& E_total,
        const arma::mat& L,
        double DTime,
        arma::vec& Phi,
        arma::vec& Y_crit
    ) const override;
    void compute_jacobian_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        arma::mat& B,
        int row_offset
    ) const override;

    /// Nominal-stress CDM: the orchestrator scales its elastic prediction by
    /// this (1 - D) factor, matching the (1 - D) * L scaling applied in
    /// tangent_contribution.
    [[nodiscard]] double stiffness_reduction() const override;

    /**
     * @brief Linearised map from the effective to the damaged stress, \f$ d\boldsymbol{\sigma} = \mathbf{Q}\,d\boldsymbol{\sigma}_{eff} \f$.
     *
     * With \f$ \boldsymbol{\sigma} = (1 - D)\,\boldsymbol{\sigma}_{eff} \f$ and, while damage grows,
     * \f$ dD = D'(Y)\,(\mathbf{S}\boldsymbol{\sigma}_{eff}) \cdot d\boldsymbol{\sigma}_{eff} \f$:
     * \f$ \mathbf{Q} = (1 - D)\,\mathbf{I} - D'(Y)\,\boldsymbol{\sigma}_{eff} \otimes \mathbf{S}\boldsymbol{\sigma}_{eff} \f$
     * (\f$ (1 - D)\,\mathbf{I} \f$ under unloading). Used by the orchestrator to couple damage
     * exactly with the other mechanisms in the tangent. Call after compute_constraints.
     * @param sigma_eff effective stress (6)
     * @return \f$ \mathbf{Q} \f$ (6x6)
     */
    [[nodiscard]] arma::mat stress_map(const arma::vec& sigma_eff) const override;

    /// @brief |max(Y_max at step start, Y(sigma_eff)) - Y_max used for D|, relative: D consistent with the stress
    [[nodiscard]] double consistency_residual(const arma::vec& sigma_eff) const override;

    /// Explicitly integrated, energy-release-typed row — excluded from the
    /// drift-guard ARMING count (see StrainMechanism::guarded_constraints);
    /// the commit checks themselves are no-ops here through the default
    /// multiplier_cap()/state_drift().
    [[nodiscard]] bool guarded_constraints() const override { return false; }
    [[nodiscard]] bool carries_multipliers() const override { return false; }   // Phi = -1, D by fixed point

    /// kappa is an M·σ product here, i.e. STRAIN-convention — unlike the
    /// stress-convention kappa of plasticity/viscoelasticity. The orchestrator's
    /// engineering-Voigt dot therefore mixes conventions for damage columns; this
    /// reproduces the established Prony/CCP-era numerics and is kept deliberately.
    /// Review as a modeling item, not a port bug. dPhi_dsigma is the base default
    /// (no column): the row is a fixed point, not an FB multiplier.
    [[nodiscard]] const arma::mat& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref) const override;

    const arma::vec& inelastic_strain() const override;

    void update(
        const arma::vec& ds,
        int offset
    ) override;

    void tangent_contribution(
        const arma::vec& sigma,
        const arma::mat& L,
        const arma::vec& Ds,
        int offset,
        arma::mat& Lt
    ) const override;

    /**
     * @brief Dissipation of the increment, the trapezoid \f$ \tfrac12 (Y_n + Y_{n+1})\,\Delta D \f$ on
     * the driving force (the quadrature of every work of the engine; exact for the LINEAR law,
     * where it telescopes to \f$ (Y_{n+1}^2 - Y_n^2)/(2 (Y_c - Y_0)) \f$), \f$ Y_n \f$ from the start
     * stress, \f$ \boldsymbol{\sigma}^n_{eff} = \boldsymbol{\sigma}^n/(1 - D_n) \f$. Nothing stored.
     */
    void compute_work(
        const arma::vec& sigma_start,
        const arma::vec& sigma,
        const arma::vec& E_start,
        const arma::vec& E_end,
        double& Wm_r,
        double& Wm_ir,
        double& Wm_d
    ) const override;
    /**
     * @brief Derivatives of the trapezoid dissipation: \f$ \tfrac12 \Delta D\,\partial Y_{n+1} +
     * \tfrac12 (Y_n + Y_{n+1})\,D'(Y)\,\partial Y_{n+1} \f$ (\f$ D'(Y) = 0 \f$ under unloading or at
     * \f$ D_c \f$), \f$ \partial Y_{n+1}/\partial\boldsymbol{\varepsilon} = (\partial\boldsymbol{\sigma}_{eff}/\partial\boldsymbol{\varepsilon})^T
     * \mathbf{S}\boldsymbol{\sigma}_{eff} \f$ and the same with \f$ \partial\boldsymbol{\sigma}_{eff}/\partial\theta \f$
     * (the @c dSdE_eff / @c dSdT_eff of the inputs). The row carries no multiplier: the multiplier
     * sensitivities are not read.
     * Call after compute_constraints (\f$ D'(Y) \f$ and the compliance of the committed state).
     */
    void dissipation_sensitivity(const DissipationInputs& in, arma::vec& dDgamma_dE, double& dDgamma_dT) const override;

    // Damage-specific methods

    /**
     * @brief Damage driving force \f$ Y = \psi_0 = \tfrac12\,\boldsymbol{\sigma}_{eff} : \mathbf{S} : \boldsymbol{\sigma}_{eff} \f$
     * @param sigma_eff Effective stress (6-component Voigt)
     * @param S Compliance tensor (6x6, the cached one is used)
     * @return Damage driving force Y
     */
    double compute_driving_force(const arma::vec& sigma_eff, const arma::mat& S) const;

    /**
     * @brief Compute damage from the history driving force
     * @param Y_eff driving force the damage law is evaluated at (history maximum)
     * @return Damage value D
     */
    double compute_damage(double Y_eff) const;

    /**
     * @brief Get the current damage value
     * @return Current damage D
     */
    double get_damage() const;

    /**
     * @brief Invalidate the cached compliance M_cached_.
     *
     * M_cached_ = inv(L) is keyed on the stiffness L passed into
     * compute_constraints — the tangent of the current elastic state — and is
     * rebuilt whenever L changes, including at every refresh of a
     * state-dependent (hyperelastic) block. This call is therefore never
     * required; it only forces the next compute_constraints to rebuild.
     */
    void invalidate_compliance_cache() noexcept { M_cached_valid_ = false; }

    // Accessors
    [[nodiscard]] DamageType damage_type() const noexcept { return damage_type_; }
    [[nodiscard]] double Y_0() const noexcept { return Y_0_; }
    [[nodiscard]] double Y_c() const noexcept { return Y_c_; }
    [[nodiscard]] double D_c() const noexcept { return D_c_; }
};

} // namespace simcoon
