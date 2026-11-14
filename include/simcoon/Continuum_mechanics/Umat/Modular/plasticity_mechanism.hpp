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
 *
 * @note **Integration and tangent.** The return mapping is a cutting plane: each FB
 * iterate adds \f$ \delta p\,\mathbf{n}_k \f$ with the flow direction
 * \f$ \mathbf{n}_k \f$ of that iterate. For von Mises with isotropic hardening
 * \f$ \mathbf{n} \f$ is fixed over the step (radial return), the update is the closest-point
 * one and the algorithmic tangent is its exact derivative. With kinematic hardening
 * (\f$ \mathbf{n} \f$ follows \f$ \boldsymbol{\sigma} - \mathbf{X} \f$) or a Hill, DFA,
 * Drucker or Tresca criterion \f$ \mathbf{n} \f$ rotates between iterates: the update is
 * path dependent and the tangent approximate (\f$ 10^{-4} \f$ to \f$ 10^{-2} \f$
 * relative, more for Tresca). Under tangent_closest_point the orchestrator integrates the
 * row by closest-point projection instead (closest_point_return_mapping): the state is
 * rebuilt from the start values at every iterate (refresh_state) and the Newton and tangent
 * use the total derivatives of the composed map (closest_point_ingredients), which makes the
 * consistent operator exact for every criterion with a flow Hessian and every hardening law.
 * The cutting-plane members above are untouched by it. Under damage the mechanism works on
 * the effective stress (strain equivalence, see DamageMechanism).
 *
 * Closest-point linearisation. With \f$ \mathbf{a}_i = (\mathbf{a}_i^n +
 * \Delta p\,\mathbf{n})/(1 + D_i\Delta p) \f$ the backstress is affine in the normal,
 * \f$ \mathbf{X} = \mathbf{X}_0 + \gamma\mathbf{T}\mathbf{n} \f$
 * (KinematicHardening::backward_euler_factors), so at fixed \f$ (\boldsymbol{\sigma},
 * \Delta p) \f$ the normal solves \f$ \mathbf{G}(\mathbf{n}) = \mathbf{n} -
 * \boldsymbol{\eta}(\boldsymbol{\sigma} - \mathbf{X}_0 - \gamma\mathbf{T}\mathbf{n}) = 0 \f$,
 * Jacobian \f$ \mathbf{A} = \mathbf{I} + \gamma\mathbf{H}\mathbf{T} \f$ (one Newton step for
 * von Mises, whose normal is scale-invariant along \f$ \mathrm{dev}\,\boldsymbol{\xi} \f$).
 * Then \f$ \mathbf{D} = \mathbf{A}^{-1}\mathbf{H} \f$, \f$ \tilde{\mathbf{n}} = \mathbf{n} -
 * \gamma\mathbf{D}^{T}\mathbf{T}\mathbf{n} \f$, \f$ K = -\mathbf{n}\cdot\mathbf{T}\mathbf{A}^{-1}
 * \boldsymbol{\beta} - R'(p) \f$, \f$ d\mathbf{n}/d\Delta p = -\mathbf{A}^{-1}
 * \mathbf{H}\mathbf{T}\boldsymbol{\beta} \f$. Limits: \f$ \gamma = 0 \f$ gives the partial
 * forms (J2 + isotropic: closest-point and cutting-plane coincide); Prager gives
 * \f$ \tilde{\mathbf{n}} = \mathbf{n} \f$, \f$ d\mathbf{n}/d\Delta p = 0 \f$, \f$ K = -(C + R') \f$
 * but \f$ \mathbf{D} \ne \mathbf{H} \f$ (the partial Hessian is inexact as soon as
 * \f$ \gamma > 0 \f$). The exact operator is symmetric for isotropic and Prager
 * hardening, not for Armstrong-Frederick / Chaboche (dynamic recovery).
 */
class PlasticityMechanism final : public StrainMechanism {
private:
    YieldType yield_type_;             ///< Yield criterion type (for deferred configuration)
    std::unique_ptr<YieldCriterion> yield_;
    std::unique_ptr<IsotropicHardening> iso_hard_;
    std::unique_ptr<KinematicHardening> kin_hard_;
    double sigma_Y_;  ///< Initial yield stress
    size_t p_index_ = 0;   ///< ivc_ index of "p", resolved by register_variables()
    size_t EP_index_ = 0;  ///< ivc_ index of "EP", resolved by register_variables()

    // Per-iteration caches populated by compute_constraints, raw engineering Voigt;
    // dPhi_dsigma() and kappa() return const-refs to their one-column views.
    mutable arma::vec flow_dir_;                ///< ∂Φ/∂σ (associated flow direction)
    mutable arma::vec kappa_;                   ///< L_ref · flow_dir_
    mutable arma::mat dPhi_column_ = arma::mat(6, 1);    ///< flow_dir_ as the 6 x 1 of the interface
    mutable arma::mat kappa_column_ = arma::mat(6, 1);   ///< kappa_ as the 6 x 1 of the interface
    mutable arma::mat::fixed<6, 6> hessian_eng_;   ///< flow Hessian at the shifted stress (dLambda_dsigma)
    mutable arma::vec::fixed<6> hessian_at_;       ///< the stress hessian_eng_ was evaluated at
    mutable bool hessian_valid_ = false;           ///< hessian_eng_ holds the Hessian of the current state at hessian_at_
    bool implicit_update_ = false;                 ///< the last state update was the backward-Euler one of refresh_state (closest point), else the cutting-plane step of update()
    mutable arma::mat X_branches_;              ///< Per-branch backstresses (6 x n, stress Voigt), built
                                                ///< ONCE per FB iteration in
                                                ///< compute_constraints and reused
                                                ///< for sigma_eff + hardening_modulus
    mutable double H_total_{0.0};               ///< Hardening modulus (iso + kin, the rate sensitivity of the hardening included)
    size_t pdot_index_ = 0;                     ///< ivc_ index of "pdot" (rate-dependent hardening only)
    bool has_pdot_ = false;                     ///< the hardening is rate dependent: a "pdot" variable follows "EP" (and the back strains)
    double theta_start_ = 0.0;                  ///< temperatures and time step of the increment (set_increment)
    double theta_end_ = 0.0;
    double DTime_ = 0.0;
    /// R and its derivatives at p for the increment Dp = p - p_start: R_rate of the hardening (rate term included)
    void hardening_stress(double p, double Dp, double& R, double& dR_dp, double& dR_dDp) const;
    std::vector<ClosestPointIngredients> cpp_cache_{ClosestPointIngredients{}};  ///< built by refresh_state, read by closest_point_ingredients

    /// Total backstress (engineering stress Voigt), summed in the X_branches_ buffer it overwrites (no allocation)
    [[nodiscard]] arma::vec total_backstress_voigt() const;

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

    /**
     * @brief Constructor configured and registered on @p props (sigma_Y, criterion,
     * isotropic, then kinematic parameters, in configure() order) — the dedicated kernels
     * build their closest-point row in one line with it.
     */
    PlasticityMechanism(
        YieldType yield_type,
        IsoHardType iso_type,
        KinHardType kin_type,
        const arma::vec& props,
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

    /**
     * @brief Reads, from @p props at @p offset (advanced past them): \f$ \sigma_Y \f$, then the
     * yield criterion parameters (unless the criterion is already configured: von Mises and Tresca
     * need none and are set in the constructor), the isotropic and the kinematic hardening parameters.
     */
    void configure(const arma::vec& props, int& offset) override;
    void register_variables() override;

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
    // callers do not need to reach into kin_hard_->total_backstress().

    /**
     * @brief Equivalent stress of (sigma - X) under the configured yield criterion.
     *
     * Dispatches through YieldCriterion to criteria.cpp::Eq_stress
     * (Mises / Tresca / Drucker / Hill / DFA / Ani).
     */
    [[nodiscard]] double equivalent_stress(const arma::vec& sigma) const;

    /// tensor2-typed overload.
    [[nodiscard]] double equivalent_stress(const tensor2& sigma) const;

    /**
     * @brief Yield function value \f$ \Phi = \sigma_{eq}(\boldsymbol{\sigma} - \mathbf{X}) - R(p) - \sigma_Y \f$,
     *        \f$ R(p) \f$ the isotropic hardening stress (IsotropicHardening), \f$ \mathbf{X} \f$ the back stress.
     *
     * Phi <= 0 inside the elastic domain; the return-mapping drives Phi to 0
     * when the trial state lies outside.
     */
    [[nodiscard]] double yield_function(const arma::vec& sigma) const;

    [[nodiscard]] double yield_function(const tensor2& sigma) const;

    /// Flow direction dPhi/dsigma at (sigma - X), strain-typed in the tensor2 overload.
    [[nodiscard]] arma::vec flow_direction(const arma::vec& sigma) const;

    [[nodiscard]] tensor2 flow_direction(const tensor2& sigma) const;

    // ========== Constitutive Computations ==========

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

    [[nodiscard]] const arma::mat& dPhi_dsigma(
        const arma::vec& sigma) const override;

    [[nodiscard]] const arma::mat& kappa(
        const arma::vec& sigma,
        double DT,
        const arma::mat& L_ref) const override;

    /// Flow Hessian dLambda/dsigma at the shifted stress (sigma - X), via the
    /// yield criterion's analytic kernel (deta_stress / ddHill / ddDFA / ddAni).
    /// nullptr for criteria without one (Tresca, Drucker) -> continuum fallback.
    [[nodiscard]] const arma::mat* dLambda_dsigma(
        const arma::vec& sigma) const override;

    /// has_flow_hessian(): von Mises, Hill, DFA, anisotropic. Tresca and Drucker have
    /// none and degrade the UMAT to the cutting-plane loop under tangent_closest_point.
    [[nodiscard]] bool supports_closest_point() const override;

    /// Backward-Euler refresh from start values (CPP contract): p = p_n + dp,
    /// EP = EP_n + dp n, back-strains via the laws' closed forms, the n <-> X
    /// coupling by a 6x6 Newton on A = I + gamma H T, then the total derivatives
    /// of the class note into closest_point_ingredients(). False when the Newton
    /// fails to converge.
    bool refresh_state(
        const arma::vec& sigma,
        const arma::vec& Ds_total,
        int offset) override;

    [[nodiscard]] const std::vector<ClosestPointIngredients>* closest_point_ingredients() const override {
        return &cpp_cache_;
    }

    /// Flow-rule residual of the committed state (tensorial strain norms,
    /// flow direction at the end stress shifted by the backstress); see
    /// StrainMechanism::state_drift.
    [[nodiscard]] double state_drift(const arma::vec& sigma) const override;

    /// Base default when the criterion's flow direction is stable within an
    /// increment, infinity (opt-out) otherwise — the classification lives on
    /// YieldCriterion::has_stable_flow_direction.
    [[nodiscard]] double drift_tolerance() const override;

    /// 1.0 = 100% plastic strain in a single increment: the accumulated
    /// plastic multiplier is strain-scale, and the first FB runaway
    /// excursion (traced: 0.83 -> 433 -> 4e13) must be caught before it
    /// saturates the hardening state, whose garbage tangent would reseed
    /// the divergence on every solver retry.
    [[nodiscard]] double multiplier_cap() const override;

    /// Committed Dp = p - p_start from the internal variable (projected by
    /// update(); see StrainMechanism::committed_multiplier).
    [[nodiscard]] double committed_multiplier() const override;

    const arma::vec& inelastic_strain() const override;

    /**
     * @brief Applies the FB correction \f$ \Delta s \f$, PROJECTED onto \f$ p \ge p_{start} \f$.
     *
     * Negative per-iteration corrections are kept, so an overshoot can walk back, but the
     * multiplier never drops below its start-of-increment value. Without the projection a tiny
     * negative \f$ \Delta p \f$ is a spurious Fischer-Burmeister root whenever \f$ \Phi > 0 \f$ is
     * large (\f$ \phi(\Phi, \Delta p) \approx \Delta p \f$ for \f$ |\Delta p| \ll \Phi \f$), reached for
     * near-singular hardening slopes such as \f$ R = k p^m \f$, \f$ m < 1 \f$, at first yield.
     */
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
     * @brief Thermodynamic split of the plastic work of the increment, the trapezoids of the
     * dedicated kernels.
     *
     * - \f$ W_m^{ir} = \tfrac12 (R_n + R_{n+1})\,\Delta p \f$: the isotropic hardening work
     *   (\f$ F_p\,\dot p \f$, IsotropicHardening), stored and never given back;
     * - \f$ W_m^{r} = \sum_i \tfrac12 (\mathbf{X}_i^n + \mathbf{X}_i^{n+1}) : \Delta\mathbf{a}_i \f$: the
     *   energy stored in the back stresses (\f$ \mathbf{X}:\dot{\mathbf{a}} \f$,
     *   KinematicHardening::stored_increment), recovered on a reverse loading;
     * - \f$ W_m^d = \tfrac12 (\boldsymbol{\sigma}_n + \boldsymbol{\sigma}_{n+1}) : \Delta\boldsymbol{\varepsilon}^p
     *   - W_m^{ir} - W_m^{r} \f$: the dissipation \f$ \Delta\gamma_{loc} \f$ the heat of the thermomechanical
     *   face is built on.
     *
     * The elastic part of the recoverable work is not here: ModularUMAT closes
     * \f$ W_m^r = W_m - W_m^{ir} - W_m^d \f$ over the whole law, which holds both.
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
     * @brief Derivatives of the dissipation of compute_work (see StrainMechanism): the flow direction
     * \f$ \boldsymbol{\Lambda} \f$ of the last compute_constraints (the direction the update moved
     * along: the last cutting-plane iterate's, as the dedicated kernels; the committed one when the
     * orchestrator re-evaluates the constraints there) and the hardening forces at the committed state,
     * the turning of the direction over the increment through the flow Hessian at the committed
     * shifted EFFECTIVE stress \f$ \boldsymbol{\xi} = \boldsymbol{\sigma}_{eff} - \mathbf{X} \f$ (criteria that have
     * one, read from dLambda_dsigma() when it was evaluated there; held fixed otherwise), applied to
     * \f$ \boldsymbol{\sigma}_n - \mathbf{X}_n + \boldsymbol{\sigma}_{n+1} - \mathbf{X}_{n+1} \f$ (the stresses of the
     * dissipation; \f$ \mathbf{H}\boldsymbol{\xi}_{n+1} = 0 \f$) and contracted with the derivatives of the
     * effective stress, and, for the back stresses, their evolution along the flow at fixed direction:
     * \f$ \partial\mathbf{a}_{i,n+1}/\partial\Delta p = \boldsymbol{\Lambda} \f$ (Prager); with recovery
     * (Armstrong-Frederick, Chaboche) \f$ (\boldsymbol{\Lambda} - D_i \mathbf{a}_{i,n+1}) / (1 + D_i \Delta p) \f$,
     * the derivative of the backward-Euler update of the closest point (second order, 3e-7 on a
     * non-proportional probe), and the secant \f$ \Delta\mathbf{a}_i/\Delta p \f$ over the explicit
     * iterates of the cutting plane (first order, 4e-5 on the same probe; the derivative of the last iterate
     * is worse there). The \f$ \boldsymbol{\sigma}:\Delta\mathbf{E}^p \f$ term is on the stress the
     * dissipation is written on (the damaged one under a damage mechanism, strain equivalence).
     */
    void dissipation_sensitivity(const DissipationInputs& in, arma::vec& dDgamma_dE, double& dDgamma_dT) const override;

    void set_increment(double theta_start, double theta_end, double DTime) override;
    /// Resets the rate variable: an increment that stays elastic leaves no plastic strain rate
    void predict(const arma::vec& E_total_end, double DTime) override;
    bool predict_multipliers(const arma::vec& sigma_eff, const arma::mat& L, arma::vec& ds, int offset) const override;
};

} // namespace simcoon
