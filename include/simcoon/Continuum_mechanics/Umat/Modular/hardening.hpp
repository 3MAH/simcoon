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
 * @file hardening.hpp
 * @brief Hardening modules for modular UMAT.
 *
 * Provides isotropic and kinematic hardening models:
 * - Isotropic: Linear, Power-law, Voce, Combined Voce
 * - Kinematic: Prager (linear), Armstrong-Frederick, Chaboche (multi-AF)
 *
 * @see PlasticityMechanism, ModularUMAT
 * @version 1.0
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/internal_variable_collection.hpp>

namespace simcoon {

// ============================================================================
// ISOTROPIC HARDENING
// ============================================================================

/**
 * @brief Types of isotropic hardening
 */
enum class IsoHardType {
    NONE = 0,           ///< No isotropic hardening
    LINEAR = 1,         ///< Linear: R = H * p
    POWER_LAW = 2,      ///< Power law: R = k * p^m
    VOCE = 3,           ///< Voce saturation: R = Q * (1 - exp(-b*p))
    COMBINED_VOCE = 4,  ///< Sum of N Voce terms
    JOHNSON_COOK = 5    ///< Johnson-Cook: (sigma_Y + B p^n)(1 + C ln(pdot/pdot0))(1 - T*^m) - sigma_Y, rate and temperature dependent
};

/**
 * @brief Base class for isotropic hardening: the hardening stress \f$ R(p) \f$ of the accumulated
 *        plastic strain \f$ p \f$.
 *
 * \f$ R \f$ is the growth of the yield stress with \f$ p \f$: the criterion of PlasticityMechanism is
 * \f$ \Phi = \sigma_{eq}(\boldsymbol{\sigma} - \mathbf{X}) - R(p) - \sigma_Y \f$, so the elastic domain has the
 * radius \f$ \sigma_Y + R(p) \f$, with \f$ R(0) = 0 \f$. It is a stress (the unit of \f$ \sigma_Y \f$), not a
 * strain. Thermodynamically it is the derivative of the hardening function \f$ F(p) \f$ of the irrecoverable
 * free energy, \f$ R = F_p = \partial F/\partial p \f$, i.e. minus the force conjugate to \f$ p \f$
 * (\f$ A_p = -F_p \f$, the notation of the theory manual and of the dedicated kernels). Its work
 * \f$ \int R\,\mathrm{d}p \f$ is stored and not given back on unloading or reverse loading: the
 * irrecoverable work \f$ W_m^{ir} \f$ of PlasticityMechanism::compute_work, excluded from the dissipation.
 *
 * Not to be confused with the heat tangents \f$ \mathbf{R}^\varepsilon = \partial r/\partial\boldsymbol{\varepsilon} \f$
 * and \f$ R^\theta = \partial r/\partial\theta \f$ of the thermomechanical face (@c drdE, @c drdT).
 */
class IsotropicHardening {
public:
    virtual ~IsotropicHardening() = default;

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Hardening stress \f$ R(p) = F_p \f$: the growth of the yield stress, \f$ \sigma_Y + R(p) \f$
     * @param p Accumulated plastic strain
     * @return Hardening stress (unit of stress)
     */
    [[nodiscard]] virtual double R(double p) const = 0;

    /**
     * @brief Hardening modulus \f$ \mathrm{d}R/\mathrm{d}p = F_{pp} \f$, the plastic modulus of a uniaxial curve
     * @param p Accumulated plastic strain
     * @return Derivative of the hardening stress (unit of stress)
     */
    [[nodiscard]] virtual double dR_dp(double p) const = 0;

    /**
     * @brief \f$ R(p) \f$ and \f$ dR/dp \f$ together, for the hot loop: one evaluation of the
     *        transcendental they share (the power or exponential); the default calls the two.
     * @param p Accumulated plastic strain
     * @param[out] R Hardening stress
     * @param[out] dR Its derivative
     */
    virtual void R_and_dR(double p, double& R, double& dR) const { R = this->R(p); dR = dR_dp(p); }

    /**
     * @brief True for a yield stress that depends on the plastic strain rate and on the temperature
     *        (JohnsonCookHardening). The mechanism then reads R_rate() in its return mapping, carries a
     *        @c pdot internal variable, and the thermomechanical face of MODUL does not serve it yet.
     */
    [[nodiscard]] virtual bool rate_dependent() const noexcept { return false; }
    /**
     * @brief The initial yield stress \f$ \sigma_Y \f$ of the mechanism, for the laws that scale it
     *        (Johnson-Cook: \f$ (\sigma_Y + B p^n) f_{\dot p} f_\theta \f$). Default: unused.
     */
    virtual void set_yield_stress(double sigma_Y) { (void)sigma_Y; }
    /**
     * @brief The increment a rate- and temperature-dependent law reads: the temperature at both ends
     *        and the time step (\f$ \Delta t \le \iota \f$ is the rate-independent evaluation, the
     *        zero-time tangent probe). Default: nothing.
     */
    virtual void set_increment(double theta_start, double theta_end, double DTime) {
        (void)theta_start; (void)theta_end; (void)DTime;
    }
    /**
     * @brief \f$ R(p, \dot p, \theta) \f$ of the implicit rate \f$ \dot p = \Delta p/\Delta t \f$ and its partial
     *        derivatives in \f$ p \f$ and in \f$ \Delta p \f$ (the second enters the multiplier Jacobian
     *        \f$ K = -(\partial R/\partial p + \partial R/\partial\Delta p) \f$). Default: R_and_dR, no rate term.
     * @param p accumulated plastic strain at the end of the increment
     * @param Dp its increment
     * @param[out] R hardening stress
     * @param[out] dR_dp \f$ \partial R/\partial p \f$ at fixed rate
     * @param[out] dR_dDp \f$ \partial R/\partial\Delta p \f$ at fixed \f$ p \f$ (the rate sensitivity)
     */
    virtual void R_rate(double p, double Dp, double& R, double& dR_dp, double& dR_dDp) const {
        (void)Dp;
        R_and_dR(p, R, dR_dp);
        dR_dDp = 0.0;
    }
    /**
     * @brief The force stored by the hardening, \f$ F_p = -A_p \f$ of the energy split
     *        (PlasticityMechanism::compute_work): \f$ R(p) \f$ for the rate-independent laws; the
     *        isothermal reference-rate part \f$ f_\theta B p^n \f$ for Johnson-Cook, whose rate and
     *        thermal multipliers of the yield stress are dissipated.
     */
    [[nodiscard]] virtual double stored_force(double p, double theta) const { (void)theta; return R(p); }
    /**
     * @brief Starting multiplier of the cutting plane for a law whose \f$ \Phi(\Delta p) \f$ has a kink
     *        (the clamped rate factor of Johnson-Cook at \f$ \dot p = \dot p_0 \f$): the Newton started at
     *        \f$ \Delta p = 0 \f$ with the quasi-static slope overshoots it and cycles. Returns the
     *        multiplier to start from (0: none) given the trial state, \f$ \bar\sigma^{tr} \f$ and
     *        \f$ \boldsymbol{\Lambda}:\mathbf{L}:\boldsymbol{\Lambda} \f$. Default: 0.
     */
    virtual double multiplier_predictor(double p, double sigma_eq_trial, double nLn) const {
        (void)p; (void)sigma_eq_trial; (void)nLn;
        return 0.0;
    }

    /**
     * @brief Get the hardening type
     * @return Hardening type enum
     */
    [[nodiscard]] virtual IsoHardType type() const = 0;

    /**
     * @brief Create an isotropic hardening instance
     * @param type Hardening type
     * @param N Number of terms (for COMBINED_VOCE)
     * @return Unique pointer to hardening instance
     */
    static std::unique_ptr<IsotropicHardening> create(IsoHardType type, int N = 1);

    /**
     * @brief Get number of props consumed by this hardening type
     * @param type The hardening type
     * @param N Number of terms (for COMBINED_VOCE)
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(IsoHardType type, int N = 1) {
        switch (type) {
            case IsoHardType::NONE:          return 0;
            case IsoHardType::LINEAR:        return 1;
            case IsoHardType::POWER_LAW:     return 2;
            case IsoHardType::VOCE:          return 2;
            case IsoHardType::COMBINED_VOCE: return 2 * N;
            case IsoHardType::JOHNSON_COOK:  return 7;
        }
        return 0;
    }
};

/**
 * @brief No isotropic hardening
 */
class NoIsotropicHardening final : public IsotropicHardening {
public:
    void configure(const arma::vec& props, int& offset) override {}
    double R(double p) const override { return 0.0; }
    double dR_dp(double p) const override { return 0.0; }
    IsoHardType type() const override { return IsoHardType::NONE; }
};

/**
 * @brief Linear isotropic hardening: R = H * p
 */
class LinearHardening final : public IsotropicHardening {
private:
    double H_;  ///< Hardening modulus
public:
    LinearHardening() : H_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override { return H_ * p; }
    double dR_dp(double p) const override { return H_; }
    IsoHardType type() const override { return IsoHardType::LINEAR; }
};

/**
 * @brief Power-law isotropic hardening: \f$ R = k\,p^m \f$
 *
 * For \f$ m < 1 \f$ the exact modulus \f$ dR/dp = k\,m\,p^{m-1} \f$ diverges
 * at the onset \f$ p \to 0^+ \f$, which destabilizes the return-map Newton.
 * The law is therefore regularized below a small cutoff \f$ p_{reg} = 10^{-6} \f$
 * by a C\f$^1\f$ quadratic blend anchored at \f$ R(0)=0 \f$ and matched in
 * value and slope at \f$ p_{reg} \f$ — the same tangent-continuation idiom as
 * the SMA volume-fraction barriers (lagrange_pow_0/1). The law is exact for
 * \f$ p \ge p_{reg} \f$ and unmodified when \f$ m \ge 1 \f$.
 */
class PowerLawHardening final : public IsotropicHardening {
private:
    double k_;  ///< Hardening coefficient
    double m_;  ///< Hardening exponent
    double a_reg_;  ///< Onset-blend linear coefficient (precomputed, m < 1 only)
    double b_reg_;  ///< Onset-blend quadratic coefficient (precomputed, m < 1 only)
public:
    PowerLawHardening() : k_(0.0), m_(1.0), a_reg_(0.0), b_reg_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    void R_and_dR(double p, double& R, double& dR) const override;
    IsoHardType type() const override { return IsoHardType::POWER_LAW; }
};

/**
 * @brief Voce saturation hardening: R = Q * (1 - exp(-b*p))
 */
class VoceHardening final : public IsotropicHardening {
private:
    double Q_;  ///< Saturation stress
    double b_;  ///< Hardening rate
public:
    VoceHardening() : Q_(0.0), b_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    void R_and_dR(double p, double& R, double& dR) const override;
    IsoHardType type() const override { return IsoHardType::VOCE; }
};

/**
 * @brief Combined Voce hardening: R = sum_i Q_i * (1 - exp(-b_i*p))
 */
class CombinedVoceHardening final : public IsotropicHardening {
private:
    int N_;           ///< Number of Voce terms
    arma::vec Q_;     ///< Saturation stresses
    arma::vec b_;     ///< Hardening rates
public:
    explicit CombinedVoceHardening(int N = 1) : N_(N), Q_(arma::zeros(N)), b_(arma::zeros(N)) {}
    void configure(const arma::vec& props, int& offset) override;
    double R(double p) const override;
    double dR_dp(double p) const override;
    void R_and_dR(double p, double& R, double& dR) const override;
    IsoHardType type() const override { return IsoHardType::COMBINED_VOCE; }
    int num_terms() const { return N_; }
};

/**
 * @brief Johnson-Cook yield stress (Johnson & Cook, 1983) as an isotropic hardening of the
 *        plastic strain rate and the temperature.
 *
 * \f[
 *   \sigma_Y + R(p, \dot p, \theta) = \left(\sigma_Y + B p^{n}\right)
 *   \left(1 + C \ln\frac{\dot p}{\dot p_0}\right)\left(1 - \theta^{*m}\right),
 *   \qquad \theta^* = \frac{\theta - \theta_{ref}}{\theta_{melt} - \theta_{ref}},
 * \f]
 * props: \f$ B, n, C, \dot p_0, m, \theta_{ref}, \theta_{melt} \f$; \f$ \sigma_Y \f$ (the A of the law) is
 * the mechanism's, handed by set_yield_stress(). The same regularizations as the EPJCK kernel
 * (plastic_johnson_cook_ccp.hpp): \f$ B p^n \f$ through PowerLawHardening (C1 onset blend for
 * \f$ n < 1 \f$), the rate factor clamped at 1 below \f$ \dot p_0 \f$ (the rate is the implicit
 * \f$ \Delta p/\Delta t \f$; \f$ \Delta t \le \iota \f$ is the quasi-static form), \f$ \theta^* \f$ clamped to
 * \f$ [0, 1) \f$ and \f$ \theta^{*m} \f$ continued below \f$ 10^{-3} \f$ by johnson_cook_softening(). The
 * kink of the clamp is met as in the kernel: multiplier_predictor() decides the branch from
 * \f$ \Phi \f$ at \f$ \Delta p = \dot p_0 \Delta t \f$ and starts the cutting plane there on the
 * logarithmic branch, whose slope \f$ C/(\dot p\,\Delta t) \f$ then stays in the Jacobian.
 * Energy: stored_force() is \f$ f_\theta B p^n \f$ (\f$ A_p = -f_\theta B p^n \f$: the rate and thermal
 * multipliers of the yield stress are dissipated, the split of the kernel).
 */
class JohnsonCookHardening final : public IsotropicHardening {
private:
    double sigma_Y_ = 0.0;      ///< A: the mechanism's initial yield stress
    PowerLawHardening power_;   ///< B p^n (C1 onset blend)
    double C_ = 0.0;            ///< rate sensitivity
    double pdot0_ = 1.0;        ///< reference plastic strain rate
    double m_ = 1.0;            ///< thermal softening exponent
    double T_ref_ = 0.0;        ///< reference temperature
    double T_melt_ = 1.0;       ///< melting temperature
    double DTime_ = 0.0;        ///< time increment (0: rate-independent evaluation)
    double f_start_ = 1.0;      ///< thermal factor 1 - theta*^m at the start temperature
    double f_end_ = 1.0;        ///< the same at the end temperature
    mutable bool log_branch_ = false;   ///< the cutting plane runs on the logarithmic branch (multiplier_predictor)
    [[nodiscard]] double thermal_factor(double theta) const;
public:
    JohnsonCookHardening() = default;
    void configure(const arma::vec& props, int& offset) override;
    bool rate_dependent() const noexcept override { return true; }
    void set_yield_stress(double sigma_Y) override { sigma_Y_ = sigma_Y; }
    void set_increment(double theta_start, double theta_end, double DTime) override;
    /// Quasi-static (rate factor 1) hardening stress at the end temperature
    double R(double p) const override;
    double dR_dp(double p) const override;
    void R_rate(double p, double Dp, double& R, double& dR_dp, double& dR_dDp) const override;
    double stored_force(double p, double theta) const override;
    double multiplier_predictor(double p, double sigma_eq_trial, double nLn) const override;
    IsoHardType type() const override { return IsoHardType::JOHNSON_COOK; }
};
// ============================================================================
// KINEMATIC HARDENING
// ============================================================================

/**
 * @brief Types of kinematic hardening
 */
enum class KinHardType {
    NONE = 0,                   ///< No kinematic hardening
    PRAGER = 1,                 ///< Linear Prager: X = (2/3)*C*alpha
    ARMSTRONG_FREDERICK = 2,    ///< Single AF: dX = (2/3)*C*dep - D*X*dp
    CHABOCHE = 3                ///< Multiple AF terms
};

/**
 * @brief Base class for kinematic hardening
 *
 * Hardening models are components of a PlasticityMechanism: they register
 * their back-strains into the OWNING MECHANISM's internal-variable
 * collection (passed explicitly), they do not own state themselves. Since
 * every mechanism owns a separate collection, the bare keys ("a", "a_0", ...)
 * can never collide across mechanisms.
 */
class KinematicHardening {
public:
    virtual ~KinematicHardening() = default;

    /**
     * @brief Configure from props array
     * @param props Material properties vector
     * @param offset Current offset in props (will be updated)
     */
    virtual void configure(const arma::vec& props, int& offset) = 0;

    /**
     * @brief Register internal variables for this hardening model
     * @param ivc The owning mechanism's internal variable collection; every other member taking a
     *        collection must get this same one (read by the indices resolved here, debug-checked).
     */
    virtual void register_variables(InternalVariableCollection& ivc) = 0;

    /**
     * @brief Compute the per-branch backstresses \f$ \mathbf{X}_i = \tfrac{2}{3}
     * C_i \mathbf{a}_i \f$, one column per branch, engineering stress Voigt
     * (the back-strain carries the doubled shear, the back stress does not).
     *
     * This is the ONE place the backstresses are built from the back-strains.
     * The hot caller (PlasticityMechanism::compute_constraints) invokes it once
     * per FB iteration and reuses the result for the shifted stress
     * \f$ \boldsymbol{\sigma} - \sum_i \mathbf{X}_i \f$ AND for
     * hardening_modulus(n, X_i) — the branch backstresses are only valid within
     * that iteration (the back-strains move in update()).
     *
     * @param ivc Internal variable collection
     * @param X_i Output: 6 x num_backstresses(), branch order fixed
     */
    void compute_backstresses(const InternalVariableCollection& ivc, arma::mat& X_i) const;

    /**
     * @brief Kinematic hardening modulus from PRECOMPUTED branch backstresses
     * (same order as compute_backstresses)
     * @param n Plastic flow direction, engineering strain Voigt (doubled shear)
     * @param X_i Branch backstresses from compute_backstresses()
     * @return Exact slope \f$ \sum_i \tfrac{2}{3} C_i \langle \mathbf{n},
     *         \mathbf{n} \rangle - D_i (\mathbf{X}_i : \mathbf{n}) \f$, the contractions
     *         on the engineering components with their shear weights
     */
    virtual double hardening_modulus(const arma::vec& n,
                                     const arma::mat& X_i) const = 0;

    // ---- Convenience forms (cold paths: work computation, tests). Hot loops use
    // the two virtuals above directly. ----

    /// Total backstress X = sum(X_i), stress-typed.
    [[nodiscard]] tensor2 total_backstress(const InternalVariableCollection& ivc) const {
        arma::mat X_i;
        compute_backstresses(ivc, X_i);
        arma::vec::fixed<6> X(arma::fill::zeros);
        for (arma::uword i = 0; i < X_i.n_cols; ++i) {
            X += X_i.col(i);
        }
        return stress(X);
    }

    /// Hardening modulus computed from the collection directly.
    [[nodiscard]] double hardening_modulus(const tensor2& n,
                                           const InternalVariableCollection& ivc) const {
        arma::mat X_i;
        compute_backstresses(ivc, X_i);
        return hardening_modulus(n.voigt(), X_i);
    }

    /**
     * @brief Update internal variables given plastic multiplier increment
     * @param dp Plastic multiplier increment
     * @param n Plastic flow direction, engineering strain Voigt
     * @param ivc Internal variable collection
     */
    virtual void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) = 0;

    /**
     * @brief Backward-Euler refresh of the back-strains FROM THE START VALUES
     * (closest-point / tangent_mode == tangent_closest_point contract).
     *
     * Unlike update() (incremental, frozen flow), this solves the implicit
     * update \f$ \mathbf{a} = \mathbf{a}_n + \Delta p\,
     * (\mathbf{n} - D\,\mathbf{a}) \f$ in closed form
     * \f$ \mathbf{a} = (\mathbf{a}_n + \Delta p\,\mathbf{n})
     * / (1 + D\,\Delta p) \f$ so the CPP Newton can re-evaluate the state at
     * every iterate. @p n is the flow normal, engineering strain Voigt. Default: no
     * state, nothing to do.
     */
    virtual void refresh_state(double dp, const arma::vec& n,
                               InternalVariableCollection& ivc) const {
        (void)dp; (void)n; (void)ivc;
    }

    /**
     * @brief Linearisation of the backward-Euler branch maps at the refreshed
     * state (closest-point integrator).
     *
     * With \f$ \mathbf{a}_i = (\mathbf{a}_i^n + \Delta p\,\mathbf{n})
     * / (1 + D_i \Delta p) \f$ and \f$ \mathbf{X} = \sum_i \tfrac{2}{3} C_i\,
     * \mathbf{T}\mathbf{a}_i \f$ (\f$ \mathbf{T} \f$: engineering to tensorial
     * shear), the total backstress is affine in the flow normal,
     * \f$ \mathbf{X} = \mathbf{X}_0(\Delta p) + \gamma\,\mathbf{T}\mathbf{n} \f$, with
     * \f[ \gamma = \Delta p \sum_i \frac{\tfrac{2}{3} C_i}{1 + D_i\Delta p}, \qquad
     *      \boldsymbol{\beta} = \sum_i \frac{\tfrac{2}{3} C_i\,(\mathbf{n} - D_i\mathbf{a}_i)}{1 + D_i\Delta p}
     *      = \frac{\partial \mathbf{X}}{\partial \Delta p}\Big|_{\mathbf{n}}\ (\text{before } \mathbf{T}). \f]
     * @p alpha_i are the CURRENT (refreshed) back-strains. Default (no state):
     * \f$ \gamma = 0 \f$, \f$ \boldsymbol{\beta} = \mathbf{0} \f$.
     *
     * @param dp     Total plastic multiplier of the increment
     * @param n      Flow normal (engineering Voigt, 6)
     * @param ivc    Collection holding the refreshed back-strains
     * @param gamma  Output \f$ \gamma \f$ (stress / strain)
     * @param beta   Output \f$ \boldsymbol{\beta} \f$ (6, engineering-strain-typed like n)
     */
    virtual void backward_euler_factors(double dp, const arma::vec& n,
                                        const InternalVariableCollection& ivc,
                                        double& gamma, arma::vec& beta) const {
        (void)dp; (void)n; (void)ivc;
        gamma = 0.0;
        beta = arma::zeros(6);
    }

    /**
     * @brief Get the hardening type
     * @return Hardening type enum
     */
    [[nodiscard]] virtual KinHardType type() const = 0;

    /**
     * @brief Get number of backstress terms
     * @return Number of backstress tensors
     */
    [[nodiscard]] virtual int num_backstresses() const = 0;

    /**
     * @brief Kinematic modulus \f$ C_i \f$ of branch @p i (\f$ \mathbf{X}_i = \tfrac23 C_i \mathbf{a}_i \f$); 0 without branches
     */
    [[nodiscard]] virtual double branch_modulus(int i) const { (void)i; return 0.0; }

    /**
     * @brief Dynamic recovery parameter \f$ D_i \f$ of branch @p i (\f$ \mathrm{d}\mathbf{a}_i = (\boldsymbol{\Lambda} - D_i \mathbf{a}_i)\,\mathrm{d}p \f$); 0 for Prager and without branches
     */
    [[nodiscard]] virtual double recovery_parameter(int i) const { (void)i; return 0.0; }

    /**
     * @brief The back-strain variable \f$ \mathbf{a}_i \f$ of branch @p i in @p ivc; nullptr without branches
     */
    [[nodiscard]] virtual const InternalVariable* back_strain_variable(const InternalVariableCollection& ivc, int i) const {
        (void)ivc; (void)i;
        return nullptr;
    }

    /**
     * @brief Back stress of branch @p i for a back-strain @p a (the linear map of compute_backstresses),
     *        engineering stress Voigt
     */
    [[nodiscard]] arma::vec backstress_voigt(int i, const arma::vec& a) const;

    /**
     * @brief Energy stored in the back stresses over the increment, per unit volume,
     * \f$ \sum_i \tfrac12 (\mathbf{X}_i^n + \mathbf{X}_i^{n+1}) : \Delta\mathbf{a}_i \f$: the trapezoid in the
     * back-strains, exact for the quadratic potentials \f$ \tfrac13 C_i\,\mathbf{a}_i:\mathbf{a}_i \f$.
     * This energy is recoverable: a reverse loading (compression after tension) unloads the back stress
     * and gives it back. It is the recoverable work of PlasticityMechanism::compute_work (part of
     * \f$ W_m^r \f$); with the isotropic hardening work, what the dissipation excludes.
     * @param ivc the collection holding the back-strains (start and current values)
     * @return the stored energy increment
     */
    [[nodiscard]] double stored_increment(const InternalVariableCollection& ivc) const;

    /**
     * @brief Create a kinematic hardening instance
     * @param type Hardening type
     * @param N Number of backstress terms (for CHABOCHE)
     * @return Unique pointer to hardening instance
     */
    static std::unique_ptr<KinematicHardening> create(KinHardType type, int N = 1);

    /**
     * @brief Get number of props consumed by this hardening type
     * @param type The hardening type
     * @param N Number of backstress terms (for CHABOCHE)
     * @return Number of properties required
     */
    [[nodiscard]] static constexpr int props_count(KinHardType type, int N = 1) {
        switch (type) {
            case KinHardType::NONE:                 return 0;
            case KinHardType::PRAGER:               return 1;
            case KinHardType::ARMSTRONG_FREDERICK:   return 2;
            case KinHardType::CHABOCHE:             return 2 * N;
        }
        return 0;
    }
};

/**
 * @brief No kinematic hardening
 */
class NoKinematicHardening final : public KinematicHardening {
public:
    void configure(const arma::vec& props, int& offset) override {}
    void register_variables(InternalVariableCollection& ivc) override {}
    double hardening_modulus(const arma::vec& n,
                             const arma::mat& X_i) const override {
        return 0.0;
    }
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override {}
    KinHardType type() const override { return KinHardType::NONE; }
    int num_backstresses() const override { return 0; }
};

/**
 * @brief Linear Prager kinematic hardening: \f$ \mathbf{X} = \tfrac{2}{3} C \mathbf{a} \f$
 *
 * The back-strain \f$ \mathbf{a} \f$ is the thermodynamic internal
 * variable (stored strain-like, factor-2 shear Voigt); the backstress
 * \f$ \mathbf{X} = \tfrac{2}{3} C \mathbf{a} \f$ is its conjugate
 * generalised force and is never independently stored. The same convention
 * holds for Armstrong-Frederick and Chaboche below.
 */
class PragerHardening final : public KinematicHardening {
private:
    double C_;  ///< Kinematic hardening modulus
    std::string a_key_{"a"};  ///< IVC key for the back-strain α
    size_t a_index_ = 0;      ///< its index in the collection, resolved by register_variables
public:
    PragerHardening() : C_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    double hardening_modulus(const arma::vec& n,
                             const arma::mat& X_i) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    void refresh_state(double dp, const arma::vec& n, InternalVariableCollection& ivc) const override;
    void backward_euler_factors(double dp, const arma::vec& n, const InternalVariableCollection& ivc,
                                double& gamma, arma::vec& beta) const override;
    KinHardType type() const override { return KinHardType::PRAGER; }
    int num_backstresses() const override { return 1; }
    double branch_modulus(int) const override { return C_; }
    const InternalVariable* back_strain_variable(const InternalVariableCollection& ivc, int) const override { return &ivc.at(a_index_); }
};

/**
 * @brief Armstrong-Frederick kinematic hardening: dX = (2/3)*C*dep - D*X*dp
 */
class ArmstrongFrederickHardening final : public KinematicHardening {
private:
    double C_;  ///< Hardening parameter
    double D_;  ///< Dynamic recovery parameter
    std::string a_key_{"a"};  ///< IVC key for the back-strain α
    size_t a_index_ = 0;      ///< its index in the collection, resolved by register_variables
public:
    ArmstrongFrederickHardening() : C_(0.0), D_(0.0) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    double hardening_modulus(const arma::vec& n,
                             const arma::mat& X_i) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    void refresh_state(double dp, const arma::vec& n, InternalVariableCollection& ivc) const override;
    void backward_euler_factors(double dp, const arma::vec& n, const InternalVariableCollection& ivc,
                                double& gamma, arma::vec& beta) const override;
    KinHardType type() const override { return KinHardType::ARMSTRONG_FREDERICK; }
    int num_backstresses() const override { return 1; }
    double branch_modulus(int) const override { return C_; }
    double recovery_parameter(int) const override { return D_; }
    const InternalVariable* back_strain_variable(const InternalVariableCollection& ivc, int) const override { return &ivc.at(a_index_); }
};

/**
 * @brief Chaboche kinematic hardening: dX_i = (2/3)*C_i*dep - D_i*X_i*dp
 *
 * Multiple Armstrong-Frederick terms for capturing different hardening time scales.
 */
class ChabocheHardening final : public KinematicHardening {
private:
    int N_;                          ///< Number of backstress terms
    arma::vec C_;                    ///< Hardening parameters
    arma::vec D_;                    ///< Dynamic recovery parameters
    std::vector<std::string> a_keys_; ///< IVC keys for back-strains α_i ("a_0", "a_1", ...)
    std::vector<size_t> a_index_;     ///< their indices in the collection, resolved by register_variables
public:
    explicit ChabocheHardening(int N = 1) : N_(N), C_(arma::zeros(N)), D_(arma::zeros(N)) {}
    void configure(const arma::vec& props, int& offset) override;
    void register_variables(InternalVariableCollection& ivc) override;
    double hardening_modulus(const arma::vec& n,
                             const arma::mat& X_i) const override;
    void update(double dp, const arma::vec& n, InternalVariableCollection& ivc) override;
    void refresh_state(double dp, const arma::vec& n, InternalVariableCollection& ivc) const override;
    void backward_euler_factors(double dp, const arma::vec& n, const InternalVariableCollection& ivc,
                                double& gamma, arma::vec& beta) const override;
    KinHardType type() const override { return KinHardType::CHABOCHE; }
    int num_backstresses() const override { return N_; }
    double branch_modulus(int i) const override { return C_(i); }
    double recovery_parameter(int i) const override { return D_(i); }
    const InternalVariable* back_strain_variable(const InternalVariableCollection& ivc, int i) const override { return &ivc.at(a_index_[static_cast<size_t>(i)]); }
};

} // namespace simcoon
