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
 * @file damage_mechanism.cpp
 * @brief Implementation of DamageMechanism class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/damage_mechanism.hpp>
#include <simcoon/parameter.hpp>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace simcoon {

// ========== Constructor ==========

DamageMechanism::DamageMechanism(DamageType type)
    : StrainMechanism("damage")
    , damage_type_(type)
    , Y_0_(0.0)
    , Y_c_(1.0)
    , D_c_(0.99)
    , A_(1.0)
    , n_(1.0)
    , Y_current_(0.0)
    , D_current_(0.0)
    , dD_dY_(0.0)
    , M_cached_(arma::zeros(6, 6))
    , M_cached_valid_(false)
{
}

// ========== Configuration ==========

void DamageMechanism::configure(const arma::vec& props, int& offset) {
    // Props layout for damage:
    // props[offset]: damage_type (if not set in constructor)
    // props[offset+1]: Y_0 (damage threshold)
    // props[offset+2]: Y_c (critical damage driving force)
    // props[offset+3]: D_c (critical damage value, optional, default 0.99)
    // props[offset+4]: A or n (damage parameter, depending on type)

    // Read parameters
    Y_0_ = props(offset);
    Y_c_ = props(offset + 1);
    offset += 2;

    // Read type-specific parameters
    switch (damage_type_) {
        case DamageType::LINEAR:
            // No additional parameters needed
            // D = (Y - Y_0) / (Y_c - Y_0)
            break;

        case DamageType::EXPONENTIAL:
            // A: exponential rate parameter
            A_ = props(offset);
            offset += 1;
            break;

        case DamageType::POWER_LAW:
            // n: power law exponent
            n_ = props(offset);
            offset += 1;
            break;

        case DamageType::WEIBULL:
            // A: scale parameter, n: shape parameter
            A_ = props(offset);
            n_ = props(offset + 1);
            offset += 2;
            break;
    }

    // Validate parameters
    if (Y_c_ <= Y_0_) {
        throw std::runtime_error("DamageMechanism: Y_c must be greater than Y_0");
    }
}

void DamageMechanism::register_variables() {
    ivc_.add_scalar("D",     0.0);
    ivc_.add_scalar("Y_max", 0.0);
}

// ========== Constitutive Computations ==========

double DamageMechanism::compute_driving_force(const arma::vec& sigma, const arma::mat& /*S*/) const {
    // Y = ½ σ : M : σ — strain-energy release rate. Uses the cached
    // compliance-typed tensor4 (set in compute_constraints alongside the
    // arma compliance), so no per-call 6×6 ctor copy.
    return 0.5 * arma::dot(sigma, (M_cached_t_ * stress(sigma)).to_arma_voigt());
}

double DamageMechanism::compute_damage(double Y, double Y_max) const {
    // Use maximum of current and historical driving force
    double Y_eff = std::max(Y, Y_max);

    // No damage below threshold
    if (Y_eff <= Y_0_) {
        return 0.0;
    }

    double D = 0.0;

    switch (damage_type_) {
        case DamageType::LINEAR:
            D = (Y_eff - Y_0_) / (Y_c_ - Y_0_);
            break;

        case DamageType::EXPONENTIAL:
            D = 1.0 - std::exp(-A_ * (Y_eff - Y_0_) / (Y_c_ - Y_0_));
            break;

        case DamageType::POWER_LAW:
            D = std::pow((Y_eff - Y_0_) / (Y_c_ - Y_0_), n_);
            break;

        case DamageType::WEIBULL:
            D = 1.0 - std::exp(-std::pow((Y_eff - Y_0_) / A_, n_));
            break;
    }

    // Limit damage to critical value
    return std::min(D, D_c_);
}

double DamageMechanism::get_damage() const {
    return ivc_.get("D").scalar();
}

void DamageMechanism::compute_constraints(
    const arma::vec& sigma,
    const arma::vec& /*E_total*/,
    const arma::mat& L,
    double /*DTime*/,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    Phi.set_size(1);
    Y_crit.set_size(1);

    // Get current damage and history
    D_current_ = ivc_.get("D").scalar();
    double Y_max = ivc_.get("Y_max").scalar();

    // Compute and cache compliance (raw + typed) on first use.
    if (!M_cached_valid_) {
        M_cached_ = arma::inv(L);
        M_cached_t_ = tensor4(M_cached_, Tensor4Type::compliance);
        M_cached_valid_ = true;
    }

    // Compute current driving force
    Y_current_ = compute_driving_force(sigma, M_cached_);

    const double Y_eff = std::max(Y_current_, Y_max);

    // History-type constraint: Phi = Y - Y_max. Damage is integrated
    // EXPLICITLY (update() sets Y_max = max(Y, Y_max) and D = f(Y_max)), not
    // through the FB multiplier — so after update() Phi self-satisfies
    // (Y == Y_max under loading, Phi <= 0 under unloading). The FB row exists
    // only to carry damage into the coupled convergence check.
    Phi(0) = Y_current_ - Y_max;

    // Critical value for convergence
    Y_crit(0) = std::max(Y_0_, 1e-6);

    // Compute derivative dD/dY for tangent
    if (Y_eff > Y_0_ && D_current_ < D_c_) {
        switch (damage_type_) {
            case DamageType::LINEAR:
                dD_dY_ = 1.0 / (Y_c_ - Y_0_);
                break;

            case DamageType::EXPONENTIAL:
                dD_dY_ = (A_ / (Y_c_ - Y_0_)) * std::exp(-A_ * (Y_eff - Y_0_) / (Y_c_ - Y_0_));
                break;

            case DamageType::POWER_LAW:
                dD_dY_ = (n_ / (Y_c_ - Y_0_)) * std::pow((Y_eff - Y_0_) / (Y_c_ - Y_0_), n_ - 1.0);
                break;

            case DamageType::WEIBULL:
                dD_dY_ = (n_ / A_) * std::pow((Y_eff - Y_0_) / A_, n_ - 1.0) *
                         std::exp(-std::pow((Y_eff - Y_0_) / A_, n_));
                break;
        }
    } else {
        dD_dY_ = 0.0;
    }
}
void DamageMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    arma::mat& B,
    int row_offset
) const {
    // Unit diagonal for the history-type damage row. Phi = Y - Y_max is
    // integrated explicitly (see compute_constraints / update), so the damage
    // multiplier is not solved implicitly; a unit slope keeps the FB system
    // well-conditioned without steering the (self-satisfying) damage row.
    // Cross-mechanism coupling (plasticity <-> damage) still flows through the
    // off-diagonal dPhi_dsigma . kappa terms assembled by the orchestrator.
    B(row_offset, row_offset) = 1.0;
}

const std::vector<tensor2>& DamageMechanism::dPhi_dsigma(
    const arma::vec& sigma) const {
    // Φ = Y - Y_max with Y = 0.5 σ : M : σ → dΦ/dσ = M · σ (strain-typed).
    // M_cached_ is populated by compute_constraints (must be called first).
    dPhi_dsigma_cache_[0] = M_cached_valid_
        ? strain(arma::vec(M_cached_ * sigma))
        : tensor2::zeros(Tensor2Type::strain);
    return dPhi_dsigma_cache_;
}

const std::vector<tensor2>& DamageMechanism::kappa(
    const arma::vec& sigma, double /*DT*/, const arma::mat& /*L_ref*/) const {
    // κ^damage = ∂M/∂D · σ. For the (1-D)·L_0 model, M(D) = (1/(1-D)) · M_0,
    // so ∂M/∂D = (1/(1-D)²) · M_0. This enables weak coupling into other
    // mechanisms' Jacobian rows (they see stress perturbations from damage
    // evolution within a single FB iteration instead of only through the
    // outer-iteration stress update).
    // NB: strain-typed (an M·σ product) — see the header note on the
    // deliberate convention mix in the orchestrator's B assembly.
    if (M_cached_valid_) {
        const double D = ivc_.get("D").scalar();
        const double factor = 1.0 / ((1.0 - D) * (1.0 - D));
        kappa_cache_[0] = strain(arma::vec(factor * (M_cached_ * sigma)));
    } else {
        kappa_cache_[0] = tensor2::zeros(Tensor2Type::strain);
    }
    return kappa_cache_;
}

double DamageMechanism::stiffness_reduction() const {
    return 1.0 - ivc_.get("D").scalar();
}

arma::vec DamageMechanism::inelastic_strain() const {
    // Damage doesn't contribute a separate inelastic strain
    // It affects the stiffness instead
    return arma::zeros(6);
}

void DamageMechanism::update(
    const arma::vec& /*ds*/,
    int /*offset*/
) {
    // Explicit integration: the FB multiplier increment is unused here.
    // Advance the history and recompute damage directly from Y_max.
    double& Y_max = ivc_.get("Y_max").scalar();
    if (Y_current_ > Y_max) {
        Y_max = Y_current_;
    }
    double& D = ivc_.get("D").scalar();
    D = compute_damage(Y_current_, Y_max);
}

void DamageMechanism::tangent_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const arma::vec& Ds,
    int offset,
    arma::mat& Lt
) const {
    double D = ivc_.get("D").scalar();

    // Apply damage to tangent
    // L_damaged = (1 - D) * L
    Lt = (1.0 - D) * Lt;

    // Additional contribution from damage evolution
    // If damage is evolving, there's a coupling term
    if (dD_dY_ > simcoon::iota && D < D_c_) {
        // Compute dY/dsigma = S : sigma (using cached compliance)
        arma::vec dY_dsigma = M_cached_ * sigma;

        // Contribution: -dD/dY * sigma ⊗ dY/dsigma
        // This represents the softening due to damage evolution
        Lt -= dD_dY_ * (sigma * dY_dsigma.t());
    }
}

void DamageMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Get damage increment
    double dD = ivc_.get("D").delta_scalar();

    Wm_r = 0.0;
    Wm_ir = 0.0;
    Wm_d = 0.0;

    if (dD > simcoon::iota) {
        // Energy dissipated by damage
        // W_d = Y * dD (approximately)
        Wm_d = Y_current_ * dD;
    }
}

} // namespace simcoon
