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

double DamageMechanism::compute_driving_force(const arma::vec& sigma_eff, const arma::mat& /*S*/) const {
    // Y = -dpsi/dD = psi_0, the UNDAMAGED energy, on the effective stress: exact for a linear
    // block (M = the tangent compliance of the undamaged block, an approximation of psi_0 for a
    // hyperelastic one).
    return 0.5 * arma::dot(sigma_eff, (M_cached_t_ * stress(sigma_eff)).to_arma_voigt());
}

double DamageMechanism::compute_damage(double Y_eff) const {
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

    // Cache the compliance (raw + typed), keyed on the stiffness it inverts: L
    // is the tangent of the current elastic state, which a state-dependent
    // (hyperelastic) block moves at every refresh, while a linear block keeps
    // it fixed and pays a single inversion.
    if (!M_cached_valid_ || !arma::approx_equal(L, L_cached_, "absdiff", 0.0)) {
        L_cached_ = L;
        M_cached_ = arma::inv(L);
        M_cached_t_ = tensor4(M_cached_, Tensor4Type::compliance);
        M_cached_valid_ = true;
    }

    // Compute current driving force
    Y_current_ = compute_driving_force(sigma, M_cached_);

    const double Y_eff = std::max(Y_current_, Y_max);

    // The row is always satisfied: D is a fixed point of update() (D = f(max(Y_max at the step
    // start, Y))), checked by consistency_residual(), not an FB multiplier.
    Phi(0) = -1.0;

    // Critical value for convergence
    Y_crit(0) = 1.0;

    // dD/dY for the tangent, only while damage grows in this increment (Y beyond the history
    // maximum at the start of the step): under unloading the response is (1-D) L, no softening
    const bool loading = Y_current_ > ivc_.get("Y_max").scalar_start();
    if (loading && Y_eff > Y_0_ && D_current_ < D_c_) {
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
    const arma::vec& /*sigma*/) const {
    // Damage is not solved by the FB rows (fixed point, consistency_residual): no stress
    // coupling in the local Jacobian. The other mechanisms see the effective stress.
    static const std::vector<tensor2> none;
    return none;
}

const std::vector<tensor2>& DamageMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/, const arma::mat& /*L_ref*/) const {
    return kappa_cache_;   // zero: no multiplier carried
}

double DamageMechanism::consistency_residual(const arma::vec& sigma_eff) const {
    if (!M_cached_valid_) {
        return 0.;
    }
    const double Y = compute_driving_force(sigma_eff, M_cached_);
    const double Y_used = ivc_.get("Y_max").scalar();
    const double target = std::max(ivc_.get("Y_max").scalar_start(), Y);
    return std::abs(target - Y_used) / std::max({std::abs(target), Y_0_, 1e-12});
}

arma::mat DamageMechanism::stress_map(const arma::vec& sigma_eff) const {
    // sigma = (1 - D) sigma_eff with dD = D'(Y) (S sigma_eff) . d sigma_eff while damage grows
    arma::mat Q = stiffness_reduction() * arma::eye(6, 6);
    if (dD_dY_ > simcoon::iota && ivc_.get("D").scalar() < D_c_ && M_cached_valid_) {
        Q -= dD_dY_ * (sigma_eff * (M_cached_ * sigma_eff).t());
    }
    return Q;
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
    // Fixed point on D (the FB multiplier increment is unused): D = f(max(Y_max at step
    // start, Y of this iterate)). The history is the one of the START of the step, never the
    // running maximum over iterates, whose overshooting first iterate would freeze too much
    // damage; the orchestrator iterates until consistency_residual() vanishes.
    double& Y_max = ivc_.get("Y_max").scalar();
    Y_max = std::max(ivc_.get("Y_max").scalar_start(), Y_current_);
    double& D = ivc_.get("D").scalar();
    D = compute_damage(Y_max);
}

void DamageMechanism::tangent_contribution(
    const arma::vec& /*sigma*/,
    const arma::mat& /*L*/,
    const arma::vec& /*Ds*/,
    int /*offset*/,
    arma::mat& /*Lt*/
) const {
    // Nothing here: sigma = (1 - D) sigma_eff, so the orchestrator left-multiplies the tangent
    // of the effective (undamaged) composition by stress_map().
}

void DamageMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    const arma::vec& E_start,
    const arma::vec& E_end,
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
        // Energy dissipated by damage: Y dD, Y = psi_0 being the force conjugate to D
        Wm_d = Y_current_ * dD;
    }
}

} // namespace simcoon
