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

void DamageMechanism::register_variables(InternalVariableCollection& ivc) {
    // Register damage variable
    ivc.add_scalar("D", 0.0, false);

    // Register maximum damage driving force (history variable)
    ivc.add_scalar("Y_max", 0.0, false);
}

// ========== Constitutive Computations ==========

double DamageMechanism::compute_driving_force(const arma::vec& sigma, const arma::mat& S) const {
    // Damage driving force: Y = (1/2) * sigma : S : sigma
    // This is the strain energy release rate
    return 0.5 * arma::dot(sigma, S * sigma);
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
            // Linear damage evolution
            D = (Y_eff - Y_0_) / (Y_c_ - Y_0_);
            break;

        case DamageType::EXPONENTIAL:
            // Exponential damage evolution
            D = 1.0 - std::exp(-A_ * (Y_eff - Y_0_) / (Y_c_ - Y_0_));
            break;

        case DamageType::POWER_LAW:
            // Power law damage evolution
            D = std::pow((Y_eff - Y_0_) / (Y_c_ - Y_0_), n_);
            break;

        case DamageType::WEIBULL:
            // Weibull distribution-based damage
            D = 1.0 - std::exp(-std::pow((Y_eff - Y_0_) / A_, n_));
            break;
    }

    // Limit damage to critical value
    return std::min(D, D_c_);
}

double DamageMechanism::get_damage(const InternalVariableCollection& ivc) const {
    return ivc.get("D").scalar();
}

arma::mat DamageMechanism::damaged_stiffness(const arma::mat& L, double D) {
    return (1.0 - D) * L;
}

void DamageMechanism::compute_constraints(
    const arma::vec& sigma,
    const arma::mat& L,
    double DTime,
    const InternalVariableCollection& ivc,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    Phi.set_size(1);
    Y_crit.set_size(1);

    // Get current damage and history
    D_current_ = ivc.get("D").scalar();
    double Y_max = ivc.get("Y_max").scalar();

    // Compute and cache compliance (only on first use or if L changes)
    if (!M_cached_valid_) {
        M_cached_ = arma::inv(L);
        M_cached_valid_ = true;
    }

    // Compute current driving force
    Y_current_ = compute_driving_force(sigma, M_cached_);

    // Update Y_max if necessary
    double Y_eff = std::max(Y_current_, Y_max);

    // Compute damage based on driving force
    double D_new = compute_damage(Y_current_, Y_max);

    // Constraint: D - D_computed = 0 when loading, D - D_old = 0 when unloading
    // Use Fischer-Burmeister formulation: Phi <= 0 and dD >= 0
    //
    // For damage, the constraint is:
    // Phi = Y - Y_max (loading condition)
    // If Phi > 0, damage evolves; otherwise, elastic unloading

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

void DamageMechanism::compute_flow_directions(
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    std::map<std::string, arma::vec>& Lambda_map
) const {
    // Damage doesn't have a strain-like flow direction
    // The "direction" is the derivative of the driving force w.r.t. stress
    // dY/dsigma = S : sigma (for quadratic Y)

    // Not applicable for scalar damage - leave empty
}

void DamageMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const InternalVariableCollection& ivc,
    arma::mat& B,
    int row_offset
) const {
    // The Jacobian contribution for damage
    // dPhi/dD = 0 (Phi is independent of D directly)
    // The coupling comes through the stress

    // For the loading condition Phi = Y - Y_max:
    // B = dPhi/dDs = dY/dsigma * dsigma/dDs

    // Simplified: assume a unit stiffness for the constraint
    B(row_offset, row_offset) = 1.0;
}

arma::vec DamageMechanism::inelastic_strain(const InternalVariableCollection& ivc) const {
    // Damage doesn't contribute a separate inelastic strain
    // It affects the stiffness instead
    return arma::zeros(6);
}

void DamageMechanism::update(
    const arma::vec& ds,
    int offset,
    InternalVariableCollection& ivc
) {
    double dY = ds(offset);

    // Update Y_max if we're loading
    double& Y_max = ivc.get("Y_max").scalar();
    if (Y_current_ > Y_max) {
        Y_max = Y_current_;
    }

    // Compute and update damage
    double& D = ivc.get("D").scalar();
    D = compute_damage(Y_current_, Y_max);
}

void DamageMechanism::tangent_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const arma::vec& Ds,
    int offset,
    const InternalVariableCollection& ivc,
    arma::mat& Lt
) const {
    double D = ivc.get("D").scalar();

    // Apply damage to tangent
    // L_damaged = (1 - D) * L
    Lt = (1.0 - D) * Lt;

    // Additional contribution from damage evolution
    // If damage is evolving, there's a coupling term
    if (dD_dY_ > sim_iota && D < D_c_) {
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
    const InternalVariableCollection& ivc,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Get damage increment
    double dD = ivc.get("D").delta_scalar();

    Wm_r = 0.0;
    Wm_ir = 0.0;
    Wm_d = 0.0;

    if (dD > sim_iota) {
        // Energy dissipated by damage
        // W_d = Y * dD (approximately)
        Wm_d = Y_current_ * dD;
    }
}

} // namespace simcoon
