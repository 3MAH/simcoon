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
 * @file hardening.cpp
 * @brief Implementation of hardening classes
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/hardening.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/parameter.hpp>
#include <cmath>
#include <stdexcept>

namespace simcoon {

// ============================================================================
// ISOTROPIC HARDENING IMPLEMENTATIONS
// ============================================================================

// Factory method
std::unique_ptr<IsotropicHardening> IsotropicHardening::create(IsoHardType type, int N) {
    switch (type) {
        case IsoHardType::NONE:
            return std::make_unique<NoIsotropicHardening>();
        case IsoHardType::LINEAR:
            return std::make_unique<LinearHardening>();
        case IsoHardType::POWER_LAW:
            return std::make_unique<PowerLawHardening>();
        case IsoHardType::VOCE:
            return std::make_unique<VoceHardening>();
        case IsoHardType::COMBINED_VOCE:
            return std::make_unique<CombinedVoceHardening>(N);
        default:
            return std::make_unique<NoIsotropicHardening>();
    }
}

// LinearHardening
void LinearHardening::configure(const arma::vec& props, int& offset) {
    H_ = props(offset);
    offset += 1;
}

// PowerLawHardening
void PowerLawHardening::configure(const arma::vec& props, int& offset) {
    k_ = props(offset);
    m_ = props(offset + 1);
    offset += 2;
}

double PowerLawHardening::R(double p) const {
    if (p > sim_iota) {
        return k_ * std::pow(p, m_);
    }
    return 0.0;
}

double PowerLawHardening::dR_dp(double p) const {
    if (p > sim_iota) {
        return k_ * m_ * std::pow(p, m_ - 1.0);
    }
    return 0.0;
}

// VoceHardening
void VoceHardening::configure(const arma::vec& props, int& offset) {
    Q_ = props(offset);
    b_ = props(offset + 1);
    offset += 2;
}

double VoceHardening::R(double p) const {
    return Q_ * (1.0 - std::exp(-b_ * p));
}

double VoceHardening::dR_dp(double p) const {
    return Q_ * b_ * std::exp(-b_ * p);
}

// CombinedVoceHardening
void CombinedVoceHardening::configure(const arma::vec& props, int& offset) {
    for (int i = 0; i < N_; ++i) {
        Q_(i) = props(offset + 2 * i);
        b_(i) = props(offset + 2 * i + 1);
    }
    offset += 2 * N_;
}

double CombinedVoceHardening::R(double p) const {
    double result = 0.0;
    for (int i = 0; i < N_; ++i) {
        result += Q_(i) * (1.0 - std::exp(-b_(i) * p));
    }
    return result;
}

double CombinedVoceHardening::dR_dp(double p) const {
    double result = 0.0;
    for (int i = 0; i < N_; ++i) {
        result += Q_(i) * b_(i) * std::exp(-b_(i) * p);
    }
    return result;
}

// ============================================================================
// KINEMATIC HARDENING IMPLEMENTATIONS
// ============================================================================

// Factory method
std::unique_ptr<KinematicHardening> KinematicHardening::create(KinHardType type, int N) {
    switch (type) {
        case KinHardType::NONE:
            return std::make_unique<NoKinematicHardening>();
        case KinHardType::PRAGER:
            return std::make_unique<PragerHardening>();
        case KinHardType::ARMSTRONG_FREDERICK:
            return std::make_unique<ArmstrongFrederickHardening>();
        case KinHardType::CHABOCHE:
            return std::make_unique<ChabocheHardening>(N);
        default:
            return std::make_unique<NoKinematicHardening>();
    }
}

// PragerHardening
void PragerHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    offset += 1;
}

void PragerHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec("X", arma::zeros(6), true);  // Backstress tensor
}

arma::vec PragerHardening::total_backstress(const InternalVariableCollection& ivc) const {
    return ivc.get("X").vec();
}

arma::vec PragerHardening::alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const {
    // Linear Prager: dalpha = n (no recovery)
    return n;
}

double PragerHardening::hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const {
    // H_kin = (2/3) * C
    return (2.0 / 3.0) * C_;
}

void PragerHardening::update(double dp, const arma::vec& n, InternalVariableCollection& ivc) {
    // dX = (2/3) * C * n * dp
    arma::vec& X = ivc.get("X").vec();
    X += (2.0 / 3.0) * C_ * n * dp;
}

// ArmstrongFrederickHardening
void ArmstrongFrederickHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    D_ = props(offset + 1);
    offset += 2;
}

void ArmstrongFrederickHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec("X", arma::zeros(6), true);  // Backstress tensor
}

arma::vec ArmstrongFrederickHardening::total_backstress(const InternalVariableCollection& ivc) const {
    return ivc.get("X").vec();
}

arma::vec ArmstrongFrederickHardening::alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const {
    // AF: dalpha = n - (3/2) * D * X / C
    const arma::vec& X = ivc.get("X").vec();
    if (C_ > sim_iota) {
        return n - (1.5 * D_ / C_) * X;
    }
    return n;
}

double ArmstrongFrederickHardening::hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const {
    // H_kin = (2/3) * C - D * sum(X % n)
    const arma::vec& X = ivc.get("X").vec();
    return (2.0 / 3.0) * C_ - D_ * arma::dot(X, n);
}

void ArmstrongFrederickHardening::update(double dp, const arma::vec& n, InternalVariableCollection& ivc) {
    // dX = (2/3) * C * n * dp - D * X * dp
    arma::vec& X = ivc.get("X").vec();
    X += ((2.0 / 3.0) * C_ * n - D_ * X) * dp;
}

// ChabocheHardening
void ChabocheHardening::configure(const arma::vec& props, int& offset) {
    for (int i = 0; i < N_; ++i) {
        C_(i) = props(offset + 2 * i);
        D_(i) = props(offset + 2 * i + 1);
    }
    offset += 2 * N_;
}

void ChabocheHardening::register_variables(InternalVariableCollection& ivc) {
    // Register N backstress tensors
    for (int i = 0; i < N_; ++i) {
        ivc.add_vec("X_" + std::to_string(i), arma::zeros(6), true);
    }
}

arma::vec ChabocheHardening::total_backstress(const InternalVariableCollection& ivc) const {
    arma::vec X_total = arma::zeros(6);
    for (int i = 0; i < N_; ++i) {
        X_total += ivc.get("X_" + std::to_string(i)).vec();
    }
    return X_total;
}

arma::vec ChabocheHardening::alpha_flow(int i, const arma::vec& n, const InternalVariableCollection& ivc) const {
    // AF for term i: dalpha_i = n - (3/2) * D_i * X_i / C_i
    const arma::vec& X_i = ivc.get("X_" + std::to_string(i)).vec();
    if (C_(i) > sim_iota) {
        return n - (1.5 * D_(i) / C_(i)) * X_i;
    }
    return n;
}

double ChabocheHardening::hardening_modulus(const arma::vec& n, const InternalVariableCollection& ivc) const {
    // H_kin = sum_i [ (2/3) * C_i - D_i * (X_i . n) ]
    double H_kin = 0.0;
    for (int i = 0; i < N_; ++i) {
        const arma::vec& X_i = ivc.get("X_" + std::to_string(i)).vec();
        H_kin += (2.0 / 3.0) * C_(i) - D_(i) * arma::dot(X_i, n);
    }
    return H_kin;
}

void ChabocheHardening::update(double dp, const arma::vec& n, InternalVariableCollection& ivc) {
    // dX_i = (2/3) * C_i * n * dp - D_i * X_i * dp
    for (int i = 0; i < N_; ++i) {
        arma::vec& X_i = ivc.get("X_" + std::to_string(i)).vec();
        X_i += ((2.0 / 3.0) * C_(i) * n - D_(i) * X_i) * dp;
    }
}

} // namespace simcoon
