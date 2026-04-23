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

namespace {
// Conjugacy: X = (2/3) C α  (chapter 6 backstress↔back-strain relation).
// File-local helper so the convention has one canonical implementation across
// Prager / Armstrong-Frederick / Chaboche.
inline tensor2 backstress_t(double C, const tensor2& alpha) {
    return (2.0 / 3.0) * C * alpha;
}
}  // namespace

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
    if (p > simcoon::iota) {
        return k_ * std::pow(p, m_);
    }
    return 0.0;
}

double PowerLawHardening::dR_dp(double p) const {
    if (p > simcoon::iota) {
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
//
// Theoretical convention (chapter 6): the back-strain α is the actual
// thermodynamic internal variable; the backstress X = (2/3) C α is the
// conjugate generalised force (it is NOT independently stored). α is stored
// strain-like (factor-2 shear); rotation, contraction, and arithmetic flow
// through the typed tensor2 API.
void PragerHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    offset += 1;
}

void PragerHardening::register_variables(InternalVariableCollection& ivc) {
    a_key_ = key("a");
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
}

arma::vec PragerHardening::total_backstress(const InternalVariableCollection& ivc) const {
    return backstress_t(C_, ivc.get(a_key_).as_tensor2()).to_arma_voigt();
}

arma::vec PragerHardening::alpha_flow(int /*i*/, const arma::vec& n,
                                       const InternalVariableCollection& /*ivc*/) const {
    return n;   // Prager: no recall
}

double PragerHardening::hardening_modulus(const arma::vec& /*n*/,
                                           const InternalVariableCollection& /*ivc*/) const {
    return (2.0 / 3.0) * C_;
}

void PragerHardening::update(double dp, const arma::vec& n, InternalVariableCollection& ivc) {
    auto& a_var = ivc.get(a_key_);
    a_var.set_tensor2(a_var.as_tensor2() + dp * strain(n));
}

// ArmstrongFrederickHardening
void ArmstrongFrederickHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    D_ = props(offset + 1);
    offset += 2;
}

void ArmstrongFrederickHardening::register_variables(InternalVariableCollection& ivc) {
    a_key_ = key("a");
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
}

arma::vec ArmstrongFrederickHardening::total_backstress(const InternalVariableCollection& ivc) const {
    return backstress_t(C_, ivc.get(a_key_).as_tensor2()).to_arma_voigt();
}

arma::vec ArmstrongFrederickHardening::alpha_flow(
    int /*i*/, const arma::vec& n, const InternalVariableCollection& ivc) const {
    // dα/dp = n − D α  (AF in α-form; equivalent to dX/dp = (2/3)C n − D X)
    return (strain(n) - D_ * ivc.get(a_key_).as_tensor2()).to_arma_voigt();
}

double ArmstrongFrederickHardening::hardening_modulus(
    const arma::vec& n, const InternalVariableCollection& ivc) const {
    // H_kin = (2/3) C  −  D · (X : n)
    const tensor2 X_t = backstress_t(C_, ivc.get(a_key_).as_tensor2());
    return (2.0 / 3.0) * C_ - D_ * arma::dot(X_t.to_arma_voigt(), n);
}

void ArmstrongFrederickHardening::update(double dp, const arma::vec& n,
                                          InternalVariableCollection& ivc) {
    // dα = (n − D α) dp
    auto& a_var = ivc.get(a_key_);
    const tensor2 a_t = a_var.as_tensor2();
    a_var.set_tensor2(a_t + dp * (strain(n) - D_ * a_t));
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
    a_keys_.resize(static_cast<size_t>(N_));
    for (int i = 0; i < N_; ++i) {
        a_keys_[i] = key("a_" + std::to_string(i));
        ivc.add_vec(a_keys_[i], arma::zeros(6), true);
    }
}

arma::vec ChabocheHardening::total_backstress(const InternalVariableCollection& ivc) const {
    // X = Σ_i (2/3) C_i α_i
    tensor2 X_t = tensor2::zeros(VoigtType::strain);
    for (int i = 0; i < N_; ++i) {
        X_t += backstress_t(C_(i), ivc.get(a_keys_[i]).as_tensor2());
    }
    return X_t.to_arma_voigt();
}

arma::vec ChabocheHardening::alpha_flow(int i, const arma::vec& n,
                                         const InternalVariableCollection& ivc) const {
    // dα_i/dp = n − D_i α_i
    return (strain(n) - D_(i) * ivc.get(a_keys_[i]).as_tensor2())
        .to_arma_voigt();
}

double ChabocheHardening::hardening_modulus(const arma::vec& n,
                                              const InternalVariableCollection& ivc) const {
    // H_kin = Σ_i [ (2/3) C_i − D_i (X_i : n) ]
    double H_kin = 0.0;
    for (int i = 0; i < N_; ++i) {
        const tensor2 X_i_t = backstress_t(C_(i), ivc.get(a_keys_[i]).as_tensor2());
        H_kin += (2.0 / 3.0) * C_(i) - D_(i) * arma::dot(X_i_t.to_arma_voigt(), n);
    }
    return H_kin;
}

void ChabocheHardening::update(double dp, const arma::vec& n,
                                 InternalVariableCollection& ivc) {
    // dα_i = (n − D_i α_i) dp
    const tensor2 n_t = strain(n);
    for (int i = 0; i < N_; ++i) {
        auto& a_var = ivc.get(a_keys_[i]);
        const tensor2 a_t = a_var.as_tensor2();
        a_var.set_tensor2(a_t + dp * (n_t - D_(i) * a_t));
    }
}

} // namespace simcoon
