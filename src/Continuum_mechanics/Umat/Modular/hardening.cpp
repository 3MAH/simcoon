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
//
// X is a STRESS-like tensor: it shifts the stress in the yield function
// (σ - X) and pairs with the strain-like flow direction n in X:n. The
// back-strain α is stored strain-like (factor-2 shear); the product carries
// that convention, so the result is re-tagged as stress (via the 3x3, which
// is convention-free) — otherwise σ - X and X:n double-count the shear terms.
inline tensor2 backstress_t(double C, const tensor2& alpha) {
    return stress(arma::mat::fixed<3,3>(((2.0 / 3.0) * C * alpha).mat()));
}
}  // namespace

// ============================================================================
// ISOTROPIC HARDENING IMPLEMENTATIONS
// ============================================================================

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

// PowerLawHardening — onset regularization for m < 1 (rationale and the SMA
// lagrange_pow analogy are in the header Doxygen). Implementation notes only:
// the p <= 0 branch extends R linearly (slope a_reg_) so a transient negative
// FB iterate meets a finite, consistent slope instead of pow(negative, m).
namespace {
constexpr double p_reg_powerlaw = 1.0e-6;  ///< onset regularization cutoff (see header)
}

void PowerLawHardening::configure(const arma::vec& props, int& offset) {
    k_ = props(offset);
    m_ = props(offset + 1);
    offset += 2;
    // Precompute the onset-blend coefficients (value AND slope share them, so
    // the C1 match between R and dR_dp cannot drift).
    if (m_ < 1.0) {
        a_reg_ = k_ * std::pow(p_reg_powerlaw, m_ - 1.0) * (2.0 - m_);
        b_reg_ = k_ * std::pow(p_reg_powerlaw, m_ - 2.0) * (m_ - 1.0);
    }
}

double PowerLawHardening::R(double p) const {
    if (m_ >= 1.0) {
        return (p > 0.0) ? k_ * std::pow(p, m_) : 0.0;
    }
    if (p >= p_reg_powerlaw) {
        return k_ * std::pow(p, m_);
    }
    if (p <= 0.0) {
        return a_reg_ * p;
    }
    return a_reg_ * p + b_reg_ * p * p;
}

double PowerLawHardening::dR_dp(double p) const {
    if (m_ >= 1.0) {
        return (p > 0.0) ? k_ * m_ * std::pow(p, m_ - 1.0) : 0.0;
    }
    if (p >= p_reg_powerlaw) {
        return k_ * m_ * std::pow(p, m_ - 1.0);
    }
    if (p <= 0.0) {
        return a_reg_;
    }
    return a_reg_ + 2.0 * b_reg_ * p;
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

void PragerHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    offset += 1;
}

void PragerHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
}

void PragerHardening::compute_backstresses(const InternalVariableCollection& ivc,
                                           std::vector<tensor2>& X_i) const {
    X_i.resize(1, tensor2(VoigtType::stress));
    X_i[0] = backstress_t(C_, ivc.get(a_key_).as_tensor2());
}

double PragerHardening::hardening_modulus(const tensor2& n,
                                          const std::vector<tensor2>& /*X_i*/) const {
    // Exact slope n : dX/dp = (2/3) C <n,n> (tensorial contraction; = C for a
    // von Mises normal, <n,n> = 3/2). The former (2/3)C underestimated it by
    // C/3, driving a systematic Newton overshoot (frozen FB loop) and a wrong
    // continuum tangent via H_eff.
    return (2.0 / 3.0) * C_ * arma::accu(n.mat() % n.mat());
}

void PragerHardening::update(double dp, const tensor2& n, InternalVariableCollection& ivc) {
    auto& a_var = ivc.get(a_key_);
    a_var.set_tensor2(a_var.as_tensor2() + dp * n);
}

void PragerHardening::refresh_state(double dp, const tensor2& n,
                                    InternalVariableCollection& ivc) const {
    // Linear: alpha = alpha_n + dp n (exact backward Euler).
    auto& a_var = ivc.get(a_key_);
    a_var.set_tensor2(a_var.as_tensor2_start() + dp * n);
}

// ArmstrongFrederickHardening
void ArmstrongFrederickHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    D_ = props(offset + 1);
    offset += 2;
}

void ArmstrongFrederickHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
}

void ArmstrongFrederickHardening::compute_backstresses(
    const InternalVariableCollection& ivc, std::vector<tensor2>& X_i) const {
    X_i.resize(1, tensor2(VoigtType::stress));
    X_i[0] = backstress_t(C_, ivc.get(a_key_).as_tensor2());
}

double ArmstrongFrederickHardening::hardening_modulus(
    const tensor2& n, const std::vector<tensor2>& X_i) const {
    // Exact slope n : dX/dp with dX/dp = (2/3)C(n − D α):
    // H_kin = (2/3) C <n,n> − D (X : n)  (tensorial <n,n> = 3/2 for von Mises,
    // giving the classical C − D(X:n)). X:n as stress-Voigt · engineering-
    // strain-Voigt is the work-conjugate pairing of the full contraction.
    const double nn = arma::accu(n.mat() % n.mat());
    return (2.0 / 3.0) * C_ * nn - D_ * arma::dot(X_i[0].to_arma_voigt(), n.to_arma_voigt());
}

void ArmstrongFrederickHardening::update(double dp, const tensor2& n,
                                          InternalVariableCollection& ivc) {
    // dα = (n − D α) dp  (AF in α-form; equivalent to dX/dp = (2/3)C n − D X)
    auto& a_var = ivc.get(a_key_);
    const tensor2 a_t = a_var.as_tensor2();
    a_var.set_tensor2(a_t + dp * (n - D_ * a_t));
}

void ArmstrongFrederickHardening::refresh_state(double dp, const tensor2& n,
                                                InternalVariableCollection& ivc) const {
    // Backward Euler closed form: alpha = (alpha_n + dp n) / (1 + D dp).
    auto& a_var = ivc.get(a_key_);
    a_var.set_tensor2((1.0 / (1.0 + D_ * dp)) * (a_var.as_tensor2_start() + dp * n));
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
        a_keys_[i] = "a_" + std::to_string(i);
        ivc.add_vec(a_keys_[i], arma::zeros(6), true);
    }
}

void ChabocheHardening::compute_backstresses(const InternalVariableCollection& ivc,
                                             std::vector<tensor2>& X_i) const {
    // X_i = (2/3) C_i α_i, STRESS-typed (backstress_t re-tags through the 3x3;
    // a strain-typed carrier would put the factor-2 shear back — wrong under
    // any loading with shear).
    X_i.resize(static_cast<size_t>(N_), tensor2(VoigtType::stress));
    for (int i = 0; i < N_; ++i) {
        X_i[static_cast<size_t>(i)] = backstress_t(C_(i), ivc.get(a_keys_[i]).as_tensor2());
    }
}

double ChabocheHardening::hardening_modulus(const tensor2& n,
                                            const std::vector<tensor2>& X_i) const {
    // H_kin = Σ_i [ (2/3) C_i <n,n> − D_i (X_i : n) ] — exact slope per
    // branch (see ArmstrongFrederick).
    const arma::vec n_v = n.to_arma_voigt();
    const double nn = arma::accu(n.mat() % n.mat());
    double H_kin = 0.0;
    for (int i = 0; i < N_; ++i) {
        H_kin += (2.0 / 3.0) * C_(i) * nn
                 - D_(i) * arma::dot(X_i[static_cast<size_t>(i)].to_arma_voigt(), n_v);
    }
    return H_kin;
}

void ChabocheHardening::update(double dp, const tensor2& n,
                                 InternalVariableCollection& ivc) {
    // dα_i = (n − D_i α_i) dp
    for (int i = 0; i < N_; ++i) {
        auto& a_var = ivc.get(a_keys_[i]);
        const tensor2 a_t = a_var.as_tensor2();
        a_var.set_tensor2(a_t + dp * (n - D_(i) * a_t));
    }
}

void ChabocheHardening::refresh_state(double dp, const tensor2& n,
                                      InternalVariableCollection& ivc) const {
    // Per-branch backward Euler closed form (see ArmstrongFrederick).
    for (int i = 0; i < N_; ++i) {
        auto& a_var = ivc.get(a_keys_[i]);
        a_var.set_tensor2((1.0 / (1.0 + D_(i) * dp))
                          * (a_var.as_tensor2_start() + dp * n));
    }
}

} // namespace simcoon
