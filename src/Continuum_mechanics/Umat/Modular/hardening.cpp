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
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/parameter.hpp>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

namespace simcoon {

namespace {
// typed-path operations in order, one rounding per statement: bit-identical unless FMA crosses statements

// back-strain registered under key at index; the name is checked in debug builds only
template <typename Collection>
auto& back_strain(Collection& ivc, size_t index, const std::string& key) {
    assert(index < ivc.size() && ivc.at(index).name() == key);
    (void)key;
    return ivc.at(index);
}

// alpha and n engineering strain Voigt (doubled shear); the entrywise forms below
// alpha <- alpha_from + dp n
inline void add_flow(const arma::vec& alpha_from, double dp, const arma::vec& n, arma::vec& alpha) {
    for (int i = 0; i < 6; ++i) {
        const double t = n[i] * dp;
        alpha[i] = alpha_from[i] + t;
    }
}

// alpha <- alpha + dp (n - D alpha)
inline void add_flow_recovery(InternalVariable& a_var, double dp, double D, const arma::vec& n) {
    arma::vec& a = a_var.raw_voigt();
    for (int i = 0; i < 6; ++i) {
        const double t1 = a[i] * D;
        const double t2 = n[i] - t1;
        const double t3 = t2 * dp;
        a[i] = a[i] + t3;
    }
}

// alpha <- s (alpha_n + dp n), the Armstrong-Frederick backward-Euler closed form
inline void refresh_from_start(InternalVariable& a_var, double dp, double s, const arma::vec& n) {
    const arma::vec& a_n = a_var.raw_voigt_start();
    arma::vec& a = a_var.raw_voigt();
    for (int i = 0; i < 6; ++i) {
        const double t = n[i] * dp;
        const double u = a_n[i] + t;
        a[i] = u * s;
    }
}

// <n, n> of a strain-like n on its engineering components (shear weight 1/2)
inline double nn_strain(const arma::vec& n) {
    return contract_voigt(n.memptr(), n.memptr(), voigt_shear_weight(Tensor2Type::strain));
}
}  // namespace

// ==== ISOTROPIC HARDENING IMPLEMENTATIONS ====

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

// p <= 0: R extended linearly (slope a_reg_), so a negative FB iterate meets a finite slope
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

void PowerLawHardening::R_and_dR(double p, double& R, double& dR) const {
    if (p > 0.0 && (m_ >= 1.0 || p >= p_reg_powerlaw)) {   // one pow: dR = m R / p
        R = k_ * std::pow(p, m_);
        dR = m_ * R / p;
        return;
    }
    R = this->R(p);
    dR = dR_dp(p);
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

void VoceHardening::R_and_dR(double p, double& R, double& dR) const {
    const double e = std::exp(-b_ * p);
    R = Q_ * (1.0 - e);
    dR = Q_ * b_ * e;
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

void CombinedVoceHardening::R_and_dR(double p, double& R, double& dR) const {
    R = 0.0;
    dR = 0.0;
    for (int i = 0; i < N_; ++i) {
        const double e = std::exp(-b_(i) * p);
        R += Q_(i) * (1.0 - e);
        dR += Q_(i) * b_(i) * e;
    }
}

// ==== KINEMATIC HARDENING IMPLEMENTATIONS ====

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

arma::vec KinematicHardening::backstress_voigt(int i, const arma::vec& a) const {
    // X = (2/3) C alpha: alpha carries the doubled shear, X does not
    const double k = (2.0 / 3.0) * branch_modulus(i);
    arma::vec::fixed<6> X;
    for (int c = 0; c < 3; ++c) X[c] = k * a[c];
    for (int c = 3; c < 6; ++c) X[c] = k * (0.5 * a[c]);
    return X;
}

void KinematicHardening::compute_backstresses(const InternalVariableCollection& ivc, arma::mat& X_i) const {
    const int n = num_backstresses();
    X_i.set_size(6, n);
    for (int i = 0; i < n; ++i) {
        X_i.col(i) = backstress_voigt(i, back_strain_variable(ivc, i)->raw_voigt());
    }
}

double KinematicHardening::stored_increment(const InternalVariableCollection& ivc) const {
    double W = 0.0;
    for (int i = 0; i < num_backstresses(); ++i) {
        const InternalVariable& a = *back_strain_variable(ivc, i);
        // X linear in alpha: the trapezoid is the back stress of the mean back-strain
        W += arma::dot(backstress_voigt(i, arma::vec(0.5*(a.raw_voigt_start() + a.raw_voigt()))), a.delta_vec());
    }
    return W;
}

void PragerHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    offset += 1;
}

void PragerHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
    a_index_ = ivc.index_of(a_key_);
}

double PragerHardening::hardening_modulus(const arma::vec& n,
                                          const arma::mat& /*X_i*/) const {
    // exact slope n : dX/dp = (2/3) C <n,n> (= C for a von Mises normal, <n,n> = 3/2)
    return (2.0 / 3.0) * C_ * nn_strain(n);
}

void PragerHardening::update(double dp, const arma::vec& n, InternalVariableCollection& ivc) {
    auto& a_var = back_strain(ivc, a_index_, a_key_);
    add_flow(a_var.raw_voigt(), dp, n, a_var.raw_voigt());   // alpha + dp n
}

void PragerHardening::refresh_state(double dp, const arma::vec& n,
                                    InternalVariableCollection& ivc) const {
    // Linear: alpha = alpha_n + dp n (exact backward Euler).
    auto& a_var = back_strain(ivc, a_index_, a_key_);
    add_flow(a_var.raw_voigt_start(), dp, n, a_var.raw_voigt());
}

void PragerHardening::backward_euler_factors(double dp, const arma::vec& n,
                                             const InternalVariableCollection& /*ivc*/,
                                             double& gamma, arma::vec& beta) const {
    gamma = dp * (2.0 / 3.0) * C_;
    beta = (2.0 / 3.0) * C_ * n;
}

// ArmstrongFrederickHardening
void ArmstrongFrederickHardening::configure(const arma::vec& props, int& offset) {
    C_ = props(offset);
    D_ = props(offset + 1);
    offset += 2;
}

void ArmstrongFrederickHardening::register_variables(InternalVariableCollection& ivc) {
    ivc.add_vec(a_key_, arma::zeros(6), true);   // back-strain (strain-like)
    a_index_ = ivc.index_of(a_key_);
}

double ArmstrongFrederickHardening::hardening_modulus(
    const arma::vec& n, const arma::mat& X_i) const {
    // H_kin = (2/3) C <n,n> - D (X:n); X:n = stress-Voigt . engineering-strain-Voigt
    const arma::vec::fixed<6> X = X_i.col(0);
    return (2.0 / 3.0) * C_ * nn_strain(n) - D_ * arma::dot(X, n);
}

void ArmstrongFrederickHardening::update(double dp, const arma::vec& n,
                                          InternalVariableCollection& ivc) {
    // dα = (n − D α) dp  (AF in α-form; equivalent to dX/dp = (2/3)C n − D X)
    add_flow_recovery(back_strain(ivc, a_index_, a_key_), dp, D_, n);
}

void ArmstrongFrederickHardening::refresh_state(double dp, const arma::vec& n,
                                                InternalVariableCollection& ivc) const {
    // Backward Euler closed form: alpha = (alpha_n + dp n) / (1 + D dp).
    refresh_from_start(back_strain(ivc, a_index_, a_key_), dp, (1.0 / (1.0 + D_ * dp)), n);
}

void ArmstrongFrederickHardening::backward_euler_factors(double dp, const arma::vec& n,
                                                         const InternalVariableCollection& ivc,
                                                         double& gamma, arma::vec& beta) const {
    const double f = (2.0 / 3.0) * C_ / (1.0 + D_ * dp);
    gamma = dp * f;
    beta = f * (n - D_ * back_strain(ivc, a_index_, a_key_).raw_voigt());
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
    a_index_.resize(static_cast<size_t>(N_));
    for (int i = 0; i < N_; ++i) {
        a_keys_[i] = "a_" + std::to_string(i);
        ivc.add_vec(a_keys_[i], arma::zeros(6), true);
        a_index_[i] = ivc.index_of(a_keys_[i]);
    }
}

double ChabocheHardening::hardening_modulus(const arma::vec& n,
                                            const arma::mat& X_i) const {
    // H_kin = Σ_i [ (2/3) C_i <n,n> − D_i (X_i : n) ] — exact slope per
    // branch (see ArmstrongFrederick).
    const double nn = nn_strain(n);
    double H_kin = 0.0;
    for (int i = 0; i < N_; ++i) {
        const arma::vec::fixed<6> X = X_i.col(i);
        H_kin += (2.0 / 3.0) * C_(i) * nn - D_(i) * arma::dot(X, n);
    }
    return H_kin;
}

void ChabocheHardening::update(double dp, const arma::vec& n,
                                 InternalVariableCollection& ivc) {
    // dα_i = (n − D_i α_i) dp
    for (int i = 0; i < N_; ++i) {
        add_flow_recovery(back_strain(ivc, a_index_[i], a_keys_[i]), dp, D_(i), n);
    }
}

void ChabocheHardening::refresh_state(double dp, const arma::vec& n,
                                      InternalVariableCollection& ivc) const {
    // Per-branch backward Euler closed form (see ArmstrongFrederick).
    for (int i = 0; i < N_; ++i) {
        refresh_from_start(back_strain(ivc, a_index_[i], a_keys_[i]), dp, (1.0 / (1.0 + D_(i) * dp)), n);
    }
}

void ChabocheHardening::backward_euler_factors(double dp, const arma::vec& n,
                                               const InternalVariableCollection& ivc,
                                               double& gamma, arma::vec& beta) const {
    gamma = 0.0;
    beta = arma::zeros(6);
    for (int i = 0; i < N_; ++i) {
        const double f = (2.0 / 3.0) * C_(i) / (1.0 + D_(i) * dp);
        gamma += dp * f;
        beta += f * (n - D_(i) * back_strain(ivc, a_index_[i], a_keys_[i]).raw_voigt());
    }
}

} // namespace simcoon
