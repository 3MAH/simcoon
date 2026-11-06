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
 * @file yield_criterion.cpp
 * @brief Implementation of YieldCriterion class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/yield_criterion.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/parameter.hpp>
#include <cmath>
#include <stdexcept>

namespace simcoon {

// ========== Constructor ==========

YieldCriterion::YieldCriterion()
    : type_(YieldType::VON_MISES)
    , params_()
    , configured_(false)
{
}

// ========== Private Helpers ==========

std::string YieldCriterion::type_string() const {
    switch (type_) {
        case YieldType::VON_MISES:
            return "Mises";
        case YieldType::TRESCA:
            return "Tresca";
        case YieldType::DRUCKER:
            return "Drucker";
        case YieldType::HILL:
            return "Hill";
        case YieldType::DFA:
            return "DFA";
        case YieldType::ANISOTROPIC:
            return "Ani";
        default:
            return "Mises";
    }
}

// ========== Configuration ==========

void YieldCriterion::configure_von_mises() {
    type_ = YieldType::VON_MISES;
    params_ = arma::zeros(1);  // No parameters needed
    P_.reset();
    configured_ = true;
}

void YieldCriterion::configure_tresca() {
    type_ = YieldType::TRESCA;
    params_ = arma::zeros(1);  // No parameters needed
    P_.reset();
    configured_ = true;
}

void YieldCriterion::configure_drucker(double b, double n) {
    type_ = YieldType::DRUCKER;
    params_ = arma::vec({b, n});
    P_.reset();
    configured_ = true;
}

void YieldCriterion::configure_hill(double F, double G, double H, double L, double M, double N) {
    type_ = YieldType::HILL;
    params_ = arma::vec({F, G, H, L, M, N});
    P_ = P_Hill(params_);
    configured_ = true;
}

void YieldCriterion::configure_dfa(double F, double G, double H, double L, double M, double N, double K) {
    type_ = YieldType::DFA;
    params_ = arma::vec({F, G, H, L, M, N, K});
    P_ = P_DFA(params_);
    configured_ = true;
}

void YieldCriterion::configure_anisotropic(double P11, double P22, double P33,
                                            double P12, double P13, double P23,
                                            double P44, double P55, double P66) {
    type_ = YieldType::ANISOTROPIC;
    params_ = arma::vec({P11, P22, P33, P12, P13, P23, P44, P55, P66});
    P_ = P_Ani(params_);
    configured_ = true;
}

void YieldCriterion::configure(YieldType type, const arma::vec& props, int& offset) {
    switch (type) {
        case YieldType::VON_MISES:
            configure_von_mises();
            // No additional props consumed
            break;
        case YieldType::TRESCA:
            configure_tresca();
            // No additional props consumed
            break;
        case YieldType::DRUCKER: {
            double b = props(offset);
            double n = props(offset + 1);
            configure_drucker(b, n);
            offset += 2;
            break;
        }
        case YieldType::HILL: {
            double F = props(offset);
            double G = props(offset + 1);
            double H = props(offset + 2);
            double L = props(offset + 3);
            double M = props(offset + 4);
            double N = props(offset + 5);
            configure_hill(F, G, H, L, M, N);
            offset += 6;
            break;
        }
        case YieldType::DFA: {
            double F = props(offset);
            double G = props(offset + 1);
            double H = props(offset + 2);
            double L = props(offset + 3);
            double M = props(offset + 4);
            double N = props(offset + 5);
            double K = props(offset + 6);
            configure_dfa(F, G, H, L, M, N, K);
            offset += 7;
            break;
        }
        case YieldType::ANISOTROPIC: {
            double P11 = props(offset);
            double P22 = props(offset + 1);
            double P33 = props(offset + 2);
            double P12 = props(offset + 3);
            double P13 = props(offset + 4);
            double P23 = props(offset + 5);
            double P44 = props(offset + 6);
            double P55 = props(offset + 7);
            double P66 = props(offset + 8);
            configure_anisotropic(P11, P22, P33, P12, P13, P23, P44, P55, P66);
            offset += 9;
            break;
        }
        default:
            throw std::invalid_argument("YieldCriterion: unknown yield type");
    }
}

// ========== Yield Function Computations ==========

// the criteria.cpp kernels called by type, the quadratic ones on P_ built at configuration

double YieldCriterion::equivalent_stress(const arma::vec& sigma) const {
    if (!configured_) {
        throw std::invalid_argument("YieldCriterion: not configured");
    }
    switch (type_) {
        case YieldType::TRESCA:
            return Tresca_stress(sigma);
        case YieldType::DRUCKER:
            return Drucker_stress(sigma, params_(0), params_(1));
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC:
            return Eq_stress_P(sigma, P_);
        default:   // VON_MISES (type_string's fallback too)
            return Mises_stress(sigma);
    }
}

arma::vec YieldCriterion::flow_direction(const arma::vec& sigma) const {
    if (!configured_) {
        throw std::invalid_argument("YieldCriterion: not configured");
    }
    switch (type_) {
        case YieldType::TRESCA:
            return dTresca_stress(sigma);
        case YieldType::DRUCKER:
            return dDrucker_stress(sigma, params_(0), params_(1));
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC:
            return dEq_stress_P(sigma, P_);
        default:   // VON_MISES
            return eta_stress(sigma);
    }
}

void YieldCriterion::equivalent_stress_and_flow(const arma::vec& sigma, double& sigma_eq,
                                                arma::vec& flow) const {
    switch (type_) {
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC: {
            if (!configured_) {
                throw std::runtime_error("YieldCriterion: not configured");
            }
            // Eq_stress_P and dEq_stress_P on one P sigma: the latter is P sigma over the former
            if (arma::norm(sigma, 2) > simcoon::iota) {
                const arma::vec P_sigma = P_ * sigma;
                sigma_eq = std::pow(arma::sum(sigma % P_sigma), 0.5);
                flow = P_sigma / sigma_eq;
            } else {
                sigma_eq = 0.;
                flow = arma::zeros(6);
            }
            return;
        }
        default:
            sigma_eq = equivalent_stress(sigma);
            flow = flow_direction(sigma);
    }
}

void YieldCriterion::equivalent_stress_flow_hessian(const arma::vec& sigma, double& sigma_eq,
                                                    arma::vec& flow, arma::mat& hessian) const {
    switch (type_) {
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC: {
            if (!configured_) {
                throw std::runtime_error("YieldCriterion: not configured");
            }
            // Eq_stress_P, dEq_stress_P and ddEq_stress_P(sigma, P_) (criteria.cpp) on one P sigma
            if (arma::norm(sigma, 2) > simcoon::iota) {
                const arma::vec P_sigma = P_ * sigma;
                sigma_eq = std::pow(arma::sum(sigma % P_sigma), 0.5);
                flow = P_sigma / sigma_eq;
                if (sigma_eq < simcoon::iota) {
                    hessian = arma::zeros(6, 6);
                } else {
                    hessian = (P_ - flow * flow.t()) / sigma_eq;   // flow is ddEq_stress_P's Lambda
                }
            } else {
                sigma_eq = 0.;
                flow = arma::zeros(6);
                hessian = arma::zeros(6, 6);
            }
            return;
        }
        default:
            sigma_eq = equivalent_stress(sigma);
            flow = flow_direction(sigma);
            hessian = flow_hessian(sigma);
    }
}

arma::vec YieldCriterion::plastic_flow(const arma::vec& sigma) const {
    return flow_direction(sigma);
}

bool YieldCriterion::has_flow_hessian() const noexcept {
    switch (type_) {
        case YieldType::VON_MISES:
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC:
            return true;
        default:
            return false;
    }
}

bool YieldCriterion::has_stable_flow_direction() const noexcept {
    // Exhaustive on purpose (no default): a criterion added to YieldType
    // must make this choice explicitly — the compiler flags the omission.
    switch (type_) {
        case YieldType::VON_MISES:
        case YieldType::HILL:
        case YieldType::ANISOTROPIC:
            return true;
        case YieldType::TRESCA:   // vertex flow
        case YieldType::DRUCKER:  // pressure-sensitive
        case YieldType::DFA:      // pressure-sensitive
            return false;
    }
    return false;  // unreachable for a valid YieldType
}

arma::mat YieldCriterion::flow_hessian(const arma::vec& sigma) const {
    if (!configured_) {
        throw std::invalid_argument("YieldCriterion: not configured");
    }
    switch (type_) {
        case YieldType::VON_MISES:
            return deta_stress(sigma);
        case YieldType::HILL:
        case YieldType::DFA:
        case YieldType::ANISOTROPIC:
            return ddEq_stress_P(sigma, P_);   // ddHill_stress / ddDFA_stress / ddAni_stress
        default:
            throw std::invalid_argument(
                "YieldCriterion: no analytic flow Hessian for " + type_string());
    }
}

double YieldCriterion::equivalent_stress(const tensor2& sigma) const {
    if (!configured_) {
        throw std::invalid_argument("YieldCriterion: not configured");
    }
    if (type_ == YieldType::VON_MISES) {
        return Mises(sigma);
    }
    return equivalent_stress(sigma.to_arma_voigt());
}

tensor2 YieldCriterion::flow_direction(const tensor2& sigma) const {
    if (!configured_) {
        throw std::invalid_argument("YieldCriterion: not configured");
    }
    if (type_ == YieldType::VON_MISES) {
        return flow_normal(sigma);
    }
    return strain(flow_direction(sigma.to_arma_voigt()));
}

} // namespace simcoon
