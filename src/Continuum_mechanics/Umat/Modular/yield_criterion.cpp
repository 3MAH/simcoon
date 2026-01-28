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
    configured_ = true;
}

void YieldCriterion::configure_tresca() {
    type_ = YieldType::TRESCA;
    params_ = arma::zeros(1);  // No parameters needed
    configured_ = true;
}

void YieldCriterion::configure_drucker(double b, double n) {
    type_ = YieldType::DRUCKER;
    params_ = arma::vec({b, n});
    configured_ = true;
}

void YieldCriterion::configure_hill(double F, double G, double H, double L, double M, double N) {
    type_ = YieldType::HILL;
    params_ = arma::vec({F, G, H, L, M, N});
    configured_ = true;
}

void YieldCriterion::configure_dfa(double F, double G, double H, double L, double M, double N, double K) {
    type_ = YieldType::DFA;
    params_ = arma::vec({F, G, H, L, M, N, K});
    configured_ = true;
}

void YieldCriterion::configure_anisotropic(double P11, double P22, double P33,
                                            double P12, double P13, double P23,
                                            double P44, double P55, double P66) {
    type_ = YieldType::ANISOTROPIC;
    params_ = arma::vec({P11, P22, P33, P12, P13, P23, P44, P55, P66});
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
            throw std::runtime_error("YieldCriterion: unknown yield type");
    }
}

// ========== Yield Function Computations ==========

double YieldCriterion::equivalent_stress(const arma::vec& sigma) const {
    if (!configured_) {
        throw std::runtime_error("YieldCriterion: not configured");
    }
    return Eq_stress(sigma, type_string(), params_);
}

double YieldCriterion::equivalent_stress(const arma::vec& sigma, const arma::vec& X) const {
    if (!configured_) {
        throw std::runtime_error("YieldCriterion: not configured");
    }
    return Eq_stress(sigma - X, type_string(), params_);
}

arma::vec YieldCriterion::flow_direction(const arma::vec& sigma) const {
    if (!configured_) {
        throw std::runtime_error("YieldCriterion: not configured");
    }
    return dEq_stress(sigma, type_string(), params_);
}

arma::vec YieldCriterion::flow_direction(const arma::vec& sigma, const arma::vec& X) const {
    if (!configured_) {
        throw std::runtime_error("YieldCriterion: not configured");
    }
    return dEq_stress(sigma - X, type_string(), params_);
}

arma::vec YieldCriterion::plastic_flow(const arma::vec& sigma) const {
    // For associated flow rules, plastic flow equals flow direction
    return flow_direction(sigma);
}

arma::vec YieldCriterion::plastic_flow(const arma::vec& sigma, const arma::vec& X) const {
    // For associated flow rules, plastic flow equals flow direction
    return flow_direction(sigma, X);
}

} // namespace simcoon
