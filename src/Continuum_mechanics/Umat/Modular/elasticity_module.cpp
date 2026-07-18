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
 * @file elasticity_module.cpp
 * @brief Implementation of ElasticityModule class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/elasticity_module.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <stdexcept>

namespace simcoon {

// Enum -> convention string of the classical builders (single source of the
// parameterization math — L_iso/L_cubic/L_ortho in constitutive.cpp).
static const char* conv_string(IsoConv conv) {
    switch (conv) {
        case IsoConv::Enu:      return "Enu";
        case IsoConv::nuE:      return "nuE";
        case IsoConv::Kmu:      return "Kmu";
        case IsoConv::muK:      return "muK";
        case IsoConv::lambdamu: return "lambdamu";
        case IsoConv::mulambda: return "mulambda";
    }
    throw std::runtime_error("ElasticityModule: unknown isotropic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..5)");
}

static const char* conv_string(CubicConv conv) {
    switch (conv) {
        case CubicConv::EnuG: return "EnuG";
        case CubicConv::Cii:  return "Cii";
    }
    throw std::runtime_error("ElasticityModule: unknown cubic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..1)");
}

static const char* conv_string(OrthoConv conv) {
    switch (conv) {
        case OrthoConv::EnuG: return "EnuG";
        case OrthoConv::Cii:  return "Cii";
    }
    throw std::runtime_error("ElasticityModule: unknown orthotropic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..1)");
}

// ========== Constructor ==========

ElasticityModule::ElasticityModule()
    : type_(ElasticityType::ISOTROPIC)
    , L_(arma::zeros(6, 6))
    , M_(arma::zeros(6, 6))
    , alpha_(arma::zeros(6))
    , configured_(false)
{
}

// ========== Configuration ==========

void ElasticityModule::configure_isotropic(double C1, double C2, double alpha_scalar,
                                           IsoConv conv) {
    type_ = ElasticityType::ISOTROPIC;

    const char* str = conv_string(conv);
    L_ = L_iso(C1, C2, str);

    M_ = M_iso(C1, C2, str);

    alpha_ = arma::zeros(6);
    alpha_(0) = alpha_scalar;
    alpha_(1) = alpha_scalar;
    alpha_(2) = alpha_scalar;
    // Shear components are zero for thermal expansion

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure_cubic(double C1, double C2, double C3, double alpha_scalar,
                                       CubicConv conv) {
    type_ = ElasticityType::CUBIC;

    const char* str = conv_string(conv);
    L_ = L_cubic(C1, C2, C3, str);

    M_ = M_cubic(C1, C2, C3, str);

    alpha_ = arma::zeros(6);
    alpha_(0) = alpha_scalar;
    alpha_(1) = alpha_scalar;
    alpha_(2) = alpha_scalar;

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure_transverse_isotropic(double EL, double ET, double nuTL, double nuTT,
                                                      double GLT, double alpha_L, double alpha_T, int axis,
                                                      IsotransConv conv) {
    type_ = ElasticityType::TRANSVERSE_ISOTROPIC;

    // Single parameterization today — reject unknown codes now so a future
    // convention cannot be silently misread from old props streams.
    if (conv != IsotransConv::EnuG) {
        throw std::runtime_error(
            "ElasticityModule: unknown transverse-isotropic convention code "
            + std::to_string(static_cast<int>(conv)) + " (valid: 0)");
    }

    // Validate BEFORE calling the library builders: L_isotrans terminates the
    // process (exit(0)!) on an invalid axis — a throw here keeps the failure
    // catchable and attributable.
    if (axis < 1 || axis > 3) {
        throw std::runtime_error(
            "ElasticityModule: transverse-isotropic axis must be 1, 2 or 3, got "
            + std::to_string(axis));
    }

    L_ = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

    M_ = M_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

    alpha_ = arma::zeros(6);
    switch (axis) {
        case 1:  // x-axis is longitudinal
            alpha_(0) = alpha_L;
            alpha_(1) = alpha_T;
            alpha_(2) = alpha_T;
            break;
        case 2:  // y-axis is longitudinal
            alpha_(0) = alpha_T;
            alpha_(1) = alpha_L;
            alpha_(2) = alpha_T;
            break;
        default:  // 3: z-axis is longitudinal
            alpha_(0) = alpha_T;
            alpha_(1) = alpha_T;
            alpha_(2) = alpha_L;
            break;
    }

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure_orthotropic(double C1, double C2, double C3,
                                              double C4, double C5, double C6,
                                              double C7, double C8, double C9,
                                              double alpha1, double alpha2, double alpha3,
                                              OrthoConv conv) {
    type_ = ElasticityType::ORTHOTROPIC;

    // L_ortho/M_ortho "EnuG" slot order is (E1, E2, E3, nu12, nu13, nu23,
    // G12, G13, G23) — NOT grouped (E, nu) per direction.
    const char* str = conv_string(conv);
    L_ = L_ortho(C1, C2, C3, C4, C5, C6, C7, C8, C9, str);

    M_ = M_ortho(C1, C2, C3, C4, C5, C6, C7, C8, C9, str);

    alpha_ = arma::zeros(6);
    alpha_(0) = alpha1;
    alpha_(1) = alpha2;
    alpha_(2) = alpha3;

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure(ElasticityType type, const arma::vec& props, int& offset) {
    // Every block starts with the convention slot; the conv_string helpers
    // (and the isotrans guard) validate the code inside the configure_* call.
    switch (type) {
        case ElasticityType::ISOTROPIC: {
            // props: conv, C1, C2, alpha
            IsoConv conv = static_cast<IsoConv>(static_cast<int>(props(offset)));
            double C1 = props(offset + 1);
            double C2 = props(offset + 2);
            double alpha = props(offset + 3);
            configure_isotropic(C1, C2, alpha, conv);
            offset += 4;
            break;
        }
        case ElasticityType::CUBIC: {
            // props: conv, C1, C2, C3, alpha
            CubicConv conv = static_cast<CubicConv>(static_cast<int>(props(offset)));
            double C1 = props(offset + 1);
            double C2 = props(offset + 2);
            double C3 = props(offset + 3);
            double alpha = props(offset + 4);
            configure_cubic(C1, C2, C3, alpha, conv);
            offset += 5;
            break;
        }
        case ElasticityType::TRANSVERSE_ISOTROPIC: {
            // props: conv, EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis
            IsotransConv conv = static_cast<IsotransConv>(static_cast<int>(props(offset)));
            double EL = props(offset + 1);
            double ET = props(offset + 2);
            double nuTL = props(offset + 3);
            double nuTT = props(offset + 4);
            double GLT = props(offset + 5);
            double alpha_L = props(offset + 6);
            double alpha_T = props(offset + 7);
            int axis = static_cast<int>(props(offset + 8));
            configure_transverse_isotropic(EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis, conv);
            offset += 9;
            break;
        }
        case ElasticityType::ORTHOTROPIC: {
            // props: conv, C1..C9, alpha1, alpha2, alpha3
            OrthoConv conv = static_cast<OrthoConv>(static_cast<int>(props(offset)));
            double C1 = props(offset + 1);
            double C2 = props(offset + 2);
            double C3 = props(offset + 3);
            double C4 = props(offset + 4);
            double C5 = props(offset + 5);
            double C6 = props(offset + 6);
            double C7 = props(offset + 7);
            double C8 = props(offset + 8);
            double C9 = props(offset + 9);
            double alpha1 = props(offset + 10);
            double alpha2 = props(offset + 11);
            double alpha3 = props(offset + 12);
            configure_orthotropic(C1, C2, C3, C4, C5, C6, C7, C8, C9,
                                 alpha1, alpha2, alpha3, conv);
            offset += 13;
            break;
        }
        default:
            throw std::runtime_error("ElasticityModule: unknown elasticity type");
    }
}

// ========== Derived Quantities ==========

arma::vec ElasticityModule::thermal_strain(double DT) const {
    if (!configured_) {
        throw std::runtime_error("ElasticityModule: not configured");
    }
    return alpha_ * DT;
}

} // namespace simcoon
