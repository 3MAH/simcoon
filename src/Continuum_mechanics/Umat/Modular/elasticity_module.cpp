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

void ElasticityModule::configure_isotropic(double E, double nu, double alpha_scalar) {
    type_ = ElasticityType::ISOTROPIC;

    // Compute stiffness using existing function
    L_ = L_iso(E, nu, "Enu");

    // Compute compliance
    M_ = M_iso(E, nu, "Enu");

    // Isotropic CTE
    alpha_ = arma::zeros(6);
    alpha_(0) = alpha_scalar;
    alpha_(1) = alpha_scalar;
    alpha_(2) = alpha_scalar;
    // Shear components are zero for thermal expansion

    configured_ = true;
}

void ElasticityModule::configure_cubic(double E, double nu, double G, double alpha_scalar) {
    type_ = ElasticityType::CUBIC;

    // Compute stiffness using existing function (EnuG convention)
    L_ = L_cubic(E, nu, G, "EnuG");

    // Compute compliance
    M_ = M_cubic(E, nu, G, "EnuG");

    // Cubic CTE is isotropic
    alpha_ = arma::zeros(6);
    alpha_(0) = alpha_scalar;
    alpha_(1) = alpha_scalar;
    alpha_(2) = alpha_scalar;

    configured_ = true;
}

void ElasticityModule::configure_cubic_Cii(double C11, double C12, double C44, double alpha_scalar) {
    type_ = ElasticityType::CUBIC;

    // Compute stiffness using existing function (Cii convention)
    L_ = L_cubic(C11, C12, C44, "Cii");

    // Compute compliance
    M_ = M_cubic(C11, C12, C44, "Cii");

    // Cubic CTE is isotropic
    alpha_ = arma::zeros(6);
    alpha_(0) = alpha_scalar;
    alpha_(1) = alpha_scalar;
    alpha_(2) = alpha_scalar;

    configured_ = true;
}

void ElasticityModule::configure_transverse_isotropic(double EL, double ET, double nuTL, double nuTT,
                                                      double GLT, double alpha_L, double alpha_T, int axis) {
    type_ = ElasticityType::TRANSVERSE_ISOTROPIC;

    // Compute stiffness using existing function
    L_ = L_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

    // Compute compliance
    M_ = M_isotrans(EL, ET, nuTL, nuTT, GLT, axis);

    // Transversely isotropic CTE
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
        case 3:  // z-axis is longitudinal (default)
        default:
            alpha_(0) = alpha_T;
            alpha_(1) = alpha_T;
            alpha_(2) = alpha_L;
            break;
    }

    configured_ = true;
}

void ElasticityModule::configure_orthotropic(double E1, double E2, double E3,
                                              double nu12, double nu13, double nu23,
                                              double G12, double G13, double G23,
                                              double alpha1, double alpha2, double alpha3) {
    type_ = ElasticityType::ORTHOTROPIC;

    // Compute stiffness using existing function
    // Note: L_ortho expects EnuG convention by default
    L_ = L_ortho(E1, nu12, nu13, E2, nu23, E3, G12, G13, G23, "EnuG");

    // Compute compliance
    M_ = M_ortho(E1, nu12, nu13, E2, nu23, E3, G12, G13, G23, "EnuG");

    // Orthotropic CTE
    alpha_ = arma::zeros(6);
    alpha_(0) = alpha1;
    alpha_(1) = alpha2;
    alpha_(2) = alpha3;

    configured_ = true;
}

void ElasticityModule::configure(ElasticityType type, const arma::vec& props, int& offset) {
    switch (type) {
        case ElasticityType::ISOTROPIC: {
            // props: E, nu, alpha
            double E = props(offset);
            double nu = props(offset + 1);
            double alpha = props(offset + 2);
            configure_isotropic(E, nu, alpha);
            offset += 3;
            break;
        }
        case ElasticityType::CUBIC: {
            // props: E, nu, G, alpha
            double E = props(offset);
            double nu = props(offset + 1);
            double G = props(offset + 2);
            double alpha = props(offset + 3);
            configure_cubic(E, nu, G, alpha);
            offset += 4;
            break;
        }
        case ElasticityType::TRANSVERSE_ISOTROPIC: {
            // props: EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis
            double EL = props(offset);
            double ET = props(offset + 1);
            double nuTL = props(offset + 2);
            double nuTT = props(offset + 3);
            double GLT = props(offset + 4);
            double alpha_L = props(offset + 5);
            double alpha_T = props(offset + 6);
            int axis = static_cast<int>(props(offset + 7));
            configure_transverse_isotropic(EL, ET, nuTL, nuTT, GLT, alpha_L, alpha_T, axis);
            offset += 8;
            break;
        }
        case ElasticityType::ORTHOTROPIC: {
            // props: E1, E2, E3, nu12, nu13, nu23, G12, G13, G23, alpha1, alpha2, alpha3
            double E1 = props(offset);
            double E2 = props(offset + 1);
            double E3 = props(offset + 2);
            double nu12 = props(offset + 3);
            double nu13 = props(offset + 4);
            double nu23 = props(offset + 5);
            double G12 = props(offset + 6);
            double G13 = props(offset + 7);
            double G23 = props(offset + 8);
            double alpha1 = props(offset + 9);
            double alpha2 = props(offset + 10);
            double alpha3 = props(offset + 11);
            configure_orthotropic(E1, E2, E3, nu12, nu13, nu23, G12, G13, G23,
                                 alpha1, alpha2, alpha3);
            offset += 12;
            break;
        }
        default:
            throw std::runtime_error("ElasticityModule: unknown elasticity type");
    }
}

// ========== Derived Quantities ==========

arma::mat ElasticityModule::damaged_L(double d) const {
    if (!configured_) {
        throw std::runtime_error("ElasticityModule: not configured");
    }
    return (1.0 - d) * L_;
}

arma::vec ElasticityModule::thermal_strain(double DT) const {
    if (!configured_) {
        throw std::runtime_error("ElasticityModule: not configured");
    }
    return alpha_ * DT;
}

} // namespace simcoon
