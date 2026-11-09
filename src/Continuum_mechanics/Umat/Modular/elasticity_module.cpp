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

#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/elasticity_module.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
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
    throw std::invalid_argument("ElasticityModule: unknown isotropic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..5)");
}

static const char* conv_string(CubicConv conv) {
    switch (conv) {
        case CubicConv::EnuG: return "EnuG";
        case CubicConv::Cii:  return "Cii";
    }
    throw std::invalid_argument("ElasticityModule: unknown cubic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..1)");
}

static const char* conv_string(OrthoConv conv) {
    switch (conv) {
        case OrthoConv::EnuG: return "EnuG";
        case OrthoConv::Cii:  return "Cii";
    }
    throw std::invalid_argument("ElasticityModule: unknown orthotropic convention code "
                             + std::to_string(static_cast<int>(conv)) + " (valid: 0..1)");
}

// ========== Constructor ==========

ElasticityModule::ElasticityModule()
    : type_(ElasticityType::ISOTROPIC)
    , L_(arma::zeros(6, 6))
    , M_(arma::zeros(6, 6))
    , alpha_(arma::zeros(6))
    , hyper_potential_(HyperPotential::NEOHC)
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

    alpha_ = alpha_scalar * Ith();

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure_cubic(double C1, double C2, double C3, double alpha_scalar,
                                       CubicConv conv) {
    type_ = ElasticityType::CUBIC;

    const char* str = conv_string(conv);
    L_ = L_cubic(C1, C2, C3, str);

    M_ = M_cubic(C1, C2, C3, str);

    alpha_ = alpha_scalar * Ith();

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
        throw std::invalid_argument(
            "ElasticityModule: unknown transverse-isotropic convention code "
            + std::to_string(static_cast<int>(conv)) + " (valid: 0)");
    }

    // validate first, with the module's own message
    if (axis < 1 || axis > 3) {
        throw std::invalid_argument(
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

void ElasticityModule::configure_hyper_invariants(HyperPotential potential, const arma::vec& params,
                                                  double alpha_scalar) {
    type_ = ElasticityType::HYPER_INVARIANTS;
    hyper_potential_ = potential;
    hyper_props_ = params;

    alpha_ = alpha_scalar * Ith();

    // closed-form ground state, isotropic and exact even for an anisotropic potential (see header)
    const arma::mat I3 = arma::eye(3, 3);
    an_ = hyper_potential_anisotropy(potential, params);   // parsed ONCE, read by evaluate()
    const hyper_invariants_dW dW0 = hyper_potential_derivatives(
        potential, params, isochoric_invariants(I3, 1.0), 1.0,
        structure_tensors_push_forward(I3, an_.a0, an_.kappa_d, 1.0));
    const double K0 = dW0.dU2dJ2;
    const double mu0 = 2. * (dW0.dWdI_1_bar + dW0.dWdI_2_bar);
    L_ = L_iso(K0, mu0, "Kmu");
    M_ = M_iso(K0, mu0, "Kmu");

    refresh_tensors();
    configured_ = true;
}

void ElasticityModule::configure(ElasticityType type, const arma::vec& props, int& offset) {
    // the conv_string helpers validate the convention slot inside configure_*
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
        case ElasticityType::HYPER_INVARIANTS: {
            // props: potential, n_params, params..., alpha
            const auto potential = static_cast<HyperPotential>(static_cast<int>(props(offset)));
            const int n_params = static_cast<int>(props(offset + 1));
            if (n_params < 1) {
                throw std::invalid_argument(
                    "ElasticityModule: HYPER_INVARIANTS needs at least one potential parameter");
            }
            const arma::vec params = props.subvec(offset + 2, offset + 1 + n_params);
            const double alpha = props(offset + 2 + n_params);
            configure_hyper_invariants(potential, params, alpha);
            offset += props_count(type, n_params);
            break;
        }
        default:
            throw std::invalid_argument("ElasticityModule: unknown elasticity type");
    }
}

// ========== Derived Quantities ==========

void ElasticityModule::evaluate(const arma::vec& eps_el, int ndi, int nshr,
                                arma::vec& sigma, arma::mat& Lt) const {
    (void)nshr;
    if (!configured_) {
        throw std::invalid_argument("ElasticityModule: not configured");
    }
    if (type_ != ElasticityType::HYPER_INVARIANTS) {
        if (L_.n_rows == 6 && L_.n_cols == 6 && eps_el.n_elem == 6 ) {
            mul6(L_, eps_el, sigma);
        } else {
            sigma = L_ * eps_el;
        }
        Lt = L_;
        return;
    }
    if (ndi != 3) {   // condensed states are solved around the whole law (condensed_state.hpp)
        throw std::invalid_argument(
            "ElasticityModule: a hyperelastic block supports ndi = 3 only (got "
            + std::to_string(ndi) + ")");
    }
    // V_el = exp(eps_el) passed where F is expected: the returned box tangent is d(tau)/d(eps_el), R = I
    const arma::mat eps_t = v2t_strain(eps_el);
    arma::vec::fixed<3> e;
    arma::mat::fixed<3,3> N;
    if (!eig_sym_3x3(e, N, eps_t)) {
        throw simcoon::exception_expmat_sym("Error in expmat_sym function inside ElasticityModule::evaluate.");
    }
    const arma::vec::fixed<3> exp_e = arma::exp(e);
    const arma::mat V_el = spectral_compose_3x3(N, exp_e);
    const arma::mat b_el = spectral_compose_3x3(N, exp_e % exp_e);
    const double J_el = std::exp(arma::trace(eps_t));

    // An anisotropic potential carries reference fibre directions: push them forward with
    // V_el, the same tensor that plays F here. an_ was parsed once at configure time.
    const std::vector<arma::mat> A =
        structure_tensors_push_forward(V_el, an_.a0, an_.kappa_d, J_el);

    const hyper_invariants_dW dW =
        hyper_potential_derivatives(hyper_potential_, hyper_props_,
                                    isochoric_invariants(b_el, J_el), J_el, A);
    // tau = dpsi/deps_el per reference volume (isotropic psi); corate 3, never the solver's.
    const bool total = (kinematic_form_ == KinematicForm::total);
    if (!total && A.empty()) {
        // isotropic: closed-form box from the Lie tangent (one decomposition); Truesdell box = Lie tangent
        arma::mat Lt_Lie;
        hyper_invariants_response(dW, b_el, J_el, V_el, 4, sigma, Lt_Lie, A);
        Lt = Dtau_LieDD_2_DtauDe_stretch(Lt_Lie, e, N, v2t_stress(sigma));
        return;
    }
    hyper_invariants_response(dW, b_el, J_el, V_el, total ? corate_type_total : 3, sigma, Lt, A);
    if (total) {
        // dual stress of H_el: equal to tau when coaxial, and what makes a fibre potential exact
        sigma = t2v_stress(Kirchoff2LogStressU(v2t_stress(sigma), V_el));
    }
}

arma::vec ElasticityModule::thermal_strain(double DT) const {
    if (!configured_) {
        throw std::invalid_argument("ElasticityModule: not configured");
    }
    return alpha_ * DT;
}

} // namespace simcoon
