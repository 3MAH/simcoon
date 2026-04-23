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
 * @file plasticity_mechanism.cpp
 * @brief Implementation of PlasticityMechanism class
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/parameter.hpp>
#include <stdexcept>

namespace simcoon {

// ========== Constructor ==========

PlasticityMechanism::PlasticityMechanism(
    YieldType yield_type,
    IsoHardType iso_type,
    KinHardType kin_type,
    int N_iso,
    int N_kin
)
    : StrainMechanism("plasticity")
    , yield_type_(yield_type)
    , yield_(std::make_unique<YieldCriterion>())
    , iso_hard_(IsotropicHardening::create(iso_type, N_iso))
    , kin_hard_(KinematicHardening::create(kin_type, N_kin))
    , sigma_Y_(0.0)
    , flow_dir_(arma::zeros(6))
    , kappa_(arma::zeros(6))
{
    // Configure parameter-free yield criteria immediately.
    // Others will be configured in configure() when props are available.
    switch (yield_type) {
        case YieldType::VON_MISES:
            yield_->configure_von_mises();
            break;
        case YieldType::TRESCA:
            yield_->configure_tresca();
            break;
        default:
            break;
    }
}

// ========== Configuration ==========

void PlasticityMechanism::configure(const arma::vec& props, int& offset) {
    // Props layout:
    // 1. sigma_Y (initial yield stress)
    // 2. Yield criterion parameters (if needed)
    // 3. Isotropic hardening parameters
    // 4. Kinematic hardening parameters

    sigma_Y_ = props(offset);
    offset += 1;

    // Configure yield criterion from props if not already done
    // (von Mises and Tresca are configured in the constructor since they
    // need no parameters; others read their parameters from props here)
    if (!yield_->is_configured()) {
        yield_->configure(yield_type_, props, offset);
    }

    // Configure hardening models
    iso_hard_->configure(props, offset);
    kin_hard_->configure(props, offset);
}

void PlasticityMechanism::register_variables(InternalVariableCollection& ivc) {
    p_key_  = key("p");
    EP_key_ = key("EP");
    ivc.add_scalar(p_key_, 0.0, false);
    ivc.add_vec(EP_key_, arma::zeros(6), true);
    kin_hard_->register_variables(ivc);
}

// Override to forward the prefix to our owned kinematic-hardening instance.
void PlasticityMechanism::set_ivc_prefix(const std::string& prefix) {
    StrainMechanism::set_ivc_prefix(prefix);
    kin_hard_->set_ivc_prefix(prefix);
}

// ========== Yield-side queries ==========

double PlasticityMechanism::equivalent_stress(
    const arma::vec& sigma, const InternalVariableCollection& ivc) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->equivalent_stress(sigma);
    }
    return yield_->equivalent_stress(sigma - kin_hard_->total_backstress(ivc));
}

double PlasticityMechanism::equivalent_stress(
    const tensor2& sigma, const InternalVariableCollection& ivc) const {
    return equivalent_stress(sigma.to_arma_voigt(), ivc);
}

double PlasticityMechanism::yield_function(
    const arma::vec& sigma, const InternalVariableCollection& ivc) const {
    const double p = ivc.get(p_key_).scalar();
    return equivalent_stress(sigma, ivc) - iso_hard_->R(p) - sigma_Y_;
}

double PlasticityMechanism::yield_function(
    const tensor2& sigma, const InternalVariableCollection& ivc) const {
    return yield_function(sigma.to_arma_voigt(), ivc);
}

arma::vec PlasticityMechanism::flow_direction(
    const arma::vec& sigma, const InternalVariableCollection& ivc) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->flow_direction(sigma);
    }
    return yield_->flow_direction(sigma - kin_hard_->total_backstress(ivc));
}

tensor2 PlasticityMechanism::flow_direction(
    const tensor2& sigma, const InternalVariableCollection& ivc) const {
    return strain(flow_direction(sigma.to_arma_voigt(), ivc));
}

// ========== Constitutive Computations ==========

void PlasticityMechanism::compute_constraints(
    const arma::vec& sigma,
    const arma::vec& /*E_total*/,
    const arma::mat& L,
    double /*DTime*/,
    const InternalVariableCollection& ivc,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    // Get internal variables
    double p = ivc.get(p_key_).scalar();

    // Shifted stress (sigma - X); skip the subtraction under pure isotropic
    // hardening since total_backstress would return zeros(6).
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    double sigma_eq;
    if (has_kin) {
        const arma::vec sigma_eff = sigma - kin_hard_->total_backstress(ivc);
        sigma_eq = yield_->equivalent_stress(sigma_eff);
        flow_dir_ = yield_->flow_direction(sigma_eff);
    } else {
        sigma_eq = yield_->equivalent_stress(sigma);
        flow_dir_ = yield_->flow_direction(sigma);
    }

    const double R = iso_hard_->R(p);
    const double dR_dp = iso_hard_->dR_dp(p);

    // Yield function: Phi = sigma_eq - R - sigma_Y
    Phi.set_size(1);
    Phi(0) = sigma_eq - R - sigma_Y_;

    Y_crit.set_size(1);
    Y_crit(0) = std::max(sigma_Y_, 1.0);

    kappa_ = L * flow_dir_;
    H_total_ = dR_dp + kin_hard_->hardening_modulus(flow_dir_, ivc);
}

void PlasticityMechanism::compute_flow_directions(
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    std::map<std::string, arma::vec>& Lambda_map
) const {
    const int N_kin = kin_hard_->num_backstresses();
    const arma::vec n = (N_kin == 0)
        ? yield_->flow_direction(sigma)
        : yield_->flow_direction(sigma - kin_hard_->total_backstress(ivc));

    Lambda_map["EP"] = n;

    for (int i = 0; i < N_kin; ++i) {
        Lambda_map["a_" + std::to_string(i)] = kin_hard_->alpha_flow(i, n, ivc);
    }
}

const std::vector<arma::vec>& PlasticityMechanism::dPhi_dsigma(
    const arma::vec& /*sigma*/, const InternalVariableCollection& /*ivc*/) const {
    // Requires compute_constraints to have been called this FB iteration.
    dPhi_dsigma_cache_[0] = flow_dir_;
    return dPhi_dsigma_cache_;
}

const std::vector<arma::vec>& PlasticityMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/,
    const arma::mat& /*L_ref*/, const InternalVariableCollection& /*ivc*/) const {
    // kappa = L_ref · n, already computed by compute_constraints.
    kappa_cache_[0] = kappa_;
    return kappa_cache_;
}

void PlasticityMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const InternalVariableCollection& ivc,
    arma::mat& B,
    int row_offset
) const {
    // B = -dPhi/dsigma * L * dPhi/dsigma + K
    // where K includes isotropic and kinematic hardening contributions

    // dPhi/dsigma * L * dPhi/dsigma
    double dPhi_L_dPhi = arma::dot(flow_dir_, kappa_);

    // B(0,0) = -dPhi_L_dPhi + H_total
    // Note: H_total already includes both iso and kin contributions
    B(row_offset, row_offset) = -dPhi_L_dPhi + H_total_;
}

arma::vec PlasticityMechanism::inelastic_strain(const InternalVariableCollection& ivc) const {
    return ivc.get(EP_key_).vec();
}

void PlasticityMechanism::update(
    const arma::vec& ds,
    int offset,
    InternalVariableCollection& ivc
) {
    double dp = ds(offset);

    if (dp > simcoon::iota) {
        // Update accumulated plastic strain
        double& p = ivc.get(p_key_).scalar();
        p += dp;

        // Update plastic strain
        arma::vec& EP = ivc.get(EP_key_).vec();
        EP += dp * flow_dir_;

        // Update kinematic hardening variables
        kin_hard_->update(dp, flow_dir_, ivc);
    }
}

void PlasticityMechanism::tangent_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const arma::vec& Ds,
    int offset,
    const InternalVariableCollection& ivc,
    arma::mat& Lt
) const {
    double Dp = Ds(offset);

    if (Dp > simcoon::iota) {
        // Consistent tangent for plasticity:
        // Lt = L - (L * n) ⊗ (L * n) / H_eff
        // where H_eff = n : L : n + H_total

        double n_L_n = arma::dot(flow_dir_, kappa_);
        double H_eff = n_L_n + H_total_;

        if (std::abs(H_eff) > simcoon::iota) {
            // Lt -= (kappa ⊗ kappa) / H_eff
            Lt -= (kappa_ * kappa_.t()) / H_eff;
        }
    }
}

void PlasticityMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    const InternalVariableCollection& ivc,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Get increments
    double Dp = ivc.get(p_key_).delta_scalar();
    arma::vec DEP = ivc.get(EP_key_).delta_vec();

    // Average stress during increment
    arma::vec sigma_avg = 0.5 * (sigma_start + sigma);

    // Get backstress (if any)
    arma::vec X = kin_hard_->total_backstress(ivc);

    // Dissipated work: sigma : dEP
    Wm_d = arma::dot(sigma_avg, DEP);

    // For kinematic hardening, part of the work is stored
    // Stored work in backstress: X : dEP (approximately)
    if (kin_hard_->num_backstresses() > 0) {
        Wm_ir = arma::dot(X, DEP);
        Wm_d -= Wm_ir;
    } else {
        Wm_ir = 0.0;
    }

    // Recoverable work (elastic) - handled by ModularUMAT
    Wm_r = 0.0;
}

} // namespace simcoon
