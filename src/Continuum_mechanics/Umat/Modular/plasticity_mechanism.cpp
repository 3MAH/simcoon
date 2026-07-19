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

void PlasticityMechanism::register_variables() {
    ivc_.add_scalar("p", 0.0);
    ivc_.add_vec("EP", arma::zeros(6), true);
    kin_hard_->register_variables(ivc_);
}

// ========== Yield-side queries ==========

double PlasticityMechanism::equivalent_stress(const arma::vec& sigma) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->equivalent_stress(sigma);
    }
    return yield_->equivalent_stress(sigma - kin_hard_->total_backstress(ivc_).to_arma_voigt());
}

double PlasticityMechanism::equivalent_stress(const tensor2& sigma) const {
    return equivalent_stress(sigma.to_arma_voigt());
}

double PlasticityMechanism::yield_function(const arma::vec& sigma) const {
    const double p = ivc_.get("p").scalar();
    return equivalent_stress(sigma) - iso_hard_->R(p) - sigma_Y_;
}

double PlasticityMechanism::yield_function(const tensor2& sigma) const {
    return yield_function(sigma.to_arma_voigt());
}

arma::vec PlasticityMechanism::flow_direction(const arma::vec& sigma) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->flow_direction(sigma);
    }
    return yield_->flow_direction(sigma - kin_hard_->total_backstress(ivc_).to_arma_voigt());
}

tensor2 PlasticityMechanism::flow_direction(const tensor2& sigma) const {
    return strain(flow_direction(sigma.to_arma_voigt()));
}

// ========== Constitutive Computations ==========

void PlasticityMechanism::compute_constraints(
    const arma::vec& sigma,
    const arma::vec& /*E_total*/,
    const arma::mat& L,
    double /*DTime*/,
    arma::vec& Phi,
    arma::vec& Y_crit
) const {
    // Get internal variables
    double p = ivc_.get("p").scalar();

    // Shifted stress (sigma - X); skip the subtraction under pure isotropic
    // hardening. The branch backstresses X_i are built ONCE here and reused
    // for hardening_modulus below — they are only valid within this FB
    // iteration (the back-strains move in update()).
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    double sigma_eq;
    if (has_kin) {
        kin_hard_->compute_backstresses(ivc_, X_branches_);
        tensor2 X_t = tensor2::zeros(VoigtType::stress);
        for (const auto& x : X_branches_) {
            X_t += x;
        }
        const arma::vec sigma_eff = sigma - X_t.to_arma_voigt();
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
    H_total_ = dR_dp + (has_kin
        ? kin_hard_->hardening_modulus(strain(flow_dir_), X_branches_)
        : 0.0);
}

const std::vector<tensor2>& PlasticityMechanism::dPhi_dsigma(
    const arma::vec& /*sigma*/) const {
    // Requires compute_constraints to have been called this FB iteration.
    dPhi_dsigma_cache_[0] = strain(flow_dir_);
    return dPhi_dsigma_cache_;
}

const std::vector<tensor2>& PlasticityMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/, const arma::mat& /*L_ref*/) const {
    // kappa = L_ref · n, already computed by compute_constraints.
    kappa_cache_[0] = stress(kappa_);
    return kappa_cache_;
}

const std::vector<tensor4>* PlasticityMechanism::dLambda_dsigma(
    const arma::vec& sigma) const {
    if (!yield_->has_flow_hessian()) {
        return nullptr;
    }
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    const arma::vec sigma_eff =
        has_kin ? arma::vec(sigma - kin_hard_->total_backstress(ivc_).to_arma_voigt()) : sigma;
    hessian_cache_[0] =
        tensor4(arma::mat(yield_->flow_hessian(sigma_eff)), Tensor4Type::compliance);
    return &hessian_cache_;
}

bool PlasticityMechanism::refresh_state(
    const arma::vec& sigma,
    const arma::vec& Ds_total,
    int offset) {
    const double dp = Ds_total(offset);

    auto& p_var = ivc_.get("p");
    p_var.scalar() = p_var.scalar_start() + dp;

    // n and the back-strains couple through X(alpha): fixed point (direct for
    // pure isotropic hardening; contraction ~ dp C / sigma_eq per iteration).
    const int n_fp = (kin_hard_->num_backstresses() > 0) ? 50 : 1;
    arma::vec n = arma::zeros(6);
    for (int it = 0; it < n_fp; ++it) {
        const arma::vec X = kin_hard_->total_backstress(ivc_).to_arma_voigt();
        const arma::vec n_new = yield_->flow_direction(sigma - X);
        const double dn = arma::norm(n_new - n, 2);
        n = n_new;
        kin_hard_->refresh_state(dp, strain(n), ivc_);
        if (it > 0 && dn < 1e-14) {
            break;
        }
    }

    auto& EP_var = ivc_.get("EP");
    EP_var.raw_voigt() = EP_var.raw_voigt_start() + dp * n;
    return true;
}

void PlasticityMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    arma::mat& B,
    int row_offset
) const {
    // B = -dPhi/dsigma : kappa + K with K = dPhi/dp = -H_total (hardening
    // REDUCES Phi as Dp grows) — the plastic_isotropic_ccp convention, so the
    // mode-1 assembly's Bhat = -B = n:kappa + H_total holds exactly.

    // dPhi/dsigma * L * dPhi/dsigma (engineering Voigt dot = n : L : n)
    double dPhi_L_dPhi = arma::dot(flow_dir_, kappa_);

    B(row_offset, row_offset) = -dPhi_L_dPhi - H_total_;
}

arma::vec PlasticityMechanism::inelastic_strain() const {
    return ivc_.get("EP").raw_voigt();
}

void PlasticityMechanism::update(
    const arma::vec& ds,
    int offset
) {
    // Apply the FB correction PROJECTED onto Dp >= 0: negative per-iteration
    // corrections walk back an overshoot (gating on dp > iota froze the state
    // after any overshoot — maxiter spins with a ~1 MPa yield violation and a
    // 12x slowdown), but the multiplier never drops below its start-of-
    // increment value. Without the projection a tiny negative Dp is a SPURIOUS
    // Fischer-Burmeister root whenever Phi > 0 is large (phi(Phi, Dp) ~ Dp for
    // |Dp| << Phi) — reached in practice for near-singular hardening slopes,
    // e.g. power-law R = k p^m with m < 1 at first yield (p < 0, Phi -> +100
    // MPa, mechanism inert).
    double dp = ds(offset);

    // Update accumulated plastic strain (projected: p never below p_start)
    double& p = ivc_.get("p").scalar();
    const double p_start = ivc_.get("p").scalar_start();
    if (p + dp < p_start) {
        dp = p_start - p;
    }
    p += dp;

    // Update plastic strain
    arma::vec& EP = ivc_.get("EP").raw_voigt();
    EP += dp * flow_dir_;

    // Update kinematic hardening variables
    kin_hard_->update(dp, strain(flow_dir_), ivc_);
}

void PlasticityMechanism::tangent_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    const arma::vec& Ds,
    int offset,
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
            // Lt -= (kappa ⊗ kappa) / H_eff — raw rank-1 update (bit-identical
            // to the typed dyadic; no eng↔Mandel round trip in the hot path)
            Lt -= (kappa_ * kappa_.t()) / H_eff;
        }
    }
}

void PlasticityMechanism::compute_work(
    const arma::vec& sigma_start,
    const arma::vec& sigma,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Get increments
    double Dp = ivc_.get("p").delta_scalar();
    arma::vec DEP = ivc_.get("EP").delta_vec();

    // Average stress during increment
    arma::vec sigma_avg = 0.5 * (sigma_start + sigma);

    // Get backstress (if any)
    const arma::vec X = kin_hard_->total_backstress(ivc_).to_arma_voigt();

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
