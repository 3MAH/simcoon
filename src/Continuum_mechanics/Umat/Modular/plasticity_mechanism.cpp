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

#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/small_lu.hpp>
#include <simcoon/parameter.hpp>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace simcoon {

namespace {

// state_drift denominator floor: rejection at near-zero Dp needs a flow-rule error > tol * 1e-4 strain
constexpr double drift_scale_floor = 1.0e-4;

// HT = H.each_row() % T.t(), A = I + gamma HT in the expressions' order: bit-identical unless FMA crosses statements
constexpr double T_shear[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};   // Ir05(): engineering -> tensorial shear

void fill_HT_and_A(const arma::mat& H, double gamma, arma::mat::fixed<6,6>& HT, arma::mat::fixed<6,6>& A) {
    for (arma::uword j = 0; j < 6; ++j) {
        for (arma::uword i = 0; i < 6; ++i) {
            HT(i, j) = H(i, j) * T_shear[j];
            const double g = HT(i, j) * gamma;
            A(i, j) = ((i == j) ? 1.0 : 0.0) + g;
        }
    }
}
}  // namespace

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

PlasticityMechanism::PlasticityMechanism(
    YieldType yield_type,
    IsoHardType iso_type,
    KinHardType kin_type,
    const arma::vec& props,
    int N_iso,
    int N_kin
)
    : PlasticityMechanism(yield_type, iso_type, kin_type, N_iso, N_kin)
{
    int offset = 0;
    configure(props, offset);
    register_variables();
}

// ========== Configuration ==========

void PlasticityMechanism::configure(const arma::vec& props, int& offset) {

    sigma_Y_ = props(offset);
    offset += 1;

    if (!yield_->is_configured()) {
        yield_->configure(yield_type_, props, offset);
    }

    // Configure hardening models
    iso_hard_->configure(props, offset);
    iso_hard_->set_yield_stress(sigma_Y_);
    kin_hard_->configure(props, offset);
}

void PlasticityMechanism::set_increment(double theta_start, double theta_end, double DTime) {
    theta_start_ = theta_start;
    theta_end_ = theta_end;
    DTime_ = DTime;
    iso_hard_->set_increment(theta_start, theta_end, DTime);
}

void PlasticityMechanism::predict(const arma::vec& /*E_total_end*/, double /*DTime*/) {
    if (has_pdot_) {
        ivc_.at(pdot_index_).scalar() = 0.0;   // set again by update / refresh_state when the increment is plastic
    }
}

void PlasticityMechanism::hardening_stress(double p, double Dp, double& R, double& dR_dp, double& dR_dDp) const {
    if (has_pdot_) {
        iso_hard_->R_rate(p, Dp, R, dR_dp, dR_dDp);
    } else {
        iso_hard_->R_and_dR(p, R, dR_dp);
        dR_dDp = 0.0;
    }
}

bool PlasticityMechanism::predict_multipliers(const arma::vec& sigma_eff, const arma::mat& L, arma::vec& ds, int offset) const {
    if (!has_pdot_) return false;
    double sigma_eq;
    arma::vec::fixed<6> n;
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    if (has_kin) {
        yield_->equivalent_stress_and_flow(arma::vec::fixed<6>(sigma_eff - total_backstress_voigt()), sigma_eq, n);
    } else {
        yield_->equivalent_stress_and_flow(sigma_eff, sigma_eq, n);
    }
    const double nLn = arma::dot(n, L * n);
    const double x0 = iso_hard_->multiplier_predictor(ivc_.at(p_index_).scalar(), sigma_eq, nLn);
    if (x0 <= 0.0) return false;
    ds(offset) = x0;
    return true;
}

void PlasticityMechanism::register_variables() {
    ivc_.add_scalar("p", 0.0);
    ivc_.add_vec("EP", arma::zeros(6), true);
    kin_hard_->register_variables(ivc_);
    has_pdot_ = iso_hard_->rate_dependent();
    if (has_pdot_) {
        ivc_.add_scalar("pdot", 0.0);   // the plastic strain rate of the increment, the EPJCK layout
    }
    // resolved once: the return-mapping loop reads them through ivc_.at()
    p_index_ = ivc_.index_of("p");
    EP_index_ = ivc_.index_of("EP");
    if (has_pdot_) {
        pdot_index_ = ivc_.index_of("pdot");
    }
}

// ========== Yield-side queries ==========

arma::vec PlasticityMechanism::total_backstress_voigt() const {
    kin_hard_->compute_backstresses(ivc_, X_branches_);
    arma::vec::fixed<6> X(arma::fill::zeros);
    for (arma::uword i = 0; i < X_branches_.n_cols; ++i) {
        X += X_branches_.col(i);
    }
    return X;
}

double PlasticityMechanism::equivalent_stress(const arma::vec& sigma) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->equivalent_stress(sigma);
    }
    return yield_->equivalent_stress(sigma - total_backstress_voigt());
}

double PlasticityMechanism::equivalent_stress(const tensor2& sigma) const {
    return equivalent_stress(sigma.to_arma_voigt());
}

double PlasticityMechanism::yield_function(const arma::vec& sigma) const {
    const auto& p_var = ivc_.at(p_index_);
    double R, dR_dp, dR_dDp;
    hardening_stress(p_var.scalar(), p_var.delta_scalar(), R, dR_dp, dR_dDp);
    return equivalent_stress(sigma) - R - sigma_Y_;
}

double PlasticityMechanism::yield_function(const tensor2& sigma) const {
    return yield_function(sigma.to_arma_voigt());
}

arma::vec PlasticityMechanism::flow_direction(const arma::vec& sigma) const {
    if (kin_hard_->num_backstresses() == 0) {
        return yield_->flow_direction(sigma);
    }
    return yield_->flow_direction(sigma - total_backstress_voigt());
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
    double p = ivc_.at(p_index_).scalar();
    hessian_valid_ = false;

    // X_i built once here and reused by hardening_modulus: valid within this FB iteration only
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    double sigma_eq;
    if (has_kin) {
        const arma::vec::fixed<6> sigma_eff = sigma - total_backstress_voigt();
        yield_->equivalent_stress_and_flow(sigma_eff, sigma_eq, flow_dir_);
    } else {
        yield_->equivalent_stress_and_flow(sigma, sigma_eq, flow_dir_);
    }

    double R, dR_dp, dR_dDp;
    hardening_stress(p, ivc_.at(p_index_).delta_scalar(), R, dR_dp, dR_dDp);
    // Yield function: Phi = sigma_eq - R - sigma_Y
    Phi.set_size(1);
    Phi(0) = sigma_eq - R - sigma_Y_;

    Y_crit.set_size(1);
    Y_crit(0) = std::max(sigma_Y_, 1.0);

    if (L.n_rows == 6 && L.n_cols == 6 && flow_dir_.n_elem == 6) {
        mul6(L, flow_dir_, kappa_);
    } else {
        kappa_ = L * flow_dir_;
    }
    H_total_ = dR_dp + dR_dDp + (has_kin
        ? kin_hard_->hardening_modulus(flow_dir_, X_branches_)
        : 0.0);
}

const arma::mat& PlasticityMechanism::dPhi_dsigma(
    const arma::vec& /*sigma*/) const {
    // Requires compute_constraints to have been called this FB iteration.
    dPhi_column_.col(0) = flow_dir_;
    return dPhi_column_;
}

const arma::mat& PlasticityMechanism::kappa(
    const arma::vec& /*sigma*/, double /*DT*/, const arma::mat& /*L_ref*/) const {
    // kappa = L_ref · n, already computed by compute_constraints.
    kappa_column_.col(0) = kappa_;
    return kappa_column_;
}

const arma::mat* PlasticityMechanism::dLambda_dsigma(
    const arma::vec& sigma) const {
    if (!yield_->has_flow_hessian()) {
        return nullptr;
    }
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    const arma::vec sigma_eff =
        has_kin ? arma::vec(sigma - total_backstress_voigt()) : sigma;
    hessian_eng_ = yield_->flow_hessian(sigma_eff);
    hessian_at_ = sigma;
    hessian_valid_ = true;
    return &hessian_eng_;
}

bool PlasticityMechanism::supports_closest_point() const {
    return yield_->has_flow_hessian();
}

bool PlasticityMechanism::refresh_state(
    const arma::vec& sigma,
    const arma::vec& Ds_total,
    int offset) {
    const double dp = std::max(Ds_total(offset), 0.0);   // the integrator projects Ds >= 0 too
    const bool has_kin = kin_hard_->num_backstresses() > 0;
    implicit_update_ = true;

    auto& p_var = ivc_.at(p_index_);
    p_var.scalar() = p_var.scalar_start() + dp;

    // Newton on the n <-> X coupling (see header); flow and Hessian of the final xi feed the derivatives below
    arma::vec n(6), xi(6);
    double gamma = 0.0;
    arma::vec beta = arma::zeros(6);
    double sigma_eq = 0.0;
    arma::mat::fixed<6,6> H;
    if (!has_kin) {
        xi = sigma;
        yield_->equivalent_stress_flow_hessian(xi, sigma_eq, n, H);
    } else {
        // Start = normal of the dp-relaxed start backstress: exact for von Mises (zero steps);
        // warm-starting from the previous normal measured worse (J2 pays real steps).
        kin_hard_->refresh_state(dp, arma::vec::fixed<6>(arma::fill::zeros), ivc_);   // alpha_i^n / (1 + D_i dp)
        n = yield_->flow_direction(sigma - total_backstress_voigt());
        bool converged = false;
        for (int it = 0; it < 25; ++it) {
            kin_hard_->refresh_state(dp, n, ivc_);
            xi = sigma - total_backstress_voigt();
            arma::vec flow;
            yield_->equivalent_stress_flow_hessian(xi, sigma_eq, flow, H);
            const arma::vec G = n - flow;
            if (arma::norm(G, 2) <= 1e-12 * std::max(1.0, arma::norm(n, 2))) {
                converged = true;
                break;
            }
            kin_hard_->backward_euler_factors(dp, n, ivc_, gamma, beta);
            arma::mat::fixed<6,6> HT, A;
            fill_HT_and_A(H, gamma, HT, A);   // I + gamma H T
            small_lu<6> lu_A;   // on the stack: a LAPACK 6x6 solve per step otherwise
            if (!lu_A.factor(A)) {
                return false;
            }
            n -= lu_A.solve(G);
        }
        if (!converged) {
            return false;
        }
        kin_hard_->backward_euler_factors(dp, n, ivc_, gamma, beta);
    }

    auto& EP_var = ivc_.at(EP_index_);
    EP_var.raw_voigt() = EP_var.raw_voigt_start() + dp * n;
    if (has_pdot_) {
        ivc_.at(pdot_index_).scalar() = (DTime_ > simcoon::iota && dp > simcoon::iota) ? dp / DTime_ : 0.0;
    }
    // Total derivatives at this state (class note in the header).
    ClosestPointIngredients& out = cpp_cache_[0];
    double R_p, dR_dp, dR_dDp;
    hardening_stress(p_var.scalar(), dp, R_p, dR_dp, dR_dDp);
    dR_dp += dR_dDp;   // K = dPhi/dDp at fixed sigma: the rate sensitivity sits with the hardening modulus
    out.Phi = sigma_eq - R_p - sigma_Y_;
    out.Lambda = n;
    if (!has_kin) {
        out.dPhi_dsigma = n;
        out.dLambda_dsigma = H;
        out.dLambda_dDs = arma::zeros(6);
        out.K = -dR_dp;
        return true;
    }
    arma::mat::fixed<6,6> HT, A;
    fill_HT_and_A(H, gamma, HT, A);
    small_lu<6> lu_A;   // A = I + gamma H T on the stack: solves instead of a LAPACK inverse
    if (!lu_A.factor(A)) {
        return false;
    }
    arma::vec::fixed<6> Tn;
    for (arma::uword i = 0; i < 6; ++i) Tn(i) = T_shear[i] * n(i);
    out.dLambda_dsigma = lu_A.solve(H);
    out.dPhi_dsigma = n - gamma * mul6_At(out.dLambda_dsigma, Tn);
    out.K = -arma::dot(Tn, lu_A.solve(beta)) - dR_dp;
    out.dLambda_dDs = -lu_A.solve(arma::vec(HT * beta));
    return true;
}

double PlasticityMechanism::state_drift(const arma::vec& sigma) const {
    // Committed Dp from the (projected) internal variable, not the raw FB
    // row — update() enforces p >= p_start, the state is the truth.
    const double dp = ivc_.at(p_index_).delta_scalar();
    const arma::vec dEP = ivc_.at(EP_index_).delta_vec();
    if (dp < simcoon::iota && norm_strain(dEP) < simcoon::iota) {
        // inert increment: avoid evaluating the flow direction at a
        // possibly-zero stress
        return 0.0;
    }
    // norm_strain, not arma::norm (sqrt2-heavy shears); the floor keeps a negligible EP residue at Dp ~ 0 inert
    const arma::vec n = flow_direction(sigma);
    return norm_strain(dEP - dp * n)
         / std::max(dp * norm_strain(n), drift_scale_floor);
}

double PlasticityMechanism::drift_tolerance() const {
    return yield_->has_stable_flow_direction()
        ? StrainMechanism::drift_tolerance()
        : std::numeric_limits<double>::infinity();
}

double PlasticityMechanism::multiplier_cap() const {
    return 1.0;
}

double PlasticityMechanism::committed_multiplier() const {
    return ivc_.at(p_index_).delta_scalar();
}

void PlasticityMechanism::compute_jacobian_contribution(
    const arma::vec& sigma,
    const arma::mat& L,
    arma::mat& B,
    int row_offset
) const {
    // K = dPhi/dp = -H_total, so the mode-1 assembly Bhat = -B = n:kappa + H_total holds exactly

    // dPhi/dsigma * L * dPhi/dsigma (engineering Voigt dot = n : L : n)
    double dPhi_L_dPhi = arma::dot(flow_dir_, kappa_);

    B(row_offset, row_offset) = -dPhi_L_dPhi - H_total_;
}

const arma::vec& PlasticityMechanism::inelastic_strain() const {
    return ivc_.at(EP_index_).raw_voigt();
}

void PlasticityMechanism::update(
    const arma::vec& ds,
    int offset
) {
    // projected onto Dp >= 0 (see update() in the header)
    double dp = ds(offset);
    implicit_update_ = false;
    if (dp == 0.0) {
        return;   // the update below is exactly the identity for dp = 0
    }

    // Update accumulated plastic strain (projected: p never below p_start)
    double& p = ivc_.at(p_index_).scalar();
    const double p_start = ivc_.at(p_index_).scalar_start();
    if (p + dp < p_start) {
        dp = p_start - p;
    }
    p += dp;
    // Update plastic strain
    arma::vec& EP = ivc_.at(EP_index_).raw_voigt();
    EP += dp * flow_dir_;
    // Update kinematic hardening variables
    kin_hard_->update(dp, flow_dir_, ivc_);
    if (has_pdot_) {
        const double Dp = p - p_start;
        ivc_.at(pdot_index_).scalar() = (DTime_ > simcoon::iota && Dp > simcoon::iota) ? Dp / DTime_ : 0.0;
    }
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
        // Lt = L - (L n) x (L n) / H_eff, H_eff = n:L:n + H_total

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
    const arma::vec& E_start,
    const arma::vec& E_end,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d
) const {
    // Get increments
    double Dp = ivc_.at(p_index_).delta_scalar();
    const auto& EP_var = ivc_.at(EP_index_);
    if (Dp == 0.0) {
        // in place, before building DEP: a NaN difference counts as a move (as is_zero would)
        bool moved = false;
        for (arma::uword i = 0; i < 6 && !moved; ++i) {
            moved = (EP_var.raw_voigt()(i) - EP_var.raw_voigt_start()(i)) != 0.0;
        }
        if (!moved) {
            Wm_r = Wm_ir = Wm_d = 0.0;   // both dots below would be (signed) zeros
            return;
        }
    }
    const arma::vec DEP = EP_var.delta_vec();

    // stored in the hardening: the trapezoids in p (irrecoverable) and in the back-strains (recoverable)
    const double p_start = ivc_.at(p_index_).scalar_start(), p = ivc_.at(p_index_).scalar();
    Wm_ir = 0.5 * (iso_hard_->stored_force(p_start, theta_start_) + iso_hard_->stored_force(p, theta_end_)) * Dp;
    Wm_r = kin_hard_->stored_increment(ivc_);
    Wm_d = 0.5 * arma::dot(sigma_start + sigma, DEP) - Wm_ir - Wm_r;
}

void PlasticityMechanism::dissipation_sensitivity(const DissipationInputs& in, arma::vec& dDgamma_dE,
                                                  double& dDgamma_dT) const {
    if (has_pdot_) {
        throw std::invalid_argument("the thermomechanical face of MODUL does not serve a rate- and temperature-dependent "
                                    "hardening (Johnson-Cook) yet: its heat needs the hardening entropy; EPJCK keeps its "
                                    "dedicated thermomechanical kernel");
    }
    const double Dp = ivc_.at(p_index_).delta_scalar();
    const arma::vec DEP = ivc_.at(EP_index_).delta_vec();
    const double p_start = ivc_.at(p_index_).scalar_start(), p = ivc_.at(p_index_).scalar();

    // the direction the update moved along (flow_dir_: the last compute_constraints), X_i at both ends
    const int n_kin = kin_hard_->num_backstresses();
    arma::mat X_i_start(6, n_kin), X_i(6, n_kin);
    arma::vec::fixed<6> X_start(arma::fill::zeros), X(arma::fill::zeros);
    for (int i = 0; i < n_kin; ++i) {
        const InternalVariable& a = *kin_hard_->back_strain_variable(ivc_, i);
        X_i_start.col(i) = kin_hard_->backstress_voigt(i, a.raw_voigt_start());
        X_i.col(i) = kin_hard_->backstress_voigt(i, a.raw_voigt());
        X_start += X_i_start.col(i);
        X += X_i.col(i);
    }
    const arma::vec& Lambda = flow_dir_;

    // d(Dgamma)/d(Dp) at fixed direction: the hardening forces A_p = -R, A_a = -X of the dedicated kernels
    double R_p, dR_p;
    iso_hard_->R_and_dR(p, R_p, dR_p);
    double f_p = 0.5 * (arma::dot(in.sigma_start + in.sigma, Lambda) - iso_hard_->R(p_start) - R_p - dR_p * Dp);
    double g_X = 0.0;
    // turning of the direction over the increment on xi = sigma_eff - X, xi_{n+1} : dLambda/dxi = 0 (degree 0)
    arma::vec g_turn(6, arma::fill::zeros);
    if (Dp > 0.0 && yield_->has_flow_hessian()) {
        // the Hessian of dLambda_dsigma when it was taken at this stress (one evaluation per increment)
        const bool cached = hessian_valid_
            && std::memcmp(hessian_at_.memptr(), in.sigma_eff.memptr(), 6 * sizeof(double)) == 0;
        // H (sigma_n - X_n + sigma_n+1 - X_n+1), H xi_eff,n+1 = 0: the committed end contributes sigma - sigma_eff (zero without damage)
        const arma::vec::fixed<6> dxi = (in.sigma_start - X_start) + (in.sigma - in.sigma_eff);
        if (cached) {
            g_turn = 0.5 * Dp * (hessian_eng_ * dxi);
        } else {
            g_turn = 0.5 * Dp * (yield_->flow_hessian(arma::vec(in.sigma_eff - X)) * dxi);
        }
    }
    for (int i = 0; i < n_kin; ++i) {
        const InternalVariable& a = *kin_hard_->back_strain_variable(ivc_, i);
        const arma::vec Da = a.delta_vec();
        // d alpha_{n+1} / d Dp at fixed direction: Lambda (Prager); with recovery, (Lambda - D alpha_{n+1}) / (1 + D Dp)
        // of the backward-Euler update, the secant Da / Dp over the explicit cutting-plane iterates (measured: the
        // consistent form of the last iterate is 5x worse there); the back stress follows the flow, dxi = dsigma_eff - X(Lambda) dDp
        const double D_i = kin_hard_->recovery_parameter(i);
        const arma::vec Lambda_i = (D_i == 0.0 || Dp <= 0.0) ? Lambda
            : implicit_update_ ? arma::vec((Lambda - D_i * a.raw_voigt()) / (1.0 + D_i * Dp)) : arma::vec(Da / Dp);
        const arma::vec X_i_of_Lambda = kin_hard_->backstress_voigt(i, Lambda_i);
        f_p -= 0.5 * (arma::dot(X_i_start.col(i) + X_i.col(i), Lambda_i) + arma::dot(X_i_of_Lambda, Da));
        g_X += arma::dot(g_turn, X_i_of_Lambda);
    }
    dDgamma_dE = 0.5 * mul6_At(in.dSdE, DEP) + (f_p - g_X) * in.P_epsilon.col(0) + mul6_At(in.dSdE_eff, g_turn);
    dDgamma_dT = 0.5 * arma::dot(in.dSdT, DEP) + (f_p - g_X) * in.P_theta(0) + arma::dot(in.dSdT_eff, g_turn);
}

} // namespace simcoon
