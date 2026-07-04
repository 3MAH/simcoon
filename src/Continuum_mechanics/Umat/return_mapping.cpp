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

///@file return_mapping.cpp
///@brief Closest-point projection return mapping (doc §cpp_return_mapping). Doc in the .hpp.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {



ReturnMappingResult closest_point_return_mapping(
    const vec &sigma_tr,
    const mat &L,
    const std::vector<ReturnMechanism> &mechanisms,
    const ReturnStateHooks &hooks,
    const vec &Y_crit,
    const ReturnMappingControl &control) {

    const int N = int(mechanisms.size());
    const int maxiter = (control.maxiter > 0) ? control.maxiter : simcoon::maxiter_umat;
    const double precision = (control.precision > 0.) ? control.precision : simcoon::precision_umat;

    ReturnMappingResult r;
    r.sigma = sigma_tr;
    r.Dlambda = zeros(N);

    const double sigma_ref = std::max(norm(sigma_tr, 2), Y_crit.min());

    auto refresh = [&](const vec &sig, const vec &Dl) -> bool {
        if (hooks.update_state) return hooks.update_state(sig, Dl);
        return true;
    };

    // Per-iterate evaluations at a state-consistent point.
    vec Phi(N);
    std::vector<vec> n_l(N), Lambda_j(N), kappa(N);
    std::vector<mat> D_j(N);
    vec R_sigma(6);

    auto eval_basic = [&](const vec &sig, const vec &Dl) {
        R_sigma = sig - sigma_tr;
        for (int j = 0; j < N; j++) {
            Phi(j) = mechanisms[j].Phi(sig);
            n_l[j] = mechanisms[j].dPhi_dsigma(sig);
            Lambda_j[j] = mechanisms[j].Lambda ? mechanisms[j].Lambda(sig) : n_l[j];
            kappa[j] = L * Lambda_j[j];
            R_sigma += Dl(j) * kappa[j];
        }
    };

    if (!refresh(r.sigma, r.Dlambda)) return r;
    eval_basic(r.sigma, r.Dlambda);

    // full=false skips the expensive FD callbacks (hooks.K, flow_state_coupling,
    // dLambda_dsigma): with Dlambda = 0 the tangent assembly masks every mechanism
    // (Ds_j > iota gate) so those values cannot influence the returned operator.
    auto fill_tangent_pieces = [&](const vec &sig, const vec &Dl, bool full) {
        mat K = (full && hooks.K) ? hooks.K(sig, Dl) : zeros(N, N);
        r.Bhat_continuum = zeros(N, N);
        // Effective flux kappa_eff = L:Lambda_tilde = kappa + c (doc eq:Lambda_tilde_state):
        // the multiplier-side state chain c^j enters the consistent tangent through kappa.
        r.kappa_j.assign(kappa.begin(), kappa.end());
        if (full && hooks.flow_state_coupling) {
            const std::vector<vec> c = hooks.flow_state_coupling(sig, Dl);
            for (int j = 0; j < N && j < int(c.size()); j++) r.kappa_j[j] += c[j];
        }
        r.dPhidsigma_l.assign(n_l.begin(), n_l.end());
        r.dLambda_dsigma_l.resize(N);
        for (int j = 0; j < N; j++) {
            r.dLambda_dsigma_l[j] = (full && mechanisms[j].dLambda_dsigma)
                                        ? mechanisms[j].dLambda_dsigma(sig)
                                        : zeros(6, 6);
        }
        for (int l = 0; l < N; l++) {
            for (int j = 0; j < N; j++) {
                r.Bhat_continuum(l, j) = sum(n_l[l] % r.kappa_j[j]) - K(l, j);
            }
        }
    };

    // Elastic guard: trial state admissible -> return it (bit-identical elasticity).
    // Cheap fill: all mechanisms are inactive, the expensive FD pieces are masked out.
    if (Phi.max() <= 0.) {
        fill_tangent_pieces(r.sigma, r.Dlambda, false);
        r.converged = true;
        return r;
    }

    const mat I6 = eye(6, 6);
    // Stall thresholds for the boundary-hovering safeguards below: plain semi-smooth
    // Newton must fail to converge for this many iterations before either kicks in.
    const int snap_niter = 20;
    const int stall_niter = 25;

    for (r.niter = 0; r.niter < maxiter; r.niter++) {

        // Newton ingredients at the current (state-consistent) iterate.
        mat K = hooks.K ? hooks.K(r.sigma, r.Dlambda) : zeros(N, N);
        mat M = I6;
        for (int j = 0; j < N; j++) {
            D_j[j] = mechanisms[j].dLambda_dsigma ? mechanisms[j].dLambda_dsigma(r.sigma)
                                                  : zeros(6, 6);
            M += r.Dlambda(j) * (L * D_j[j]);
        }
        mat Minv;
        if (!inv(Minv, M)) {
            if (std::getenv("SIMCOON_CPP_DEBUG")) std::fprintf(stderr, "[CPP fail] M breakdown niter=%d\n", r.niter);
            return r;  // Jacobian breakdown
        }

        std::vector<vec> c(N, zeros(6));
        if (hooks.flow_state_coupling) c = hooks.flow_state_coupling(r.sigma, r.Dlambda);

        const vec MinvR = Minv * R_sigma;
        vec Phi_red(N);
        mat B_red(N, N);
        for (int l = 0; l < N; l++) {
            Phi_red(l) = Phi(l) - sum(n_l[l] % MinvR);
            for (int j = 0; j < N; j++) {
                B_red(l, j) = -sum(n_l[l] % (Minv * (kappa[j] + c[j]))) + K(l, j);
            }
        }

        // Combined error on the TRUE residuals (Phi, R_sigma) at the current iterate.
        const double err_fb = Fischer_Burmeister_residual(Phi, r.Dlambda, B_red, Y_crit);
        const double err_R = norm(R_sigma, 2) / sigma_ref;
        const double err = err_fb + err_R;
        r.error_history.push_back(err);
        r.error = err;
        if (err < precision) {
            fill_tangent_pieces(r.sigma, r.Dlambda, true);
            r.converged = true;
            return r;
        }
        // Stall-guarded loose acceptance: semi-smooth iterations can hover at the
        // complementarity boundary (Phi slightly negative with a vanishing multiplier and
        // R_sigma at machine precision) without ever reaching the tight tolerance. Once the
        // iteration has demonstrably stalled (>= 25 iterations), accept when the stress
        // residual is converged and the KKT violation is below 1% of the criterion scale --
        // the corresponding stress error is orders of magnitude below any FE-level tolerance.
        if ((r.niter >= stall_niter) && (err_R < precision) && (err_fb < 1.e-2)) {
            bool snapped = false;
            for (int l = 0; l < N; l++) {
                if ((Phi(l) < 0.) && (r.Dlambda(l) * fabs(B_red(l, l)) < 1.e-2 * Y_crit(l))) {
                    r.Dlambda(l) = 0.;   // final cleanup of boundary-hovering multipliers
                    snapped = true;
                }
            }
            // Keep the returned (sigma, Dlambda, internal state) triple consistent:
            // re-solve the state at the cleaned multipliers before assembling the tangent.
            if (snapped && !refresh(r.sigma, r.Dlambda)) return r;
            if (snapped) eval_basic(r.sigma, r.Dlambda);
            fill_tangent_pieces(r.sigma, r.Dlambda, true);
            r.converged = true;
            return r;
        }

        // Semi-smooth Newton update on the reduced multiplier system.
        vec Dl_fb = r.Dlambda;
        vec dDl = zeros(N);
        double err_fb_unused = 0.;
        Fischer_Burmeister_m(Phi_red, Y_crit, B_red, Dl_fb, dDl, err_fb_unused);

        vec dsigma = -Minv * R_sigma;
        for (int j = 0; j < N; j++) dsigma -= dDl(j) * (Minv * (kappa[j] + c[j]));

        // Newton-step cap (joint scaling keeps the Newton direction).
        double scale = 1.;
        const double nds = norm(dsigma, 2);
        if (nds > control.max_dsigma_frac * sigma_ref) {
            scale = control.max_dsigma_frac * sigma_ref / nds;
        }

        // Apply with backtracking on the combined error (semi-smooth iterations may be
        // locally non-monotone; only cut when the error grows by more than x2).
        const vec sigma_old = r.sigma;
        const vec Dl_old = r.Dlambda;
        int bt = 0;
        while (true) {
            r.sigma = sigma_old + scale * dsigma;
            r.Dlambda = clamp(Dl_old + scale * dDl, 0., datum::inf);
            // Deactivation snap: a mechanism sitting at the complementarity boundary
            // (admissible Phi < 0 with a vanishing multiplier) is EXACTLY inactive;
            // snapping Dlambda to zero zeroes its Fischer-Burmeister residual and
            // prevents an asymptotic grind to maxiter at onset/unloading points.
            // Threshold 1e-2*Y_crit on Dl*|B| (~1% of the criterion scale in stress terms);
            // GATED on a stalled iteration (niter >= 20) so a genuinely active mechanism
            // with a small per-increment multiplier -- whose Phi can transiently dip
            // negative under fine time-stepping -- always gets the plain Newton first
            // (Q-quadratic solves finish well within 20 iterations).
            if (r.niter >= snap_niter) {
                for (int l = 0; l < N; l++) {
                    if ((Phi(l) < 0.) && (r.Dlambda(l) * fabs(B_red(l, l)) < 1.e-2 * Y_crit(l))) {
                        r.Dlambda(l) = 0.;
                    }
                }
            }
            if (!refresh(r.sigma, r.Dlambda)) {
                if (std::getenv("SIMCOON_CPP_DEBUG")) std::fprintf(stderr, "[CPP fail] inner state niter=%d\n", r.niter);
                return r;  // inner state solve failed
            }
            eval_basic(r.sigma, r.Dlambda);
            if (r.sigma.has_nan() || Phi.has_nan()) {
                if (std::getenv("SIMCOON_CPP_DEBUG")) std::fprintf(stderr, "[CPP fail] NaN niter=%d\n", r.niter);
                return r;
            }
            const double err_new = Fischer_Burmeister_residual(Phi, r.Dlambda, B_red, Y_crit) + norm(R_sigma, 2) / sigma_ref;
            if (err_new <= 2. * err || bt >= control.max_backtrack) break;
            scale *= 0.5;
            bt++;
        }
    }

    // maxiter reached: r.converged stays false; caller applies the tnew_dt step-cut.
    if (std::getenv("SIMCOON_CPP_DEBUG")) {
        std::fprintf(stderr, "[CPP fail] maxiter: err=%.3e Phi_max=%.3e |Dl|=%.3e |R|=%.3e\n",
                     r.error, Phi.max(), norm(r.Dlambda, 2), norm(R_sigma, 2));
    }
    return r;
}

ReturnMappingResult closest_point_return_mapping(
    const vec &sigma_tr,
    const mat &L,
    const ReturnMechanism &mechanism,
    const std::function<bool(const vec &, double)> &update_state,
    const std::function<double(const vec &, double)> &K_scalar,
    double Y_crit,
    const std::function<vec(const vec &, double)> &flow_state_coupling,
    const ReturnMappingControl &control) {

    std::vector<ReturnMechanism> mechs = {mechanism};
    ReturnStateHooks hooks;
    if (update_state) {
        hooks.update_state = [&update_state](const vec &sig, const vec &Dl) {
            return update_state(sig, Dl(0));
        };
    }
    hooks.K = [&K_scalar](const vec &sig, const vec &Dl) {
        mat K(1, 1);
        K(0, 0) = K_scalar ? K_scalar(sig, Dl(0)) : 0.;
        return K;
    };
    if (flow_state_coupling) {
        hooks.flow_state_coupling = [&flow_state_coupling](const vec &sig, const vec &Dl) {
            return std::vector<vec>{flow_state_coupling(sig, Dl(0))};
        };
    }
    vec Y(1);
    Y(0) = Y_crit;
    return closest_point_return_mapping(sigma_tr, L, mechs, hooks, Y, control);
}

ContinuumTangent cpp_consistent_tangent(const ReturnMappingResult &r, const mat &L) {
    if (!r.converged || r.kappa_j.empty()) {
        // Precondition not met: return the elastic operator as a safe placeholder.
        // P_epsilon must be SIZED (zeros): thermomechanical UMATs read P_epsilon[l]
        // unconditionally on this path.
        ContinuumTangent ct;
        ct.Lt = L;
        const uword Nm = std::max(r.Dlambda.n_elem, uword(1));
        ct.invBhat = zeros(Nm, Nm);
        ct.P_epsilon.assign(Nm, zeros(6));
        return ct;
    }
    return assemble_algorithmic_tangent(r.Bhat_continuum, r.kappa_j, r.dPhidsigma_l,
                                        r.Dlambda, L, r.dLambda_dsigma_l);
}

} // namespace simcoon
