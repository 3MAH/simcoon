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
///@brief Closest-point projection return mapping. Contract and algorithm in the .hpp.

#include <cmath>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

namespace {

// Fischer-Burmeister residual of (Phi, Dl) with the |diag(B)| scaling of Fischer_Burmeister_m:
// the function returns the residual of the ENTERING iterate before stepping, so copies are
// handed in and the step discarded (N <= 3: one tiny solve).
double fb_merit(const vec &Phi, const vec &Y_crit, const mat &B, const vec &Dl) {
    vec Dl_copy = Dl, dDl;
    double err = 0.;
    Fischer_Burmeister_m(Phi, Y_crit, B, Dl_copy, dDl, err);
    return err;
}

} // namespace

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

    // full=false skips the state callbacks (hooks.K, flow_state_coupling, dLambda_dsigma): with
    // Dlambda = 0 the tangent assembly masks every mechanism (Ds_j > iota gate), so their values
    // cannot reach the returned operator.
    auto fill_tangent_pieces = [&](const vec &sig, const vec &Dl, bool full) {
        mat K = (full && hooks.K) ? hooks.K(sig, Dl) : zeros(N, N);
        r.Bhat_continuum = zeros(N, N);
        // kappa_eff = kappa + c: the multiplier-side state chain enters the tangent through the flux.
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

    // Elastic guard: admissible trial state -> return it (bit-identical elasticity across modes).
    if (Phi.max() <= 0.) {
        fill_tangent_pieces(r.sigma, r.Dlambda, false);
        r.converged = true;
        return r;
    }

    const mat I6 = eye(6, 6);
    const int snap_niter = 20;   // plain Newton gets this many iterations before the snap may act

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
        if (!inv(Minv, M)) return r;   // condensation breakdown: caller step-cuts

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

        // Merit on the TRUE residuals (Phi, R_sigma) at the current iterate.
        const double err = fb_merit(Phi, Y_crit, B_red, r.Dlambda) + norm(R_sigma, 2) / sigma_ref;
        r.error_history.push_back(err);
        r.error = err;
        if (err < precision) {
            fill_tangent_pieces(r.sigma, r.Dlambda, true);
            r.converged = true;
            return r;
        }

        // Semi-smooth Newton step on the reduced multiplier system.
        vec Dl_fb = r.Dlambda;
        vec dDl = zeros(N);
        double err_fb_unused = 0.;
        Fischer_Burmeister_m(Phi_red, Y_crit, B_red, Dl_fb, dDl, err_fb_unused);

        vec dsigma = -Minv * R_sigma;
        for (int j = 0; j < N; j++) dsigma -= dDl(j) * (Minv * (kappa[j] + c[j]));

        // Backtracking on the merit: semi-smooth iterations are locally non-monotone, so the
        // step is halved only when the merit grows by more than x2.
        const vec sigma_old = r.sigma;
        const vec Dl_old = r.Dlambda;
        double scale = 1.;
        int bt = 0;
        while (true) {
            r.sigma = sigma_old + scale * dsigma;
            r.Dlambda = clamp(Dl_old + scale * dDl, 0., datum::inf);
            // Deactivation snap (see the header): exact inactivity for a row hovering at the
            // complementarity boundary, gated on a stalled iteration.
            if (r.niter >= snap_niter) {
                for (int l = 0; l < N; l++) {
                    if ((Phi(l) < 0.) && (r.Dlambda(l) * fabs(B_red(l, l)) < 1.e-2 * Y_crit(l))) {
                        r.Dlambda(l) = 0.;
                    }
                }
            }
            if (!refresh(r.sigma, r.Dlambda)) return r;   // inner state solve failed
            eval_basic(r.sigma, r.Dlambda);
            if (r.sigma.has_nan() || Phi.has_nan()) return r;
            const double err_new = fb_merit(Phi, Y_crit, B_red, r.Dlambda) + norm(R_sigma, 2) / sigma_ref;
            if (err_new <= 2. * err || bt >= control.max_backtrack) break;
            scale *= 0.5;
            bt++;
        }
    }

    // maxiter reached: converged stays false, the caller applies the step cut.
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
        // Precondition not met: elastic operator as a sized placeholder (thermomechanical
        // callers read P_epsilon[l] unconditionally on the step-cut path).
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
