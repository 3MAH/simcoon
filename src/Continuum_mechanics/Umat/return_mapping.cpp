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
#include <stdexcept>
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

    // Tangent ingredients at the converged iterate from the pieces the Newton just built
    // (kappa_eff = kappa + c: the multiplier-side state chain enters through the flux).
    auto fill_tangent_pieces = [&](const mat &K, const std::vector<mat> &D, const std::vector<vec> &c) {
        r.Bhat_continuum.set_size(N, N);
        r.kappa_j.assign(kappa.begin(), kappa.end());
        for (int j = 0; j < N; j++) r.kappa_j[j] += c[j];
        r.dPhidsigma_l.assign(n_l.begin(), n_l.end());
        r.dLambda_dsigma_l = D;
        for (int l = 0; l < N; l++) {
            for (int j = 0; j < N; j++) {
                r.Bhat_continuum(l, j) = sum(n_l[l] % r.kappa_j[j]) - K(l, j);
            }
        }
    };

    // Elastic guard: admissible trial state -> return it (bit-identical elasticity across
    // modes). With Dlambda = 0 the tangent assembly masks every mechanism, so the state
    // callbacks are not needed.
    if (Phi.max() <= 0.) {
        fill_tangent_pieces(zeros(N, N), std::vector<mat>(N, zeros(6, 6)), std::vector<vec>(N, zeros(6)));
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

        std::vector<vec> c = hooks.flow_state_coupling ? hooks.flow_state_coupling(r.sigma, r.Dlambda)
                                                       : std::vector<vec>(N, zeros(6));

        const vec MinvR = Minv * R_sigma;
        std::vector<vec> Minv_kc(N);   // M^-1 (kappa_j + c_j), for B_red and the stress step
        for (int j = 0; j < N; j++) Minv_kc[j] = Minv * (kappa[j] + c[j]);
        vec Phi_red(N);
        mat B_red(N, N);
        for (int l = 0; l < N; l++) {
            Phi_red(l) = Phi(l) - sum(n_l[l] % MinvR);
            for (int j = 0; j < N; j++) {
                B_red(l, j) = -sum(n_l[l] % Minv_kc[j]) + K(l, j);
            }
        }

        // Merit on the TRUE residuals (Phi, R_sigma) at the current iterate.
        if (!B_red.is_finite()) return r;
        const double err = Fischer_Burmeister_residual(Phi, Y_crit, B_red, r.Dlambda) + norm(R_sigma, 2) / sigma_ref;
        r.error_history.push_back(err);
        r.error = err;
        if (err < precision) {
            fill_tangent_pieces(K, D_j, c);
            r.converged = true;
            return r;
        }

        // Semi-smooth Newton step on the reduced multiplier system.
        vec Dl_fb = r.Dlambda;
        vec dDl = zeros(N);
        double err_fb_unused = 0.;
        try {
            Fischer_Burmeister_m(Phi_red, Y_crit, B_red, Dl_fb, dDl, err_fb_unused);
        } catch (const std::exception &) {
            return r;   // LAPACK failure on the reduced system: non-convergence, the caller step-cuts
        }

        vec dsigma = -MinvR;
        for (int j = 0; j < N; j++) dsigma -= dDl(j) * Minv_kc[j];

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
            // A trial point whose inner state solve fails, or that is not finite, is a rejected
            // step like a merit blow-up: halve it (the inner Newton of a backstress row has no
            // solution past a dp where the relaxed backstress overtakes the stress).
            bool trial_ok = refresh(r.sigma, r.Dlambda);
            if (trial_ok) {
                eval_basic(r.sigma, r.Dlambda);
                trial_ok = r.sigma.is_finite() && Phi.is_finite();
            }
            const double err_new = trial_ok
                ? Fischer_Burmeister_residual(Phi, Y_crit, B_red, r.Dlambda) + norm(R_sigma, 2) / sigma_ref
                : datum::inf;
            if (trial_ok && err_new <= 2. * err) break;
            if (bt >= control.max_backtrack) {
                if (!trial_ok) return r;   // no admissible step within the backtracking budget
                break;
            }
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
