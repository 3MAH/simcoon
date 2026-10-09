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
#include <memory>
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

    // Elastic block: Lc is its tangent at the iterate (the given L for a linear block).
    mat Lc = L;
    const bool nonlinear_block = static_cast<bool>(hooks.elastic_response);
    auto eval_basic = [&](const vec &sig, const vec &Dl) {
        for (int j = 0; j < N; j++) {
            Phi(j) = mechanisms[j].Phi(sig);
            n_l[j] = mechanisms[j].dPhi_dsigma(sig);
            Lambda_j[j] = mechanisms[j].Lambda ? mechanisms[j].Lambda(sig) : n_l[j];
        }
        if (nonlinear_block) {
            vec eps_el = hooks.eps_el_tr;
            for (int j = 0; j < N; j++) eps_el -= Dl(j) * Lambda_j[j];
            vec sigma_F;
            hooks.elastic_response(eps_el, sigma_F, Lc);
            R_sigma = sig - sigma_F;
            for (int j = 0; j < N; j++) kappa[j] = Lc * Lambda_j[j];
        } else {
            R_sigma = sig - sigma_tr;
            for (int j = 0; j < N; j++) {
                kappa[j] = L * Lambda_j[j];
                R_sigma += Dl(j) * kappa[j];
            }
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
            M += r.Dlambda(j) * (Lc * D_j[j]);
        }
        mat Minv;
        if (!inv(Minv, M)) return r;   // condensation breakdown: caller step-cuts

        std::vector<vec> c(N, zeros(6));
        if (hooks.flow_state_coupling) {
            const std::vector<vec> dLambda = hooks.flow_state_coupling(r.sigma, r.Dlambda);
            for (int j = 0; j < N && j < int(dLambda.size()); j++) c[j] = Lc * dLambda[j];
        }

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

        // Semi-smooth Newton step on the reduced multiplier system, over the active set: a row
        // with Phi < 0, no multiplier and a vanishing diagonal (a reorientation surface at zero
        // effective stress) has an identically zero FB row and would make the solve singular;
        // its step is zero anyway.
        const uvec active = find((Phi >= 0.) + (r.Dlambda > 0.) + (abs(B_red.diag()) >= simcoon::iota));
        vec dDl = zeros(N);
        if (active.n_elem > 0) {
            vec Dl_fb = r.Dlambda(active);
            vec dDl_a(active.n_elem);
            double err_fb_unused = 0.;
            try {
                Fischer_Burmeister_m(Phi_red(active), Y_crit(active), B_red(active, active), Dl_fb, dDl_a, err_fb_unused);
            } catch (const std::exception &) {
                return r;   // LAPACK failure on the reduced system: non-convergence, the caller step-cuts
            }
            dDl(active) = dDl_a;
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

CuttingPlaneResult cutting_plane_return_mapping(
    vec &Phi,
    vec &Y_crit,
    vec &Ds_total,
    const CuttingPlaneHooks &hooks,
    const ReturnMappingControl &control,
    int iter0) {

    const int maxiter = (control.maxiter > 0) ? control.maxiter : simcoon::maxiter_umat;
    const double precision = (control.precision > 0.) ? control.precision : simcoon::precision_umat;
    const uword N = Phi.n_elem;

    CuttingPlaneResult r;
    mat B = zeros(N, N);
    vec ds = zeros(N);
    double error = 1.0;
    int iter = iter0;
    while (iter < maxiter && error > precision) {
        hooks.jacobian(B);
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_total, ds, error);
        hooks.update(ds);
        hooks.refresh();
        if (hooks.consistency) error += hooks.consistency();
        ++iter;
        if (iter < maxiter && error > precision) hooks.constraints(Phi, Y_crit);
    }
    r.converged = error <= precision;
    r.niter = iter;
    r.error = error;
    return r;
}

void finite_difference_total_derivatives(
    std::vector<ReturnMechanism> &mechanisms,
    ReturnStateHooks &hooks,
    const std::function<bool(const vec &, const vec &)> &refresh,
    const std::function<void()> &save,
    const std::function<void()> &restore,
    double h_sigma_rel,
    double h_Dlambda) {

    const int N = int(mechanisms.size());
    // The multipliers of the current iterate: the stress derivatives are taken at fixed Dlambda,
    // which the mechanism callbacks do not receive.
    auto Dl_cur = std::make_shared<vec>(zeros(N));
    hooks.update_state = [refresh, Dl_cur](const vec &sig, const vec &Dl) {
        *Dl_cur = Dl;
        return refresh(sig, Dl);
    };
    // Evaluate a quantity at (sig, Dl) with the state re-solved there, then put the state back.
    auto probe = [=, &mechanisms](const vec &sig, const vec &Dl, const std::function<vec()> &quantity) {
        save();
        refresh(sig, Dl);
        const vec out = quantity();
        restore();
        return out;
    };
    for (int j = 0; j < N; j++) {
        mechanisms[j].dPhi_dsigma = [=, &mechanisms](const vec &sig) {
            const double h = h_sigma_rel*(norm(sig, 2) + 1.);
            vec g(6);
            for (int c = 0; c < 6; c++) {
                vec sp = sig, sm = sig;
                sp(c) += h;
                sm(c) -= h;
                g(c) = (probe(sp, *Dl_cur, [&]{ return vec{mechanisms[j].Phi(sp)}; })(0)
                      - probe(sm, *Dl_cur, [&]{ return vec{mechanisms[j].Phi(sm)}; })(0))/(2.*h);
            }
            return g;
        };
        mechanisms[j].dLambda_dsigma = [=, &mechanisms](const vec &sig) {
            const double h = h_sigma_rel*(norm(sig, 2) + 1.);
            mat D(6, 6);
            for (int c = 0; c < 6; c++) {
                vec sp = sig, sm = sig;
                sp(c) += h;
                sm(c) -= h;
                D.col(c) = (probe(sp, *Dl_cur, [&]{ return mechanisms[j].Lambda(sp); })
                          - probe(sm, *Dl_cur, [&]{ return mechanisms[j].Lambda(sm); }))/(2.*h);
            }
            return D;
        };
    }
    // Multiplier derivatives: one-sided at Dlambda = 0 (the state is undefined for Dlambda < 0).
    auto bracket = [=](const vec &Dl, int j, vec &Dp, vec &Dm, double &den) {
        Dp = Dl; Dm = Dl;
        Dp(j) += h_Dlambda;
        Dm(j) = std::max(Dl(j) - h_Dlambda, 0.);
        den = Dp(j) - Dm(j);
    };
    hooks.K = [=, &mechanisms](const vec &sig, const vec &Dl) {
        mat K(N, N);
        for (int j = 0; j < N; j++) {
            vec Dp, Dm; double den;
            bracket(Dl, j, Dp, Dm, den);
            for (int l = 0; l < N; l++) {
                K(l, j) = (probe(sig, Dp, [&]{ return vec{mechanisms[l].Phi(sig)}; })(0)
                         - probe(sig, Dm, [&]{ return vec{mechanisms[l].Phi(sig)}; })(0))/den;
            }
        }
        return K;
    };
    hooks.flow_state_coupling = [=, &mechanisms](const vec &sig, const vec &Dl) {
        // Dl_j dLambda^j/dDl_j = d(flux)/dDl_j - Lambda^j, flux = sum_k Dl_k Lambda^k
        auto flux = [&](const vec &D) {
            vec f = zeros(6);
            for (int k = 0; k < N; k++) f += D(k)*mechanisms[k].Lambda(sig);
            return f;
        };
        std::vector<vec> c(N);
        for (int j = 0; j < N; j++) {
            vec Dp, Dm; double den;
            bracket(Dl, j, Dp, Dm, den);
            const vec dflux = (probe(sig, Dp, [&]{ return flux(Dp); }) - probe(sig, Dm, [&]{ return flux(Dm); }))/den;
            c[j] = dflux - mechanisms[j].Lambda(sig);
        }
        return c;
    };
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
