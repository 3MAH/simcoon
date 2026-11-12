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
///@brief Return-mapping integrators (closest point, cutting plane). Contracts and algorithms in the .hpp.

#include <cmath>
#include <stdexcept>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Simulation/Maths/small_lu.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

namespace {
int resolved_maxiter(const ReturnMappingControl &c) { return (c.maxiter > 0) ? c.maxiter : simcoon::maxiter_umat; }
double resolved_precision(const ReturnMappingControl &c) { return (c.precision > 0.) ? c.precision : simcoon::precision_umat; }

// per-row buffers of the Newton, sized once per call: the callbacks write into them
struct RowScratch {
    vec::fixed<6> n;          // criterion gradient n^l
    vec::fixed<6> Lambda;     // flow direction
    vec::fixed<6> kappa;      // L : Lambda
    vec::fixed<6> c;          // L : (multiplier-side state chain)
    vec::fixed<6> Minv_kc;    // M^-1 (kappa + c)
    mat::fixed<6,6> D;        // dLambda/dsigma
};

// M += a A B for 6x6 operands, in place: no 6x6 temporary on the heap
void add_scaled_product(mat::fixed<6,6> &M, const double a, const mat &A, const mat &B) {
    for (int k = 0; k < 6; k++) {
        for (int m = 0; m < 6; m++) {
            const double b = a*B(m, k);
            if (b == 0.) continue;
            for (int i = 0; i < 6; i++) M(i, k) += A(i, m)*b;
        }
    }
}
}

ReturnMappingResult closest_point_return_mapping(
    const vec &sigma_tr,
    const mat &L,
    const std::vector<ReturnMechanism> &mechanisms,
    const ReturnStateHooks &hooks,
    const vec &Y_crit,
    const ReturnMappingControl &control) {

    const int N = int(mechanisms.size());
    const int maxiter = resolved_maxiter(control);
    const double precision = resolved_precision(control);

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
    std::vector<RowScratch> row(N);
    mat K(N, N), dLdDs(6, N);
    vec::fixed<6> R_sigma;

    // Elastic block: Lc is its tangent at the iterate (the given L for a linear block).
    mat Lc = L;
    const bool nonlinear_block = static_cast<bool>(hooks.elastic_response);
    vec::fixed<6> eps_el;
    vec sigma_F;
    // criterion gradients only at the Newton top (eval_gradients): with FD callbacks each is a full state re-solve
    auto eval_basic = [&](const vec &sig, const vec &Dl) {
        for (int j = 0; j < N; j++) {
            Phi(j) = mechanisms[j].Phi(sig);
            if (mechanisms[j].Lambda) mechanisms[j].Lambda(sig, row[j].Lambda);
            else mechanisms[j].dPhi_dsigma(sig, row[j].Lambda);
        }
        if (nonlinear_block) {
            eps_el = hooks.eps_el_tr;
            for (int j = 0; j < N; j++) eps_el -= Dl(j) * row[j].Lambda;
            hooks.elastic_response(eps_el, sigma_F, Lc);
            R_sigma = sig - sigma_F;
            for (int j = 0; j < N; j++) row[j].kappa = Lc * row[j].Lambda;
        } else {
            R_sigma = sig - sigma_tr;
            for (int j = 0; j < N; j++) {
                row[j].kappa = L * row[j].Lambda;
                R_sigma += Dl(j) * row[j].kappa;
            }
        }
    };

    auto eval_gradients = [&](const vec &sig) {
        for (int j = 0; j < N; j++) {
            if (mechanisms[j].Lambda) mechanisms[j].dPhi_dsigma(sig, row[j].n);
            else row[j].n = row[j].Lambda;
        }
    };

    if (!refresh(r.sigma, r.Dlambda)) return r;
    eval_basic(r.sigma, r.Dlambda);

    // Tangent ingredients at the converged iterate from the pieces the Newton just built
    // (kappa_eff = kappa + c: the multiplier-side state chain enters through the flux).
    auto fill_tangent_pieces = [&]() {
        r.Bhat_continuum.set_size(N, N);
        r.kappa_j.set_size(6, N);
        r.dPhidsigma_l.set_size(6, N);
        r.dLambda_dsigma_l.set_size(6, 6 * N);
        vec::fixed<6> kappa_c;
        for (int j = 0; j < N; j++) {
            kappa_c = row[j].kappa + row[j].c;
            r.kappa_j.col(j) = kappa_c;
            r.dPhidsigma_l.col(j) = row[j].n;
            r.dLambda_dsigma_l.cols(6 * j, 6 * j + 5) = row[j].D;
            for (int l = 0; l < N; l++) {
                r.Bhat_continuum(l, j) = sum(row[l].n % kappa_c) - K(l, j);
            }
        }
    };

    // admissible trial: return it; Dlambda = 0 masks every mechanism, so no state callbacks
    if (Phi.max() <= 0.) {
        for (int j = 0; j < N; j++) {   // masked by Dlambda = 0 in the assembly: the gradients are not evaluated
            row[j].n.zeros();
            row[j].D.zeros();
            row[j].c.zeros();
        }
        K.zeros();
        fill_tangent_pieces();
        r.converged = true;
        return r;
    }

    const int snap_niter = 20;   // plain Newton gets this many iterations before the snap may act
    r.error_history.reserve(8);  // a converging solve takes a handful of iterations

    // M factorized on the stack (small_lu)
    mat::fixed<6,6> M;
    small_lu<6> lu_M;
    vec Phi_red(N), dDl(N), Dl_old(N), Dl_fb, dDl_a;
    mat B_red(N, N);
    vec::fixed<6> MinvR, dsigma, sigma_old;
    uvec active(N);
    ncp_step local_step(control.ncp);   // Fischer-Burmeister, or the min step with its fallback: the same root
    for (r.niter = 0; r.niter < maxiter; r.niter++) {

        // Newton ingredients at the current (state-consistent) iterate.
        eval_gradients(r.sigma);
        if (hooks.K) hooks.K(r.sigma, r.Dlambda, K);
        else K.zeros();
        M.eye();
        for (int j = 0; j < N; j++) {
            if (mechanisms[j].dLambda_dsigma) mechanisms[j].dLambda_dsigma(r.sigma, row[j].D);
            else row[j].D.zeros();
            add_scaled_product(M, r.Dlambda(j), Lc, row[j].D);
        }
        if (!lu_M.factor(M)) return r;   // condensation breakdown: caller step-cuts
        if (hooks.flow_state_coupling) {
            hooks.flow_state_coupling(r.sigma, r.Dlambda, dLdDs);
            for (int j = 0; j < N; j++) row[j].c = Lc * dLdDs.col(j);
        } else {
            for (int j = 0; j < N; j++) row[j].c.zeros();
        }
        MinvR = R_sigma;
        lu_M.solve_in_place(MinvR.memptr());
        for (int j = 0; j < N; j++) {
            row[j].Minv_kc = row[j].kappa + row[j].c;
            lu_M.solve_in_place(row[j].Minv_kc.memptr());
        }
        for (int l = 0; l < N; l++) {
            Phi_red(l) = Phi(l) - sum(row[l].n % MinvR);
            for (int j = 0; j < N; j++) {
                B_red(l, j) = -sum(row[l].n % row[j].Minv_kc) + K(l, j);
            }
        }

        // Merit on the TRUE residuals (Phi, R_sigma) at the current iterate.
        if (!B_red.is_finite()) return r;
        const double err = Fischer_Burmeister_residual(Phi, Y_crit, B_red, r.Dlambda) + norm(R_sigma, 2) / sigma_ref;
        r.error_history.push_back(err);
        r.error = err;
        if (err < precision) {
            fill_tangent_pieces();
            r.converged = true;
            return r;
        }

        // active set excludes rows with Phi < 0, no multiplier and a vanishing diagonal (zero FB row, singular solve)
        uword n_active = 0;
        for (int l = 0; l < N; l++) {
            if (Phi(l) >= 0. || r.Dlambda(l) > 0. || std::abs(B_red(l, l)) >= simcoon::iota) active(n_active++) = l;
        }
        dDl.zeros();
        if (n_active > 0) {
            double err_fb_unused = 0.;
            try {
                if (int(n_active) == N) {   // every row active: no gather
                    Dl_fb = r.Dlambda;
                    local_step(Phi_red, Y_crit, B_red, Dl_fb, dDl, err_fb_unused);
                } else {
                    const uvec act = active.head(n_active);
                    Dl_fb = r.Dlambda(act);
                    dDl_a.set_size(n_active);
                    local_step(Phi_red(act), Y_crit(act), B_red(act, act), Dl_fb, dDl_a, err_fb_unused);
                    dDl(act) = dDl_a;
                }
            } catch (const std::exception &) {
                return r;   // LAPACK failure on the reduced system: non-convergence, the caller step-cuts
            }
        }

        dsigma = -MinvR;
        for (int j = 0; j < N; j++) dsigma -= dDl(j) * row[j].Minv_kc;

        // Backtracking on the merit: semi-smooth iterations are locally non-monotone, so the
        // step is halved only when the merit grows by more than x2.
        sigma_old = r.sigma;
        Dl_old = r.Dlambda;
        double scale = 1.;
        int bt = 0;
        while (true) {
            r.sigma = sigma_old + scale * dsigma;
            for (int l = 0; l < N; l++) r.Dlambda(l) = std::max(Dl_old(l) + scale * dDl(l), 0.);   // clamp at 0
            // Deactivation snap (see the header): exact inactivity for a row hovering at the
            // complementarity boundary, gated on a stalled iteration.
            if (r.niter >= snap_niter) {
                for (int l = 0; l < N; l++) {
                    if ((Phi(l) < 0.) && (r.Dlambda(l) * fabs(B_red(l, l)) < 1.e-2 * Y_crit(l))) {
                        r.Dlambda(l) = 0.;
                    }
                }
            }
            // failed inner state solve or non-finite trial = rejected step: halve it
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
    hooks.K = [&K_scalar](const vec &sig, const vec &Dl, mat &K) {
        K(0, 0) = K_scalar ? K_scalar(sig, Dl(0)) : 0.;
    };
    if (flow_state_coupling) {
        hooks.flow_state_coupling = [&flow_state_coupling](const vec &sig, const vec &Dl, mat &dLdDs) {
            dLdDs.col(0) = flow_state_coupling(sig, Dl(0));
        };
    }
    vec Y(1);
    Y(0) = Y_crit;
    return closest_point_return_mapping(sigma_tr, L, mechs, hooks, Y, control);
}

bool cutting_plane_return_mapping(
    vec &Phi,
    vec &Y_crit,
    vec &Ds_total,
    const CuttingPlaneHooks &hooks,
    const ReturnMappingControl &control,
    int iter0) {

    const int maxiter = resolved_maxiter(control);
    const double precision = resolved_precision(control);
    const uword N = Phi.n_elem;

    mat B = zeros(N, N);
    vec ds = zeros(N);
    double error = 1.0;
    int iter = iter0;
    ncp_step local_step(control.ncp);
    while (iter < maxiter && error > precision) {
        hooks.jacobian(B);
        local_step(Phi, Y_crit, B, Ds_total, ds, error);
        hooks.update(ds);
        hooks.refresh();
        if (hooks.consistency) error += hooks.consistency();
        ++iter;
        if (iter < maxiter && error > precision) hooks.constraints(Phi, Y_crit);
    }
    return error <= precision;
}

ContinuumTangent cpp_consistent_tangent(const ReturnMappingResult &r, const mat &L) {
    if (!r.converged || r.kappa_j.is_empty()) {
        // Precondition not met: elastic operator as a sized placeholder (thermomechanical
        // callers read P_epsilon[l] unconditionally on the step-cut path).
        ContinuumTangent ct;
        ct.Lt = L;
        ct.L_tilde = L;
        const uword Nm = std::max(r.Dlambda.n_elem, uword(1));
        ct.invBhat = zeros(Nm, Nm);
        ct.P_epsilon.zeros(6, Nm);
        return ct;
    }
    return assemble_algorithmic_tangent(r.Bhat_continuum, r.kappa_j, r.dPhidsigma_l,
                                        r.Dlambda, L, r.dLambda_dsigma_l);
}

} // namespace simcoon
