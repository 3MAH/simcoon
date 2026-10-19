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

///@file damage_LLD_0.cpp
///@brief Ladeveze-Le Dantec damage model with coupled transverse-shear plasticity in effective
///       stress (LLDM0); the equations and the integration scheme are in damage_LLD_0.hpp

#include <iostream>
#include <math.h>
#include <string>
#include <array>
#include <vector>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

namespace {

// d22 at transverse failure (Y_t > Y_22_u): held there, d12 keeps its own law
constexpr double d_broken = 0.99;

// The fibre frame (fibres along 1), damaged: the compliance of the 22-33 block and of G23 scaled
// by 1/(1 - d22), of G12 and G13 by 1/(1 - d12); nuTL is the minor Poisson ratio, as in L_isotrans.
// The stiffness is the inverse of the 11-22-33 block plus the reciprocal shears.
struct damaged_elasticity {
    double EL, ET, nuTL, nuTT, GLT;

    mat L(const double &d22, const double &d12) const {
        const double ET_d = ET*(1. - d22);
        mat::fixed<3, 3> S;
        S(0, 0) = 1./EL;
        S(0, 1) = S(1, 0) = S(0, 2) = S(2, 0) = -nuTL/ET;
        S(1, 1) = S(2, 2) = 1./ET_d;
        S(1, 2) = S(2, 1) = -nuTT/ET_d;
        mat L = zeros(6, 6);
        L.submat(0, 0, 2, 2) = inv(S);
        L(3, 3) = L(4, 4) = GLT*(1. - d12);
        L(5, 5) = ET_d/(2.*(1. + nuTT));
        return L;
    }
    // Lambda_j = dS/dd_j : sigma, the strain the damage releases at fixed stress
    vec::fixed<6> Lambda_22(const vec &sigma, const double &d22) const {
        const double c = 1./(ET*(1. - d22)*(1. - d22));
        vec::fixed<6> v = arma::fill::zeros;
        v(1) = c*(sigma(1) - nuTT*sigma(2));
        v(2) = c*(sigma(2) - nuTT*sigma(1));
        v(5) = c*2.*(1. + nuTT)*sigma(5);
        return v;
    }
    vec::fixed<6> Lambda_12(const vec &sigma, const double &d12) const {
        const double c = 1./(GLT*(1. - d12)*(1. - d12));
        vec::fixed<6> v = arma::fill::zeros;
        v(3) = c*sigma(3);
        v(4) = c*sigma(4);
        return v;
    }
};

// The rotation that brings the fibre axis (1, 2 or 3) onto 1: local = R^T global R
mat fibre_frame(const int &axis) {
    switch (axis) {
        case 1: return eye(3, 3);
        case 2: return mat({{0., -1., 0.}, {1., 0., 0.}, {0., 0., 1.}});   // e1 -> e2
        case 3: return mat({{0., 0., -1.}, {0., 1., 0.}, {1., 0., 0.}});   // e1 -> e3
        default:
            throw std::invalid_argument("LLDM0: props(0) is the fibre axis, 1, 2 or 3; got " + std::to_string(axis));
    }
}

// The damage driving force of the criterion (unilateral in 22 and 33) and its derivatives
struct driving_force {
    double Y22 = 0., Y33 = 0., Y12 = 0., Y13 = 0., Y_t = 0., Y_ts = 0.;
    vec::fixed<6> dYts_dsigma = arma::fill::zeros;
    double dYts_dd22 = 0., dYts_dd12 = 0.;

    driving_force(const vec &sigma, const double &d22, const double &d12, const double &ET, const double &GLT, const double &b) {
        const double s22 = Macaulay_p(sigma(1)), s33 = Macaulay_p(sigma(2));
        const double c22 = ET*(1. - d22)*(1. - d22), c12 = GLT*(1. - d12)*(1. - d12);
        Y22 = 0.5*s22*s22/c22;
        Y33 = 0.5*s33*s33/c22;
        Y12 = 0.5*sigma(3)*sigma(3)/c12;
        Y13 = 0.5*sigma(4)*sigma(4)/c12;
        Y_t = sqrt(Y22 + Y33);
        Y_ts = sqrt(Y12 + Y13 + b*(Y22 + Y33));
        if (Y_ts > simcoon::limit) {
            const double h = 0.5/Y_ts;
            dYts_dsigma(1) = h*b*s22/c22;
            dYts_dsigma(2) = h*b*s33/c22;
            dYts_dsigma(3) = h*sigma(3)/c12;
            dYts_dsigma(4) = h*sigma(4)/c12;
            dYts_dd22 = h*b*2.*(Y22 + Y33)/(1. - d22);
            dYts_dd12 = h*2.*(Y12 + Y13)/(1. - d12);
        }
    }
};

}  // namespace

void umat_damage_LLD_0(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode) {

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(ndi);
    UNUSED(nshr);
    UNUSED(tangent_mode);   // continuum tangent (see the header)

    double Tinit = statev(0);
    double d_22 = statev(1);
    double d_12 = statev(2);
    double p_ts = statev(3);
    vec EP = statev.subvec(4, 9);

    // Rotation of the tensorial internal variables, then everything in the fibre frame
    EP = rotate_strain(EP, DR);
    const int axis = static_cast<int>(props(0));
    const mat R_f = fibre_frame(axis);
    const bool rotated = (axis != 1);
    const vec Etot_f = rotated ? rotate_strain(Etot, R_f.t()) : Etot;
    const vec DEtot_f = rotated ? rotate_strain(DEtot, R_f.t()) : DEtot;
    vec sigma = rotated ? rotate_stress(stress, R_f.t()) : stress;
    if (rotated)
        EP = rotate_strain(EP, R_f.t());

    const damaged_elasticity C{props(1), props(2), props(3), props(4), props(5)};
    const double alphaL = props(6);
    const double alphaT = props(7);
    const double Y_12_0 = props(8);
    const double Y_12_c = props(9);
    const double Y_22_0 = props(10);
    const double Y_22_c = props(11);
    const double Y_22_u = props(12);
    const double b = props(13);
    const double A_ts = props(14);
    const double sigma_ts_0 = props(15);
    const double m_ts = props(16);      // hardening R(p) = beta p^m
    const double beta_ts = props(17);

    // Lagrange penalty keeping d below 1 (lagrange_pow_1)
    const double c_lambda = 1.E-2;
    const double p0_lambda = 1.E-1;
    const double n_lambda = 1.0;
    const double alpha_lambda = 0.;

    if (start) {
        Tinit = T;
        d_22 = 0.;
        d_12 = 0.;
        p_ts = 10.*simcoon::limit;
        EP = zeros(6);
        sigma = zeros(6);
    }

    const mat L0 = C.L(0., 0.);   // undamaged: the effective-stress plasticity
    vec alpha = alphaT*Ith();
    alpha(0) += alphaL - alphaT;

    auto hardening = [&](const double &p) { return (p > simcoon::limit) ? beta_ts*pow(p, m_ts) : simcoon::iota; };
    auto dhardening = [&](const double &p) { return (p > simcoon::limit) ? m_ts*beta_ts*pow(p, m_ts - 1.) : simcoon::iota; };

    // Start of the increment, for the work split
    const vec sigma_start = sigma;
    const vec EP_start = EP;
    const double Hp_start = hardening(p_ts);
    const vec d_start = {d_22, d_12};
    const vec Yth_start = {0.5*sum(sigma_start % C.Lambda_22(sigma_start, d_22)), 0.5*sum(sigma_start % C.Lambda_12(sigma_start, d_12))};
    bool broken_22 = (d_22 >= d_broken);

    const vec theta = {0., A_ts, A_ts, 1., 1., 0.};   // the diagonal of Theta

    // ---- 1. Plasticity in effective stress (independent of the damage) ----
    vec Eel = Etot_f + DEtot_f - alpha*(T + DT - Tinit) - EP;
    vec Lambda_p = zeros(6);
    vec Phi_p = zeros(1), Dp = zeros(1), dp = zeros(1);
    const vec Y_pcrit = {sigma_ts_0};
    mat denom_p = zeros(1, 1);
    double error = 1.;
    int compteur = 0;
    for (compteur = 0; (compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat); compteur++) {
        const vec sigma_eff_ts = theta % (L0*Eel);
        Phi_p(0) = Mises_stress(sigma_eff_ts) - hardening(p_ts) - sigma_ts_0;
        Lambda_p = theta % eta_stress(sigma_eff_ts);   // associated: the gradient of the criterion in effective stress
        denom_p(0, 0) = -sum(Lambda_p % (L0*Lambda_p)) - dhardening(p_ts);
        Fischer_Burmeister_m(Phi_p, Y_pcrit, denom_p, Dp, dp, error);
        p_ts += dp(0);
        EP += dp(0)*Lambda_p;
        Eel = Etot_f + DEtot_f - alpha*(T + DT - Tinit) - EP;
    }
    if (compteur == simcoon::maxiter_umat)
        tnew_dt = 0.2;
    Lambda_p = theta % eta_stress(vec(theta % (L0*Eel)));

    // ---- 2. Damage at fixed plastic strain: d22 and d12 driven by the same Y_ts ----
    const double Yc[2] = {Y_22_c, Y_12_c};
    const double Y0[2] = {Y_22_0, Y_12_0};
    vec d = {d_22, d_12};
    vec Dd_all = zeros(2);
    mat L_t = C.L(d(0), d(1));
    sigma = L_t*Eel;
    // dPhi_l/dd_j at fixed strain: through sigma = L~(d) Eel and explicitly through Y_ts
    auto dPhi_dd = [&](const driving_force &Y, const std::array<vec::fixed<6>, 2> &Lambda_d, const int &l, const int &j) {
        const double dYts_dd = (j == 0) ? Y.dYts_dd22 : Y.dYts_dd12;
        double J = ((Y.Y_ts > Y0[l]) ? (-sum(Y.dYts_dsigma % (L_t*Lambda_d[j])) + dYts_dd)/Yc[l] : 0.);
        if (l == j)
            J -= dlagrange_pow_1(d(l), c_lambda, p0_lambda, n_lambda, alpha_lambda) + 1.;
        return J;
    };
    error = 1.;
    for (compteur = 0; (compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat); compteur++) {
        const driving_force Y(sigma, d(0), d(1), C.ET, C.GLT, b);
        if (!broken_22 && Y.Y_t > Y_22_u) {   // transverse failure
            Dd_all(0) += d_broken - d(0);
            d(0) = d_broken;
            broken_22 = true;
            L_t = C.L(d(0), d(1));
            sigma = L_t*Eel;
            continue;
        }
        const std::vector<int> mech = broken_22 ? std::vector<int>{1} : std::vector<int>{0, 1};
        const int n = mech.size();
        const std::array<vec::fixed<6>, 2> Lambda_d = {C.Lambda_22(sigma, d(0)), C.Lambda_12(sigma, d(1))};
        vec Phi_d(n), Y_dcrit = ones(n), Dd(n), dd(n);
        mat denom_d(n, n);
        for (int a = 0; a < n; a++) {
            const int l = mech[a];
            Phi_d(a) = Macaulay_p(Y.Y_ts - Y0[l])/Yc[l] - lagrange_pow_1(d(l), c_lambda, p0_lambda, n_lambda, alpha_lambda) - d(l);
            Dd(a) = Dd_all(l);
            for (int c = 0; c < n; c++)
                denom_d(a, c) = dPhi_dd(Y, Lambda_d, l, mech[c]);
        }
        Fischer_Burmeister_m(Phi_d, Y_dcrit, denom_d, Dd, dd, error);
        for (int a = 0; a < n; a++) {
            d(mech[a]) += dd(a);
            Dd_all(mech[a]) = Dd(a);
        }
        L_t = C.L(d(0), d(1));
        sigma = L_t*Eel;
    }
    if (compteur == simcoon::maxiter_umat)
        tnew_dt = 0.2;
    d_22 = d(0);
    d_12 = d(1);

    // ---- 3. Tangent of the two-stage update (plasticity, then damage) ----
    // dEel/deps = I - Lambda_p (x) L0 Lambda_p / h_p when plasticity is active (the direction held
    // fixed: continuum); the active damage mechanisms then solve dPhi_l = 0 with
    // dsigma = L~ dEel - L~ Lambda_j dd_j
    mat A = eye(6, 6);
    if (Dp(0) > simcoon::iota) {
        const double h_p = sum(Lambda_p % (L0*Lambda_p)) + dhardening(p_ts);
        A -= Lambda_p*(L0*Lambda_p).t()/h_p;
    }
    Lt = L_t*A;
    {
        const driving_force Y(sigma, d(0), d(1), C.ET, C.GLT, b);
        const std::array<vec::fixed<6>, 2> Lambda_d = {C.Lambda_22(sigma, d(0)), C.Lambda_12(sigma, d(1))};
        std::vector<int> act;
        for (int l = 0; l < 2; l++)
            if (Dd_all(l) > simcoon::iota && !(l == 0 && broken_22) && Y.Y_ts > Y0[l])
                act.push_back(l);
        const int n = act.size();
        if (n > 0) {
            mat J(n, n), G(n, 6);
            for (int a = 0; a < n; a++) {
                for (int c = 0; c < n; c++)
                    J(a, c) = dPhi_dd(Y, Lambda_d, act[a], act[c]);
                G.row(a) = (A.t()*(L_t*Y.dYts_dsigma)).t()/Yc[act[a]];
            }
            const mat dd_deps = -solve(J, G);
            for (int a = 0; a < n; a++)
                Lt -= (L_t*Lambda_d[act[a]])*dd_deps.row(a);
        }
    }

    // ---- 4. Work split: stored elastic energy, stored hardening, dissipation ----
    const vec DEP = EP - EP_start;
    const double Hp = hardening(p_ts);
    const vec Yth = {0.5*sum(sigma % C.Lambda_22(sigma, d_22)), 0.5*sum(sigma % C.Lambda_12(sigma, d_12))};
    const vec Dd_inc = vec({d_22, d_12}) - d_start;
    const double Dgamma_d = 0.5*sum((Yth_start + Yth) % Dd_inc);
    const double Dgamma_p = 0.5*sum((sigma_start + sigma) % DEP) - 0.5*(Hp_start + Hp)*Dp(0);
    Wm += 0.5*sum((sigma_start + sigma) % DEtot_f);
    Wm_r += 0.5*sum((sigma_start + sigma) % (DEtot_f - DEP)) - Dgamma_d;
    Wm_ir += 0.5*(Hp_start + Hp)*Dp(0);
    Wm_d += Dgamma_p + Dgamma_d;

    // Back to the global frame
    if (rotated) {
        stress = rotate_stress(sigma, R_f);
        Lt = rotate_stiffness(Lt, R_f);
        L = rotate_stiffness(L_t, R_f);
        EP = rotate_strain(EP, R_f);
    }
    else {
        stress = sigma;
        L = L_t;
    }

    statev(0) = Tinit;
    statev(1) = d_22;
    statev(2) = d_12;
    statev(3) = p_ts;
    statev.subvec(4, 9) = EP;
}

} //namespace simcoon
