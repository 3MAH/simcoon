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

///@file chaboche_cpp.cpp
///@brief Shared Chaboche-family CPP local solve. Doc in the .hpp.

#include <cmath>
#include <algorithm>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/chaboche_cpp.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

ReturnMappingResult chaboche_cpp_return_mapping(
    const vec &sigma_tr,
    const mat &L,
    const std::function<double(const vec &)> &crit,
    const std::function<vec(const vec &)> &flow,
    const double &sigmaY, const double &Q, const double &b,
    const double &C_1, const double &D_1, const double &C_2, const double &D_2,
    const double &p_n, const vec &a_1n, const vec &a_2n,
    const vec &X_1n, const vec &X_2n, const double &Hp_n,
    double &p, double &Hp, double &dHpdp,
    vec &a_1, vec &a_2, vec &X_1, vec &X_2, vec &X) {

    // Backward-Euler state at (sig, Dl): fixed point on X; updates the referenced state.
    auto solve_state = [&](const vec &sig, double Dl) {
        for (int it_fp = 0; it_fp < 20; it_fp++) {
            const vec eta_xi = flow(sig - X);
            const vec a1l = (a_1n + Dl*eta_xi)/(1. + Dl*D_1);
            const vec a2l = (a_2n + Dl*eta_xi)/(1. + Dl*D_2);
            const vec X1l = X_1n + (2./3.)*C_1*((a1l - a_1n) % Ir05());
            const vec X2l = X_2n + (2./3.)*C_2*((a2l - a_2n) % Ir05());
            const vec Xn = X1l + X2l;
            const double dX = norm(Xn - X, 2);
            a_1 = a1l; a_2 = a2l; X_1 = X1l; X_2 = X2l; X = Xn;
            if (dX < 1e-12*(norm(Xn, 2) + 1.)) break;
        }
        Hp = (Hp_n + b*Q*Dl)/(1. + b*Dl);
        dHpdp = b*(Q - Hp);
    };
    // Guarded probe of Phi at (sig, Dl) with the inner state solved, captures restored.
    auto probe_Phi = [&](const vec &sig, double Dl) {
        const vec Xs1 = X_1, Xs2 = X_2, Xs = X, as1 = a_1, as2 = a_2;
        const double Hs = Hp, dHs = dHpdp;
        solve_state(sig, Dl);
        const double v = crit(sig - X) - Hp - sigmaY;
        X_1 = Xs1; X_2 = Xs2; X = Xs; a_1 = as1; a_2 = as2;
        Hp = Hs; dHpdp = dHs;
        return v;
    };
    // With the inner-consistent state map (X, Hp) = V_hat(sigma, Dl), ALL derivative
    // callbacks must be TOTAL derivatives of the composed map (partial-only inputs
    // leave an O(dX/dsigma) error in the consistent tangent).
    ReturnMechanism mech;
    mech.Phi            = [&](const vec &sig) { return crit(sig - X) - Hp - sigmaY; };
    mech.dPhi_dsigma    = [&](const vec &sig) {
        const double hfd = 1.e-5*(norm(sig, 2) + 1.);
        const double Dl = p - p_n;
        vec g(6);
        for (int c6 = 0; c6 < 6; c6++) {
            vec sp = sig, sm = sig;
            sp(c6) += hfd;
            sm(c6) -= hfd;
            g(c6) = (probe_Phi(sp, Dl) - probe_Phi(sm, Dl))/(2.*hfd);
        }
        return g;
    };
    mech.Lambda         = [&](const vec &sig) { return flow(sig - X); };
    mech.dLambda_dsigma = [&](const vec &sig) {
        const double hfd = 1.e-5*(norm(sig, 2) + 1.);
        const double Dl = p - p_n;
        const vec Xsave1 = X_1, Xsave2 = X_2, Xsave = X, asave1 = a_1, asave2 = a_2;
        const double Hpsave = Hp, dHsave = dHpdp;
        mat D(6, 6);
        for (int c6 = 0; c6 < 6; c6++) {
            vec sp = sig, sm = sig;
            sp(c6) += hfd;
            sm(c6) -= hfd;
            solve_state(sp, Dl);
            const vec ep = flow(sp - X);
            solve_state(sm, Dl);
            const vec em = flow(sm - X);
            D.col(c6) = (ep - em)/(2.*hfd);
        }
        X_1 = Xsave1; X_2 = Xsave2; X = Xsave; a_1 = asave1; a_2 = asave2;
        Hp = Hpsave; dHpdp = dHsave;
        return D;
    };
    return closest_point_return_mapping(sigma_tr, L, mech,
             [&](const vec &sig, double Dl) {
                 p = p_n + Dl;
                 solve_state(sig, Dl);
                 return true;
             },
             // TOTAL dPhi/dDlambda at fixed sigma (Voce + both backstress chains).
             [&](const vec &sig, double Dl) {
                 const double hDl = 1.e-6*(fabs(Dl) + 1.e-8);
                 const double Dlp = Dl + hDl;
                 const double Dlm = std::max(Dl - hDl, 0.);
                 return (probe_Phi(sig, Dlp) - probe_Phi(sig, Dlm))/(Dlp - Dlm);
             },
             sigmaY,
             // multiplier-side state chain c = Dl*L*dLambda/dDl (FD of the consistent map):
             // required for the exact consistent tangent with backstress (doc eq:Lambda_tilde_state).
             [&](const vec &sig, double Dl) -> vec {
                 const vec Xs1 = X_1, Xs2 = X_2, Xs = X, as1 = a_1, as2 = a_2;
                 const double Hs = Hp, dHs = dHpdp;
                 const double hDl = 1.e-6*(fabs(Dl) + 1.e-8);
                 const double Dlp = Dl + hDl;
                 const double Dlm = std::max(Dl - hDl, 0.);
                 solve_state(sig, Dlp);
                 const vec ep = flow(sig - X);
                 solve_state(sig, Dlm);
                 const vec em = flow(sig - X);
                 X_1 = Xs1; X_2 = Xs2; X = Xs; a_1 = as1; a_2 = as2;
                 Hp = Hs; dHpdp = dHs;
                 return vec(Dl*(L*((ep - em)/(Dlp - Dlm))));
             });
}

} // namespace simcoon
