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

#include <armadillo>
#include <string>
#include <stdexcept>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>

using namespace arma;

namespace simcoon {

int default_nshr(const int &ndi) {
    return (ndi == 1) ? 0 : ((ndi == 2) ? 1 : 3);
}

condensed_state make_condensed_state(const int &ndi, const int &nshr) {

    const int ns = (nshr < 0) ? default_nshr(ndi) : nshr;
    if (ndi < 1 || ndi > 3 || ns > 3) {
        throw std::invalid_argument("invalid stress state (ndi = " + std::to_string(ndi) + ", nshr = " + std::to_string(nshr)
                                    + "): ndi is 1, 2 or 3 and nshr 0 to 3");
    }
    condensed_state cs;
    cs.K.set_size(ndi + ns);
    int nk = 0;
    for (int i = 0; i < ndi; i++) cs.K(nk++) = i;
    for (int i = 0; i < ns; i++) cs.K(nk++) = 3 + i;
    if (ndi < 3) {
        cs.C.set_size(6 - ndi - ns);
        int nc = 0;
        for (int i = 0; i < 6; i++) {
            if (!any(cs.K == static_cast<uword>(i))) cs.C(nc++) = i;
        }
    }
    return cs;
}

bool condensed_newton_step(const condensed_state &cs, vec &DEtot, const vec &sigma, const mat &Lt, bool &converged) {

    const vec res = sigma(cs.C);
    if (!res.is_finite()) return false;
    const mat L_CC = Lt(cs.C, cs.C);
    // relative to the stress, with a floor at the stress of a 1e-14 strain (zero load)
    converged = abs(res).max() <= 1.e-10*abs(sigma).max() + 1.e-14*abs(L_CC).max();
    if (converged) return true;
    vec dE;
    if (!solve(dE, L_CC, res, solve_opts::no_approx) || !dE.is_finite()) return false;
    DEtot(cs.C) -= dE;
    return true;
}

void condense_tangent(const condensed_state &cs, mat &Lt) {

    if (!cs.active()) return;
    const mat Q = Lt(cs.K, cs.K) - Lt(cs.K, cs.C)*solve(Lt(cs.C, cs.C), Lt(cs.C, cs.K));
    Lt.zeros();
    Lt(cs.K, cs.K) = Q;
}

void condense_tangent_T(const condensed_state &cs, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT) {

    if (!cs.active()) return;
    // dSdT and drdE are 6-entry vectors stored 6x1 or 1x6: work on their entries, keep the shape
    const vec qT = vectorise(dSdT), rE = vectorise(drdE);
    const mat X = solve(dSdE(cs.C, cs.C), join_rows(dSdE(cs.C, cs.K), qT(cs.C)));   // [L_CC^-1 L_CK | L_CC^-1 dSdT_C]
    const uword nK = cs.K.n_elem;
    const mat X_E = X.cols(0, nK - 1);
    const vec X_T = X.col(nK);
    const mat Q = dSdE(cs.K, cs.K) - dSdE(cs.K, cs.C)*X_E;
    const vec qT_K = qT(cs.K) - dSdE(cs.K, cs.C)*X_T;
    const vec rE_K = rE(cs.K) - (rE(cs.C).t()*X_E).t();
    const double rT = drdT(0, 0) - dot(rE(cs.C), X_T);
    dSdE.zeros();
    dSdE(cs.K, cs.K) = Q;
    dSdT.zeros();
    drdE.zeros();
    for (uword k = 0; k < nK; k++) {
        dSdT(cs.K(k)) = qT_K(k);
        drdE(cs.K(k)) = rE_K(k);
    }
    drdT(0, 0) = rT;
}

} //namespace simcoon
