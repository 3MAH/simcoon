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

///@file linear_viscoelastic.cpp
///@brief Closed-form backward-Euler step of the linear viscoelastic rheologies (see the header)

#include <vector>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

namespace {

// B_i^-1 = (I + dt H^-1 L)^-1 and C_i = dt (H + dt L)^-1; dt = 0 gives I and 0 (inactive branch)
void branch_operators(const mat &L, const mat &H, const double &DTime, mat &B_inv, mat &C) {
    const mat HdtL_inv = inv(H + DTime*L);
    B_inv = HdtL_inv*H;
    C = DTime*HdtL_inv;
}

}  // namespace

LinearViscoStep kelvin_series_step(const mat &L0, const vector<mat> &L_i, const vector<mat> &H_i,
                                   const vector<vec> &EV_i_start, const vec &eps_e,
                                   const vec &alpha, const double &DTime) {
    const size_t N = L_i.size();
    vector<mat> C(N);
    vector<vec> EV_relaxed(N);   // B_i^-1 EV_i,n
    mat sumC = zeros(6,6);
    vec rhs = eps_e;
    for (size_t i=0; i<N; i++) {
        mat B_inv;
        branch_operators(L_i[i], H_i[i], DTime, B_inv, C[i]);
        EV_relaxed[i] = B_inv*EV_i_start[i];
        sumC += C[i];
        rhs -= EV_relaxed[i];
    }
    const mat S_inv = inv(eye(6,6) + L0*sumC);

    LinearViscoStep st;
    st.sigma = S_inv*(L0*rhs);
    st.dSdE = S_inv*L0;
    st.dSdE = 0.5*(st.dSdE + st.dSdE.t());   // (L0^-1 + sum C)^-1: symmetric up to round-off
    st.dSdT = -st.dSdE*alpha;
    st.EV_i.resize(N);
    st.dEVdE_i.resize(N);
    st.dEVdT_i.resize(N);
    for (size_t i=0; i<N; i++) {
        st.EV_i[i] = EV_relaxed[i] + C[i]*st.sigma;
        st.dEVdE_i[i] = C[i]*st.dSdE;
        st.dEVdT_i[i] = C[i]*st.dSdT;
    }
    return st;
}

LinearViscoStep maxwell_parallel_step(const mat &L0, const vector<mat> &L_i, const vector<mat> &H_i,
                                      const vector<vec> &EV_i_start, const vec &eps_e,
                                      const vec &alpha, const double &DTime) {
    const size_t N = L_i.size();
    LinearViscoStep st;
    st.sigma = L0*eps_e;
    st.dSdE = L0;
    st.EV_i.resize(N);
    st.dEVdE_i.resize(N);
    st.dEVdT_i.resize(N);
    for (size_t i=0; i<N; i++) {
        mat B_inv, C;
        branch_operators(L_i[i], H_i[i], DTime, B_inv, C);
        const mat G = C*L_i[i];
        st.EV_i[i] = B_inv*EV_i_start[i] + G*eps_e;
        st.dEVdE_i[i] = G;
        st.dEVdT_i[i] = -G*alpha;
        st.sigma -= L_i[i]*st.EV_i[i];
        st.dSdE -= L_i[i]*G;
    }
    st.dSdE = 0.5*(st.dSdE + st.dSdE.t());   // L0 - sum L C L: symmetric up to round-off
    st.dSdT = -st.dSdE*alpha;
    return st;
}

} //namespace simcoon
