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

///@file eig_sym3.cpp
///@brief Fixed-size spectral decomposition of symmetric 3x3 matrices (cyclic Jacobi, Kopp 2008)
///@version 1.0

#include <cmath>
#include <utility>
#include <armadillo>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

bool eig_sym_3x3(vec::fixed<3> &lam, mat::fixed<3,3> &Q, const mat::fixed<3,3> &A_in) {

    // Working copy, symmetrised from the upper triangle (the lower one is never read)
    double A[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            A[i][j] = A_in(i, j);
            A[j][i] = A_in(i, j);
            if (!std::isfinite(A[i][j])) return false;
        }
    }
    double w[3] = {A[0][0], A[1][1], A[2][2]};
    double V[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    bool converged = false;
    for (int sweep = 0; sweep < 50; sweep++) {
        const double so = std::fabs(A[0][1]) + std::fabs(A[0][2]) + std::fabs(A[1][2]);
        if (so == 0.) {
            converged = true;
            break;
        }
        // threshold pivoting in the first sweeps, then every non-negligible pair
        const double thresh = (sweep < 4) ? 0.2*so/9. : 0.;
        for (int p = 0; p < 2; p++) {
            for (int q = p + 1; q < 3; q++) {
                const double apq = A[p][q];
                const double g = 100.*std::fabs(apq);
                if (sweep > 4 && std::fabs(w[p]) + g == std::fabs(w[p]) && std::fabs(w[q]) + g == std::fabs(w[q])) {
                    A[p][q] = A[q][p] = 0.;   // below the precision of both diagonal terms
                }
                else if (std::fabs(apq) > thresh) {
                    const double h = w[q] - w[p];
                    double t;
                    if (std::fabs(h) + g == std::fabs(h)) {
                        t = apq/h;
                    }
                    else {
                        const double theta = 0.5*h/apq;
                        t = 1./(std::fabs(theta) + std::sqrt(1. + theta*theta));
                        if (theta < 0.) t = -t;
                    }
                    const double c = 1./std::sqrt(1. + t*t);
                    const double s = t*c;
                    const double z = t*apq;
                    A[p][q] = A[q][p] = 0.;
                    w[p] -= z;
                    w[q] += z;
                    for (int r = 0; r < 3; r++) {
                        if (r == p || r == q) continue;
                        const double arp = A[r][p], arq = A[r][q];
                        A[r][p] = A[p][r] = c*arp - s*arq;
                        A[r][q] = A[q][r] = s*arp + c*arq;
                    }
                    for (int r = 0; r < 3; r++) {
                        const double vrp = V[r][p], vrq = V[r][q];
                        V[r][p] = c*vrp - s*vrq;
                        V[r][q] = s*vrp + c*vrq;
                    }
                }
            }
        }
    }
    if (!converged) return false;

    // ascending order, as arma::eig_sym
    int idx[3] = {0, 1, 2};
    if (w[idx[0]] > w[idx[1]]) std::swap(idx[0], idx[1]);
    if (w[idx[1]] > w[idx[2]]) std::swap(idx[1], idx[2]);
    if (w[idx[0]] > w[idx[1]]) std::swap(idx[0], idx[1]);
    for (int k = 0; k < 3; k++) {
        lam(k) = w[idx[k]];
        for (int r = 0; r < 3; r++) Q(r, k) = V[r][idx[k]];
    }
    return true;
}

mat::fixed<3,3> spectral_compose_3x3(const mat::fixed<3,3> &Q, const vec::fixed<3> &f) {
    mat::fixed<3,3> out;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
            double x = 0.;
            for (int k = 0; k < 3; k++) x += Q(i, k)*f(k)*Q(j, k);
            out(i, j) = x;
            out(j, i) = x;
        }
    }
    return out;
}

bool sqrtmat_sympd_3x3(mat::fixed<3,3> &out, const mat::fixed<3,3> &A) {
    vec::fixed<3> lam;
    mat::fixed<3,3> Q;
    if (!eig_sym_3x3(lam, Q, A) || lam(0) < 0.) return false;
    for (int k = 0; k < 3; k++) lam(k) = std::sqrt(lam(k));
    out = spectral_compose_3x3(Q, lam);
    return true;
}

bool logmat_sympd_3x3(mat::fixed<3,3> &out, const mat::fixed<3,3> &A) {
    vec::fixed<3> lam;
    mat::fixed<3,3> Q;
    if (!eig_sym_3x3(lam, Q, A) || !(lam(0) > 0.)) return false;
    for (int k = 0; k < 3; k++) lam(k) = std::log(lam(k));
    out = spectral_compose_3x3(Q, lam);
    return true;
}

bool expmat_sym_3x3(mat::fixed<3,3> &out, const mat::fixed<3,3> &A) {
    vec::fixed<3> lam;
    mat::fixed<3,3> Q;
    if (!eig_sym_3x3(lam, Q, A)) return false;
    for (int k = 0; k < 3; k++) lam(k) = std::exp(lam(k));
    out = spectral_compose_3x3(Q, lam);
    return true;
}

} //namespace simcoon
