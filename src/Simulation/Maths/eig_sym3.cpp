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
///@brief Fixed-size spectral decomposition of symmetric 3x3 matrices (Scherzinger-Dohrmann 2008, cyclic Jacobi fallback)
///@version 1.0

#include <algorithm>
#include <cmath>
#include <utility>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

namespace {

constexpr double eig_relative_gap = 1.e-2;

bool eig_jacobi_3x3(vec::fixed<3> &lam, mat::fixed<3,3> &Q, const mat::fixed<3,3> &A_in) {

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

// Scherzinger & Dohrmann (2008): the most distinct eigenvalue of the deviator, its eigenvector from
// two columns of S - eta1 I, then one Jacobi rotation of the 2x2 restriction on the complement
bool eig_deviatoric_3x3(vec::fixed<3> &lam, mat::fixed<3,3> &Q, const mat::fixed<3,3> &A_in) {
    const double m = (A_in(0,0) + A_in(1,1) + A_in(2,2))/3.;
    double S[3][3];
    double sc = 0.;
    for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++) {
            S[i][j] = S[j][i] = A_in(i, j) - (i == j ? m : 0.);
            sc = std::max(sc, std::fabs(S[i][j]));
        }
    const double isc = 1./sc;
    if (!std::isfinite(isc)) return false;   // subnormal deviator
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) S[i][j] *= isc;
    const double J2 = 0.5*(S[0][0]*S[0][0] + S[1][1]*S[1][1] + S[2][2]*S[2][2]) + S[0][1]*S[0][1] + S[0][2]*S[0][2] + S[1][2]*S[1][2];
    const double J3 = S[0][0]*(S[1][1]*S[2][2] - S[1][2]*S[1][2]) - S[0][1]*(S[0][1]*S[2][2] - S[1][2]*S[0][2])
                    + S[0][2]*(S[0][1]*S[1][2] - S[1][1]*S[0][2]);
    const double r = std::sqrt(J2/3.);
    const double c3a = std::max(-1., std::min(1., 0.5*J3/(r*r*r)));
    const double alpha = std::acos(c3a)/3.;
    const double eta1 = 2.*r*((alpha < simcoon::pi/6.) ? std::cos(alpha) : std::cos(alpha + 2.*simcoon::pi/3.));

    double R[3][3], nr[3];
    int imax = 0;
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) R[k][i] = S[i][k] - (i == k ? eta1 : 0.);
        nr[k] = R[k][0]*R[k][0] + R[k][1]*R[k][1] + R[k][2]*R[k][2];
        if (nr[k] > nr[imax]) imax = k;
    }
    double s1[3], s2[3], v1[3], t[2][3], nt[2];
    const double in1 = 1./std::sqrt(nr[imax]);
    for (int i = 0; i < 3; i++) s1[i] = R[imax][i]*in1;
    for (int a = 0, k = 0; k < 3; k++) {
        if (k == imax) continue;
        const double d = s1[0]*R[k][0] + s1[1]*R[k][1] + s1[2]*R[k][2];
        for (int i = 0; i < 3; i++) t[a][i] = R[k][i] - d*s1[i];
        nt[a] = t[a][0]*t[a][0] + t[a][1]*t[a][1] + t[a][2]*t[a][2];
        ++a;
    }
    const int b = (nt[0] >= nt[1]) ? 0 : 1;
    const double in2 = 1./std::sqrt(nt[b]);
    for (int i = 0; i < 3; i++) s2[i] = t[b][i]*in2;
    v1[0] = s1[1]*s2[2] - s1[2]*s2[1];
    v1[1] = s1[2]*s2[0] - s1[0]*s2[2];
    v1[2] = s1[0]*s2[1] - s1[1]*s2[0];

    double Ss1[3], Ss2[3];
    for (int i = 0; i < 3; i++) {
        Ss1[i] = S[i][0]*s1[0] + S[i][1]*s1[1] + S[i][2]*s1[2];
        Ss2[i] = S[i][0]*s2[0] + S[i][1]*s2[1] + S[i][2]*s2[2];
    }
    const double a11 = s1[0]*Ss1[0] + s1[1]*Ss1[1] + s1[2]*Ss1[2];
    const double a22 = s2[0]*Ss2[0] + s2[1]*Ss2[1] + s2[2]*Ss2[2];
    const double a12 = s1[0]*Ss2[0] + s1[1]*Ss2[1] + s1[2]*Ss2[2];
    double c = 1., sn = 0., tt = 0.;
    if (a12 != 0.) {
        const double theta = 0.5*(a22 - a11)/a12;
        tt = 1./(std::fabs(theta) + std::sqrt(1. + theta*theta));
        if (theta < 0.) tt = -tt;
        c = 1./std::sqrt(1. + tt*tt);
        sn = tt*c;
    }
    double v2[3], v3[3];
    for (int i = 0; i < 3; i++) {
        v2[i] = c*s1[i] - sn*s2[i];
        v3[i] = sn*s1[i] + c*s2[i];
    }
    const double w[3] = {eta1, a11 - tt*a12, a22 + tt*a12};
    const double *V[3] = {v1, v2, v3};
    int idx[3] = {0, 1, 2};
    if (w[idx[0]] > w[idx[1]]) std::swap(idx[0], idx[1]);
    if (w[idx[1]] > w[idx[2]]) std::swap(idx[1], idx[2]);
    if (w[idx[0]] > w[idx[1]]) std::swap(idx[0], idx[1]);
    for (int k = 0; k < 3; k++) {
        lam(k) = m + sc*w[idx[k]];
        for (int i = 0; i < 3; i++) Q(i, k) = V[idx[k]][i];
    }
    return true;
}

}  // namespace

bool eig_sym_3x3(vec::fixed<3> &lam, mat::fixed<3,3> &Q, const mat::fixed<3,3> &A) {
    for (int i = 0; i < 3; i++)
        for (int j = i; j < 3; j++)
            if (!std::isfinite(A(i, j))) return false;
    const int n_off = (A(0,1) != 0.) + (A(0,2) != 0.) + (A(1,2) != 0.);
    if (n_off <= 1) return eig_jacobi_3x3(lam, Q, A);   // one rotation at most, diagonal exact
    if (!eig_deviatoric_3x3(lam, Q, A)) return eig_jacobi_3x3(lam, Q, A);
    const double lmax = std::max(std::fabs(lam(0)), std::fabs(lam(2)));
    const double lmin = std::min(std::fabs(lam(0)), std::min(std::fabs(lam(1)), std::fabs(lam(2))));
    if (lmin < eig_relative_gap*lmax) return eig_jacobi_3x3(lam, Q, A);   // relative accuracy of the small eigenvalue
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
