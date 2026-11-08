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

///@file objective_rates.cpp
///@brief A set of function that help to define different quantities, depending on a selected objective rate
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <Fastor/Fastor.h>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Maths/eig_sym3.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

// Helper: copy arma 3x3 to Fastor tensor
// symmetric=true for stress/strain/C (memcpy sufficient), false for F/invF (needs transpose)
static inline Fastor::Tensor<double,3,3> to_fastor2(const mat &m, bool symmetric = true) {
    return arma_to_fastor2(mat::fixed<3,3>(m), symmetric);
}

// Helper: copy arma 6x6 Voigt to Fastor 3x3x3x3 (stiffness convention)
static inline Fastor::Tensor<double,3,3,3,3> to_fastor4(const mat &m) {
    return voigt_to_fastor4(mat::fixed<6,6>(m));
}

// Helper: symmetric 4th-order identity I_ijkl = 0.5*(δ_ik δ_jl + δ_il δ_jk)
static inline Fastor::Tensor<double,3,3,3,3> sym_identity_4() {
    Fastor::Tensor<double,3,3,3,3> I;
    I.zeros();
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
        I(i,j,i,j) += 0.5;
        I(i,j,j,i) += 0.5;
    }
    return I;
}

// Helper: 4th-order pull-back DSDE_LHMN = invF_Li invF_Hj invF_Mk invF_Nl C_ijkl
static inline Fastor::Tensor<double,3,3,3,3> pullback_4(
    const Fastor::Tensor<double,3,3,3,3> &C,
    const Fastor::Tensor<double,3,3> &invF)
{
    return push_forward_4(C, invF);
}

typedef mat::fixed<3,3> m33;

// Voigt index of the pair (i,j), and the pair (vi, vj) of a Voigt index
static const int vm[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
static const int vi[6] = {0, 1, 2, 0, 0, 1}, vj[6] = {0, 1, 2, 1, 2, 2};

// Fixed-size 3x3 products and Voigt-native contractions: the same products summed in the same order
// as the armadillo / Fastor expressions they replace (bit-identical under -ffp-contract=on; the
// tolerance-based tests remain the guarantee).
static inline m33 mm(const m33 &A, const m33 &B) {          // A*B
    m33 C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(0,j), x1 = B(1,j), x2 = B(2,j);
        for (int i = 0; i < 3; ++i) C(i,j) = A(i,0)*x0 + A(i,1)*x1 + A(i,2)*x2;
    }
    return C;
}
static inline m33 mtm(const m33 &A, const m33 &B) {         // A^T*B
    m33 C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(0,j), x1 = B(1,j), x2 = B(2,j);
        for (int i = 0; i < 3; ++i) C(i,j) = A(0,i)*x0 + A(1,i)*x1 + A(2,i)*x2;
    }
    return C;
}
static inline m33 mmt(const m33 &A, const m33 &B) {         // A*B^T
    m33 C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(j,0), x1 = B(j,1), x2 = B(j,2);
        for (int i = 0; i < 3; ++i) C(i,j) = A(i,0)*x0 + A(i,1)*x1 + A(i,2)*x2;
    }
    return C;
}

// alpha*X entrywise: the scaling armadillo applies after a product carrying a scalar factor
static inline m33 alpha_times(const double alpha, const m33 &X) {
    m33 Y;
    for (int e = 0; e < 9; ++e) Y(e) = alpha*X(e);
    return Y;
}

// 2nd-order centered velocity gradient L = (2/dt)(F1-F0)(F0+F1)^-1, as the armadillo expression
// (the inverse throws std::runtime_error on a singular F0+F1)
static inline m33 velocity_gradient(const mat &F0, const mat &F1, const double &DTime) {
    const m33 Fsum_inv = inv(F0+F1);
    return alpha_times(2./DTime, mm(m33(F1-F0), Fsum_inv));
}

// Helper: push-forward of a stiffness-convention Voigt 6x6, C'_isrp = F_iL F_sJ F_rM F_pN C_LJMN,
// as push_forward_4(voigt_to_fastor4(L), F): the four one-index contractions run in the same order
// over the same values, but only for the entries a Voigt result reads. C is symmetric in (M,N) exactly, and so are the first two partial
// contractions; only the pairs i <= s of the result are read. 198 contracted entries instead of
// 324, no 3x3x3x3 temporaries.
// R[is][r][p] = C'(i,s,r,p) for the Voigt pair is = (i <= s), unaveraged.
static void push_forward_voigt_core(double R[6][3][3], const mat::fixed<6,6> &L, const mat::fixed<3,3> &F) {
    // T1(i,J,MN) = sum_L F(i,L) C(L,J,M,N)
    double T1[3][3][6];
    for (int i = 0; i < 3; ++i)
    for (int J = 0; J < 3; ++J)
    for (int MN = 0; MN < 6; ++MN) {
        double sum = 0.0;
        for (int l = 0; l < 3; ++l)
            sum += F(i,l) * L(vm[l][J], MN);
        T1[i][J][MN] = sum;
    }
    // T2(is,MN) = sum_J F(s,J) T1(i,J,MN), i <= s
    double T2[6][6];
    for (int is = 0; is < 6; ++is)
    for (int MN = 0; MN < 6; ++MN) {
        const int i = vi[is], s = vj[is];
        double sum = 0.0;
        for (int J = 0; J < 3; ++J)
            sum += F(s,J) * T1[i][J][MN];
        T2[is][MN] = sum;
    }
    // T3(is,r,N) = sum_M F(r,M) T2(is,M,N)
    double T3[6][3][3];
    for (int is = 0; is < 6; ++is)
    for (int r = 0; r < 3; ++r)
    for (int N = 0; N < 3; ++N) {
        double sum = 0.0;
        for (int M = 0; M < 3; ++M)
            sum += F(r,M) * T2[is][vm[M][N]];
        T3[is][r][N] = sum;
    }
    // R(is,r,p) = sum_N F(p,N) T3(is,r,N)
    for (int is = 0; is < 6; ++is)
    for (int r = 0; r < 3; ++r)
    for (int p = 0; p < 3; ++p) {
        double sum = 0.0;
        for (int N = 0; N < 3; ++N)
            sum += F(p,N) * T3[is][r][N];
        R[is][r][p] = sum;
    }
}

// The Voigt average of the minor pair (r,p), as fastor4_to_voigt
static mat::fixed<6,6> voigt_average(const double R[6][3][3]) {
    mat::fixed<6,6> out;
    for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J)
            out(I, J) = 0.5 * (R[I][vi[J]][vj[J]] + R[I][vj[J]][vi[J]]);
    return out;
}

// fastor4_to_voigt(push_forward_4(voigt_to_fastor4(L), F))
static mat::fixed<6,6> push_forward_voigt(const mat::fixed<6,6> &L, const mat::fixed<3,3> &F) {
    double R[6][3][3];
    push_forward_voigt_core(R, L, F);
    return voigt_average(R);
}

// The Jaumann correction of apply_jaumann_correction on the entries a Voigt result reads:
// R[is][p][r] += sign*corr(i,s,p,r), the same terms added in the same order. tau_f is indexed as
// the Fastor copy of tau (to_fastor2 with the symmetric memcpy: tau_f(a,b) = tau(b,a)).
static void add_jaumann_correction(double R[6][3][3], const mat::fixed<3,3> &tau_f, const bool add) {
    const double sign = add ? 0.5 : -0.5;
    for (int is = 0; is < 6; ++is) {
        const int i = vi[is], s = vj[is];
        for (int p = 0; p < 3; ++p)
        for (int r = 0; r < 3; ++r) {
            double corr = 0.0;
            if (i == r) corr += tau_f(p,s);
            if (i == p) corr += tau_f(r,s);
            if (s == p) corr += tau_f(i,r);
            if (s == r) corr += tau_f(i,p);
            R[is][p][r] += sign * corr;
        }
    }
}

// R[is][p][r] = A(i,s,p,r) of voigt_to_fastor4(L), for the Voigt pairs is = (i <= s)
static void voigt_expand(double R[6][3][3], const mat::fixed<6,6> &L) {
    for (int is = 0; is < 6; ++is)
        for (int p = 0; p < 3; ++p)
            for (int r = 0; r < 3; ++r)
                R[is][p][r] = L(is, vm[p][r]);
}

// The Voigt 6x6 of R[is][p][r] for the pairs p <= r (a result that keeps both minor symmetries exactly)
static mat::fixed<6,6> voigt_collapse(const double R[6][3][3]) {
    mat::fixed<6,6> out;
    for (int is = 0; is < 6; ++is)
        for (int p = 0; p < 3; ++p)
            for (int r = p; r < 3; ++r)
                out(is, vm[p][r]) = R[is][p][r];
    return out;
}

// The B-correction of apply_B_correction on the entries a Voigt result reads (pairs i <= j):
// R[ij][k][l] = A(i,j,k,l) + sign*(sum1 + sum2), with B - I and the sums as there; A and B are
// stiffness-convention Voigt 6x6, tau_f indexed as the Fastor copy of tau.
static void add_B_correction(double R[6][3][3], const mat::fixed<6,6> &A, const mat::fixed<6,6> &B, const mat::fixed<3,3> &tau_f, const bool add) {
    // B - I, I the symmetric 4th-order identity of sym_identity_4 (entries 0, 0.5 or 1)
    auto BmI = [&](const int a, const int b, const int c, const int d) {
        double I = 0.;
        if (c == a && d == b) I += 0.5;
        if (c == b && d == a) I += 0.5;
        return B(vm[a][b], vm[c][d]) - I;
    };
    const double sign = add ? 1.0 : -1.0;
    for (int ij = 0; ij < 6; ++ij) {
        const int i = vi[ij], j = vj[ij];
        for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int p = 0; p < 3; ++p) {
                sum1 += BmI(i,p,k,l) * tau_f(p,j);
                sum2 += tau_f(i,p) * BmI(j,p,k,l);
            }
            double r = A(ij, vm[k][l]);
            r += sign * (sum1 + sum2);
            R[ij][k][l] = r;
        }
    }
}

// The Fastor copy of a symmetric 3x3 (to_fastor2 with the symmetric memcpy), as an arma matrix
static inline mat::fixed<3,3> fastor_sym_copy(const mat &tau) {
    return mat::fixed<3,3>(tau).t();
}

static inline bool is_exactly_symmetric(const mat &X) {
    return X(0,1) == X(1,0) && X(0,2) == X(2,0) && X(1,2) == X(2,1);
}

// Helper: B-correction for objective rate conversions
// result_ijkl = A_ijkl + (B_ipkl - I_ipkl)*τ_pj + τ_ip*(B_jpkl - I_jpkl)
static inline Fastor::Tensor<double,3,3,3,3> apply_B_correction(
    const Fastor::Tensor<double,3,3,3,3> &A,
    const Fastor::Tensor<double,3,3,3,3> &B,
    const Fastor::Tensor<double,3,3> &tau,
    bool add)
{
    auto I = sym_identity_4();
    Fastor::Tensor<double,3,3,3,3> BmI = B - I;
    Fastor::Tensor<double,3,3,3,3> result = A;
    double sign = add ? 1.0 : -1.0;

    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l) {
        double sum1 = 0.0, sum2 = 0.0;
        for (int p = 0; p < 3; ++p) {
            sum1 += BmI(i,p,k,l) * tau(p,j);
            sum2 += tau(i,p) * BmI(j,p,k,l);
        }
        result(i,j,k,l) += sign * (sum1 + sum2);
    }
    return result;
}

// Helper: Jaumann correction
// result_ispr = A_ispr ± 0.5*(τ_ps δ_ir + τ_rs δ_ip + τ_ir δ_sp + τ_ip δ_sr)
static inline Fastor::Tensor<double,3,3,3,3> apply_jaumann_correction(
    const Fastor::Tensor<double,3,3,3,3> &A,
    const Fastor::Tensor<double,3,3> &tau,
    bool add)
{
    Fastor::Tensor<double,3,3,3,3> result = A;
    double sign = add ? 0.5 : -0.5;

    for (int i = 0; i < 3; ++i)
    for (int s = 0; s < 3; ++s)
    for (int p = 0; p < 3; ++p)
    for (int r = 0; r < 3; ++r) {
        double corr = 0.0;
        if (i == r) corr += tau(p,s);
        if (i == p) corr += tau(r,s);
        if (s == p) corr += tau(i,r);
        if (s == r) corr += tau(i,p);
        result(i,s,p,r) += sign * corr;
    }
    return result;
}

// Spectral decomposition of B = F F^T (eigenvalues bi, eigenvectors the columns of Bi) with the
// fixed-size 3x3 solver; exception_eig_sym on failure, as arma::eig_sym did.
static inline void eig_B(vec::fixed<3> &bi, mat::fixed<3,3> &Bi, const mat &F, const char *where) {
    const mat::fixed<3,3> F_(F);
    if (!eig_sym_3x3(bi, Bi, F_*F_.t())) {
        throw simcoon::exception_eig_sym(std::string("Error in eig_sym function inside ") + where + ".");
    }
}

// Spin coefficients of the logarithmic rate (Xiao-Bruhns-Meyers), f(r) of r = b_i/b_j
static inline double log_spin_f(const double &bi, const double &bj) {
    const double s = bi/bj - 1.;
    if (fabs(s) > 1.e-4) {
        return (1.+(bi/bj))/(1.-(bi/bj))+2./log(bi/bj);
    }
    // Taylor expansion: f(r) = -s/6 + s^2/12 + O(s^3) where s = r - 1
    return s*(-1./6. + s/12.);
}

// g(r) = (1 - sqrt r)/(1 + sqrt r), r = b_i/b_j
static inline double log_spin_g(const double &bi, const double &bj) {
    const double s = bi/bj - 1.;
    if (fabs(s) > 1.e-4) {
        return (1.-(pow(bi/bj,0.5)))/(1.+(pow(bi/bj,0.5)));
    }
    // Taylor expansion: g(r) = -s/4 + 3s^2/16 + O(s^3) where s = r - 1
    return s*(-1./4. + 3.*s/16.);
}

// t/sinh(t), t = 1/2 ln(b_i/b_j): the A^R (log_R) Daleckii-Krein kernel
static inline double kappa_R(const double &bi, const double &bj) {
    const double t = 0.5*log(bi/bj);            // ln(lambda_i/lambda_j); 0 on the diagonal
    if (fabs(t) > 1.e-4) {
        return t/sinh(t);
    }
    return 1. - t*t/6. + 7.*t*t*t*t/360.;       // Taylor of t/sinh(t)
}

// N_1 and N_2 of the logarithmic rates on the eigenbasis (bi, Bi) of B
static inline void log_spin_N1_N2(mat &N_1, mat &N_2, const vec::fixed<3> &bi, const mat::fixed<3,3> &Bi, const mat &D) {
    mat::fixed<3,3> f, g;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            f(i, j) = (i != j) ? log_spin_f(bi(i), bi(j)) : 0.;
            g(i, j) = (i != j) ? log_spin_g(bi(i), bi(j)) : 0.;
        }
    }
    N_1 = spectral_congruence(Bi, f, D);
    N_2 = spectral_congruence(Bi, g, D);
}

void Jaumann(mat &DR, mat &D, mat &W, const double &DTime, const mat &F0, const mat &F1) {
    mat L = zeros(3,3);   // DTime <= iota: D = W = 0, DR = I (rate undefined, see finite_rotation)
    if(DTime > simcoon::iota) {
        try {
            // 2nd-order midpoint L = (2/dt)(F1-F0)(F0+F1)^-1, shared by every rate so D = sym(L) agrees
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Jaumann (L).");
        }          
    }   

    //decomposition of L
    D = 0.5*(L+L.t());
    W = 0.5*(L-L.t());
    
    DR = Hughes_Winget(W, DTime);
}
    
void Green_Naghdi(mat &DR, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    //Green-Naghdi
    mat L = zeros(3,3);
    if(DTime > simcoon::iota) {    
        try {
            // Same 2nd-order centered velocity gradient as the other rate functions (see Jaumann),
            // so that D=sym(L) matches across all rates.
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Green_Naghdi (L).");
        }          
    }   
    
    //decomposition of L
    D = 0.5*(L+L.t());
    // DR = R1 R0^T exactly; Hughes_Winget(Omega) == DR, no reconstruction needed
    finite_rotation(F0, F1, DTime, DR, Omega);
}

namespace {

// one eigen-solve of C0 and one of C1
void logarithmic_R_kinematics(mat &DR, mat &D, mat &Omega, vec::fixed<3> &b1, mat::fixed<3,3> &Nb, const double &DTime, const mat &F0, const mat &F1) {
    mat L = zeros(3,3);
    if(DTime > simcoon::iota) {    
        try {
            // 2nd-order centered velocity gradient (see Jaumann); D=sym(L) is shared by all rates.
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic_R (L).");
        }          
    }   
    
    //decomposition of L
    D = 0.5*(L+L.t());
    // Exact relative polar rotation + inverse-Cayley spin (see Green_Naghdi):
    // the exact DR is what makes the log_R tangent transport equivariant.
    mat::fixed<3,3> R0, R1, N0, N1;
    vec::fixed<3> c0;
    RU_spectral(R0, c0, N0, F0);
    RU_spectral(R1, b1, N1, F1);
    DR = mmt(R1, R0);
    Nb = mm(R1, N1);   // B1 = R1 C1 R1^T

    if (DTime <= simcoon::iota) {
        Omega = zeros(3,3);
        return;
    }
    // Inverse Cayley, as finite_rotation
    mat::fixed<3,3> DRpI = DR + eye(3,3), DRpI_inv;
    if (!inv(DRpI_inv, DRpI)) {
        throw simcoon::exception_inv("Error in inv function inside finite_rotation.");
    }
    Omega = alpha_times(2./DTime, mm(m33(DR-eye(3,3)), DRpI_inv));
}

// De = A^R : D on the eigenbasis (b, Nb) of B, by congruence (t/sinh t kernel)
mat::fixed<3,3> A_R_D(const vec::fixed<3> &b, const mat::fixed<3,3> &Nb, const mat &D) {
    mat::fixed<3,3> K;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            K(i, j) = (i != j) ? kappa_R(b(i), b(j)) : 1.;   // kappa_R(b, b) = 1 exactly
    return spectral_congruence(Nb, K, D);
}

}  // namespace

void logarithmic_R(mat &DR, mat &N_1, mat &N_2, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    logarithmic_R_kinematics(DR, D, Omega, bi, Bi, DTime, F0, F1);
    log_spin_N1_N2(N_1, N_2, bi, Bi, D);
}

void logarithmic_F(mat &DF, mat &N_1, mat &N_2, mat &D, mat &L, const double &DTime, const mat &F0, const mat &F1) {
    L = zeros(3,3);   // DTime <= iota: D = 0, DF = I (L is an output; never leave it stale)
    if(DTime > simcoon::iota) {
        try {
            // 2nd-order centered velocity gradient (see Jaumann): D=sym(L) is shared by all rates,
            // and exp(L*dt) approximates F1 F0^-1 used below for the transport DF.
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic_F (L).");
        }
    }

    //decomposition of L
    D = 0.5*(L+L.t());

    //Logarithmic
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F1, "logarithmic_F");
    log_spin_N1_N2(N_1, N_2, bi, Bi, D);
    
    DF = Hughes_Winget(L, DTime);
}

void Truesdell(mat &DF, mat &D, mat &L, const double &DTime, const mat &F0, const mat &F1) {
    L = zeros(3,3);   // DTime <= iota: D = 0, DF = I (L is an output; never leave it stale)
    if(DTime > simcoon::iota) {
        try {
            // 2nd-order centered velocity gradient (see Jaumann): D=sym(L) is shared by all rates,
            // and exp(L*dt) approximates F1 F0^-1 used below for the transport DF.
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Truesdell (L).");
        }          
    }      

    //Note that The "spin" is actually L (spin for rigid frames of reference, "flot" for Truesdell)    
    D = 0.5*(L+L.t());

    // The convected transport is exact: DF = F1 F0^-1, so the closed-form Almansi increment of
    // Delta_log_strain_corate recovers e_A(F1) exactly.
    try {
        DF = F1*inv(F0);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Truesdell (DF).");
    }
}

mat Hughes_Winget(const mat &Omega, const double &DTime) {
    mat I = eye(3,3);
    mat DR;
    try {
        DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Hughes_Winget.");
    }
    return DR;
}

mat get_BBBB(const mat &F1) {
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F1, "get_BBBB");
    mat::fixed<6,6> BBBB(fill::zeros);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if (i!=j) {
                add_linearop_eigsym(BBBB, log_spin_f(bi(i), bi(j)), Bi.col(i), Bi.col(j));
            }
        }
    }
    return BBBB;
}

mat get_BBBB_GN(const mat &F1) {
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F1, "get_BBBB_GN");
    mat::fixed<6,6> BBBB(fill::zeros);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if (i!=j) {
                const double f_z = (sqrt(bi(j)) - sqrt(bi(i)))/(sqrt(bi(j)) + sqrt(bi(i)));
                add_linearop_eigsym(BBBB, f_z, Bi.col(i), Bi.col(j));
            }
        }
    }
    return BBBB;
}

// A^R: log_R-frame strain-concentration tensor, De = A^R:D (full doc in objective_rates.hpp).
// t/sinh(t) geometric-mean Daleckii-Krein kernel on the eigenbasis of B = F F^T.
mat A_R(const mat &F) {
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F, "A_R");
    mat::fixed<6,6> AR(fill::zeros);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            add_linearop_eigsym(AR, kappa_R(bi(i), bi(j)), Bi.col(i), Bi.col(j));
        }
    }
    AR.rows(3,5) *= 2.0;    // tensor -> engineering strain-concentration convention: A^R(I)=I, read De with v2t_strain
    return AR;
}

// kernel t*coth(t): unit diagonal, so A^F:D integrates to ln V like A^R:D
mat A_F(const mat &F) {
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F, "A_F");
    mat::fixed<6,6> AF(fill::zeros);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            double t = (i != j) ? 0.5*log(bi(i)/bi(j)) : 0.;   // log(b/b) = 0 exactly
            double tcoth;
            if (fabs(t) > 1.e-4) {
                tcoth = t/tanh(t);
            } else {
                tcoth = 1. + t*t/3. - t*t*t*t/45.;       // Taylor of t*coth(t)
            }
            double c = tcoth;     // t*coth(t): F-frame strain-concentration kernel (diagonal -> 1, recovers ln V)
            add_linearop_eigsym(AF, c, Bi.col(i), Bi.col(j));
        }
    }
    AF.rows(3,5) *= 2.0;    // tensor -> engineering strain-concentration convention: A^F(I)=I, read De with v2t_strain
    return AF;
}

void logarithmic(mat &DR, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    mat L = zeros(3,3);

    if(DTime > simcoon::iota) {    
        try {
            // 2nd-order centered velocity gradient (see Jaumann); D=sym(L) is shared by all rates.
            L = velocity_gradient(F0, F1, DTime);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic (L).");
        }          
    }    
        
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    
    //Logarithmic
    vec::fixed<3> bi;
    mat::fixed<3,3> Bi;
    eig_B(bi, Bi, F1, "logarithmic");
    mat::fixed<3,3> f;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            f(i, j) = (i != j) ? log_spin_f(bi(i), bi(j)) : 0.;
    const mat::fixed<3,3> N = spectral_congruence(Bi, f, D);
    Omega = W + N;
    DR = Hughes_Winget(Omega, DTime);
}

mat Delta_log_strain(const mat &D, const mat &Omega, const double &DTime) {
    const m33 DR = Hughes_Winget(Omega, DTime), D_(D);
    return 0.5*(D_+mmt(mm(DR, D_), DR))*DTime;
}

mat Delta_log_strain_F(const mat &D, const mat &L, const double &DTime) {
    // DF is non-orthogonal: the rotated term is DF*D*inv(DF), not DF*D*DF^T
    const mat DF = Hughes_Winget(L, DTime);
    return 0.5*(D+(DF*D*inv(DF)))*DTime;
}

// Corate spin dispatch (full doc in objective_rates.hpp): set DR + rate D + spin/L Omega for the
// chosen objective rate. Single source of truth for the solver's control_type ladders.
void corate_kinematics(const int &corate_type, mat &DR, mat &D, mat &Omega, const mat &F0, const mat &F1, const double &DTime) {
    mat N_1 = zeros(3,3), N_2 = zeros(3,3);
    vec::fixed<3> b1;
    mat::fixed<3,3> Nb;
    switch (corate_type) {
        case 0: Jaumann(DR, D, Omega, DTime, F0, F1); break;
        case 1: Green_Naghdi(DR, D, Omega, DTime, F0, F1); break;
        case 2: logarithmic(DR, D, Omega, DTime, F0, F1); break;
        case 3: logarithmic_R_kinematics(DR, D, Omega, b1, Nb, DTime, F0, F1); break;   // DR = R-rotation (N_1, N_2 unused)
        case 4: Truesdell(DR, D, Omega, DTime, F0, F1); break;                 // DR = DF, Omega = L
        case 5: logarithmic_F(DR, N_1, N_2, D, Omega, DTime, F0, F1); break;   // DR = DF, Omega = L
        case corate_type_total: Jaumann(DR, D, Omega, DTime, F0, F1); DR.eye(3,3); Omega.zeros(3,3); break;   // strain space: no frame increment
        default: break;
    }
}

// One-step kinematics plus the strain increment, sharing the spectral work (doc in objective_rates.hpp)
void corate_kinematics(const int &corate_type, mat &DR, mat &D, mat &Omega, mat &De, const mat &F0, const mat &F1, const double &DTime) {
    if (corate_type == 3) {   // log_R: the eigen-solve of C1 gives R1, the eigenbasis of B1 and A^R
        vec::fixed<3> b1;
        mat::fixed<3,3> Nb;
        logarithmic_R_kinematics(DR, D, Omega, b1, Nb, DTime, F0, F1);
        const mat::fixed<3,3> Dh = A_R_D(b1, Nb, D);
        const m33 DR_(DR);
        De = 0.5*(Dh + mmt(mm(DR_, Dh), DR_))*DTime;
        return;
    }
    corate_kinematics(corate_type, DR, D, Omega, F0, F1, DTime);
    De = Delta_log_strain_corate(F0, F1, DR, D, Omega, DTime, corate_type);
}

namespace {

// 1/2 ln of the SPD matrix A (a Cauchy-Green tensor)
mat::fixed<3,3> half_logmat(const mat::fixed<3,3> &A, const char *where) {
    mat::fixed<3,3> H;
    if (!logmat_sympd_3x3(H, A)) {
        throw simcoon::exception_logmat_sympd(std::string("Error in logmat_sympd function inside ") + where + ".");
    }
    return 0.5*H;
}

}  // namespace

// Dispatch documented in objective_rates.hpp.
mat Delta_log_strain_corate(const mat &F0, const mat &F1, const mat &DR, const mat &D, const mat &Omega, const double &DTime, const int &corate_type) {
    if (corate_type == corate_type_total) {   // total formulation: H = 1/2 ln C is a state function of F, reference components
        const mat::fixed<3,3> F0_(F0), F1_(F1);
        return half_logmat(F1_.t()*F1_, "Delta_log_strain_corate") - half_logmat(F0_.t()*F0_, "Delta_log_strain_corate");
    }
    if (corate_type == 5) {   // log_F: convected A^F:D rate (Omega carries the velocity gradient L)
        return Delta_log_strain_F(v2t_strain(A_F(F1)*t2v_strain(D)), Omega, DTime);
    }
    if (corate_type == 2) {   // XBM: exact closed-form spatial log-strain difference -> etot = ln V1
        const mat::fixed<3,3> F0_(F0), F1_(F1);
        const mat::fixed<3,3> lnV0 = half_logmat(F0_*F0_.t(), "Delta_log_strain_corate");
        const mat::fixed<3,3> lnV1 = half_logmat(F1_*F1_.t(), "Delta_log_strain_corate");
        return lnV1 - DR*lnV0*DR.t();                          // orthogonal (DR^T) transport
    }
    // log_R spin is carried by sv_M->DR: the solver must not also apply DR_N (double count)
    // frame increment is DR itself (inverse Cayley of Omega)
    if (corate_type == 3) {   // log_R: A^R:D
        vec::fixed<3> b1;
        mat::fixed<3,3> Nb;
        eig_B(b1, Nb, F1, "Delta_log_strain_corate");
        const mat::fixed<3,3> Dh = A_R_D(b1, Nb, D);
        return 0.5*(Dh + DR*Dh*DR.t())*DTime;
    }
    if (corate_type == 4) {   // Truesdell: e_A(F1) = DF^-T e_A(F0) DF^-1 + De, exactly, with DR = DF
        return 0.5*(eye(3,3) - inv_sympd(DR*DR.t()));
    }
    return Delta_log_strain(D, Omega, DTime);   // Jaumann / GN
}

double Delta_work_conjugacy(const vec &tau_start, const vec &tau_start_tr, const vec &tau, const vec &Detot, const mat &F0, const mat &F1, const int &corate_type) {
    if (!work_correction_applies(corate_type))
        return 0.;
    // fixed size and the bool inv(): no heap, no throw (sim.umat calls it in a parallel region)
    mat::fixed<3,3> Fsum_inv;
    if (!inv(Fsum_inv, mat::fixed<3,3>(F1 + F0)))
        return 0.;
    const mat::fixed<3,3> LDt = alpha_times(2., mm(m33(F1 - F0), Fsum_inv));   // midpoint velocity gradient times DTime
    // both stresses in the lab frame: a start stress carried to the end frame against a midpoint
    // D would leave a first-order [W, tau] : D error, even for isotropic elasticity
    return 0.5*dot(tau_start + tau, t2v_strain(0.5*(LDt + LDt.t())))
         - 0.5*dot(tau_start_tr + tau, Detot);
}

// Exact log-box <-> material tangent maps (contract: DtauDe_2_DSDE / DSDE_2_DtauDe)

namespace {

// Daleckii-Krein in congruence form: dY = N (G % N^T dX N) N^T
struct dk_spectral {
    arma::vec::fixed<3> lam;      // eigenvalues of C (SPD)
    arma::mat::fixed<3,3> N;      // eigenvectors of C (columns)
    arma::mat::fixed<3,3> G_h;    // DK coefficients of h(C) = 1/2 ln C
    arma::mat::fixed<3,3> G_iU;   // DK coefficients of C^{-1/2}
    arma::mat::fixed<3,3> U;      // C^{1/2} (built when asked: only the inverse map reads it)
    arma::mat::fixed<3,3> invU;   // C^{-1/2}
    arma::mat::fixed<3,3> R;      // polar rotation F U^{-1}
    arma::mat::fixed<3,3> tau_hat;  // corotational Kirchhoff stress R^T tau R

    void build(const arma::mat &F, const arma::mat &tau, const bool &with_U = true) {
        const arma::mat::fixed<3,3> F_(F);
        const bool success_eig_sym = eig_sym_3x3(lam, N, F_.t()*F_);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym(
                "Error in eig_sym function inside the exact DtauDe<->DSDE map.");
        }
        // near-singular F: throw exception_inv, the solver step-cut signal, rather than NaN tangents
        if (!lam.is_finite() || lam.min() <= simcoon::iota) {
            throw simcoon::exception_inv(
                "Near-singular F inside the exact DtauDe<->DSDE map (eigenvalue of C <= iota).");
        }
        G_h  = gmat([](double l) { return 0.5 * std::log(l); },
                    [](double l) { return 0.5 / l; });
        G_iU = gmat([](double l) { return 1.0 / std::sqrt(l); },
                    [](double l) { return -0.5 * std::pow(l, -1.5); });
        if (with_U) U = spectral_compose_3x3(N, sqrt(lam));
        invU = spectral_compose_3x3(N, 1.0 / sqrt(lam));
        R = mm(F_, invU);
        tau_hat = mm(mtm(R, arma::mat::fixed<3,3>(tau)), R);
    }

    // f' at the midpoint under the 1e-4 relative coincidence guard (as get_BBBB / A_R)
    template <typename Ff, typename Fdf>
    arma::mat gmat(Ff f, Fdf df) const {
        const arma::vec fl = {f(lam(0)), f(lam(1)), f(lam(2))};
        arma::mat G(3, 3);
        for (int a = 0; a < 3; ++a) {
            G(a, a) = df(lam(a));
            for (int b = a + 1; b < 3; ++b) {
                G(a, b) = (std::abs(lam(a) / lam(b) - 1.0) > 1.0e-4)
                    ? (fl(a) - fl(b)) / (lam(a) - lam(b))
                    : df(0.5 * (lam(a) + lam(b)));
                G(b, a) = G(a, b);
            }
        }
        return G;
    }

};

// Fixed-size Voigt helpers of the 6-column tangent maps below: per column, 3x3 products on the
// stack and the 6x6 box applied by hand (a BLAS gemv call per column otherwise).

// engineering-Voigt unit strain e_j as a tensor (v2t_strain of the j-th basis vector)
inline m33 unit_strain_tensor(const int j) {
    m33 e(arma::fill::zeros);
    const double v = (j < 3) ? 1. : 0.5;
    e(vi[j], vj[j]) = v;
    e(vj[j], vi[j]) = v;
    return e;
}

// y = Lt t2v_strain(e)
inline void box_apply(const arma::mat::fixed<6,6> &Lt, const m33 &e, double y[6]) {
    const double x[6] = {e(0,0), e(1,1), e(2,2), e(0,1) + e(1,0), e(0,2) + e(2,0), e(1,2) + e(2,1)};
    for (int i = 0; i < 6; ++i) {
        double yi = 0.;
        for (int k = 0; k < 6; ++k) yi += Lt(i, k)*x[k];
        y[i] = yi;
    }
}

inline m33 v2t_stress_arr(const double v[6]) {
    m33 S;
    S(0,0) = v[0]; S(1,1) = v[1]; S(2,2) = v[2];
    S(0,1) = S(1,0) = v[3];
    S(0,2) = S(2,0) = v[4];
    S(1,2) = S(2,1) = v[5];
    return S;
}

// column j of a 6x6 Voigt map, as t2v_stress(S)
inline void set_stress_col(arma::mat &M, const int j, const m33 &S) {
    M(0, j) = S(0,0); M(1, j) = S(1,1); M(2, j) = S(2,2);
    M(3, j) = 0.5*(S(0,1) + S(1,0));
    M(4, j) = 0.5*(S(0,2) + S(2,0));
    M(5, j) = 0.5*(S(1,2) + S(2,1));
}

}  // namespace

namespace {

// S = 2 N (G % That) N^T, tau = R N (kappa % That) N^T R^T, That = N^T T N, G_ab = f[c_a, c_b],
// f = 1/2 ln, kappa = 2 lam_a lam_b G_ab = t/sinh(t) (A_R kernel)
struct logU_spectral {
    arma::vec::fixed<3> c;              // eigenvalues of C (SPD)
    arma::mat::fixed<3,3> N;            // eigenvectors of C (columns)
    arma::mat::fixed<3,3> R;            // polar rotation F U^{-1}
    arma::mat::fixed<3,3> kappa;        // t/sinh(t)
    arma::mat::fixed<3,3> G;            // first divided differences of 1/2 ln
    double Gam[3][3][3];                // second divided differences of 1/2 ln

    // f[x, y]: no difference of logs, hence no cancellation at any gap (Taylor below 1e-4 as A_R)
    static double kappa_of(const double &x, const double &y) {
        const double t = 0.5*std::log(x/y);
        return (std::abs(t) > 1.e-4) ? t/std::sinh(t) : 1. - t*t/6. + 7.*t*t*t*t/360.;
    }
    static double dd1(const double &x, const double &y) {
        return kappa_of(x, y)/(2.*std::sqrt(x*y));
    }
    // f[x,y,z]: separated -> difference over the widest gap; clustered (< 1e-4) -> Taylor about the mean to h_4
    static double dd2(double x, double y, double z) {
        if (x > y) std::swap(x, y);
        if (y > z) std::swap(y, z);
        if (x > y) std::swap(x, y);
        const double m = (x + y + z)/3.;
        if (z - x > 1.e-4*m) {
            return (dd1(x, y) - dd1(y, z))/(x - z);
        }
        const double dx = (x - m)/m, dy = (y - m)/m, dz = (z - m)/m;
        const double s2 = dx*dx + dy*dy + dz*dz;
        const double s3 = dx*dx*dx + dy*dy*dy + dz*dz*dz;
        return -(1. + s2/4. - 2.*s3/15. + s2*s2/12.)/(4.*m*m);
    }

    // tangent = false: the stress maps need kappa only
    void build(const arma::mat &F, const bool &tangent) {
        const arma::mat::fixed<3,3> F_(F);
        const bool success_eig_sym = eig_sym_3x3(c, N, F_.t()*F_);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside the total-formulation (strain-space) maps.");
        }
        // same recoverable failure as dk_spectral: a step cut, not NaN stresses
        if (!c.is_finite() || c.min() <= simcoon::iota) {
            throw simcoon::exception_inv("Near-singular F inside the total-formulation (strain-space) maps (eigenvalue of C <= iota).");
        }
        R = F_*spectral_compose_3x3(N, 1.0/sqrt(c));
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                kappa(a, b) = (a != b) ? kappa_of(c(a), c(b)) : 1.;   // kappa_of(c, c) = 1 exactly
                if (!tangent) continue;
                G(a, b) = kappa(a, b)/(2.*std::sqrt(c(a)*c(b)));
            }
        }
        if (!tangent) return;
        // dd2 sorts its arguments, so Gam is the same, to the bit, for every permutation of
        // (a, b, k): the 10 index multisets are evaluated once each, instead of 27 calls
        for (int a = 0; a < 3; ++a) {
            for (int b = a; b < 3; ++b) {
                for (int k = b; k < 3; ++k) {
                    const double g = dd2(c(a), c(b), c(k));
                    Gam[a][b][k] = Gam[a][k][b] = Gam[b][a][k] = Gam[b][k][a] = Gam[k][a][b] = Gam[k][b][a] = g;
                }
            }
        }
    }

    // T in the eigenbasis of C, from the Kirchhoff stress: entrywise, amplification sinh(t)/t only
    arma::mat That_of_tau(const arma::mat &tau) const {
        return (N.t()*R.t()*tau*R*N)/kappa;
    }

    // D2f(C)[T, dC] in the eigenbasis
    arma::mat::fixed<3,3> geometric(const arma::mat &That, const arma::mat &dCp) const {
        arma::mat::fixed<3,3> X;
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double x = 0.;
                for (int k = 0; k < 3; ++k) {
                    x += Gam[a][k][b]*(That(a, k)*dCp(k, b) + dCp(a, k)*That(k, b));
                }
                X(a, b) = x;
            }
        }
        return X;
    }
};

// dS/dE = P^T Lt P + T:L of the total formulation. The box Lt = dT/dH acts in reference
// components, so unlike DtauDe_2_DSDE there is no R conjugation.
mat logU_box_2_DSDE(const mat &Lt, const mat &F, const mat &tau) {
    logU_spectral sp;
    sp.build(F, true);
    const mat That = sp.That_of_tau(tau);
    const arma::mat::fixed<6,6> Lt_(Lt);
    mat DSDE(6, 6);
    vec ej(6);
    for (int j = 0; j < 6; ++j) {
        ej.zeros(); ej(j) = 1.;
        const mat dCp = sp.N.t()*(2.*v2t_strain(ej))*sp.N;
        const mat dH = sp.N*(sp.G % dCp)*sp.N.t();
        double y[6];
        box_apply(Lt_, dH, y);   // Lt t2v_strain(dH) by hand (a BLAS gemv call per column otherwise)
        const mat dTp = sp.N.t()*v2t_stress_arr(y)*sp.N;
        DSDE.col(j) = t2v_stress(2.*sp.N*(sp.G % dTp + sp.geometric(That, dCp))*sp.N.t());
    }
    return DSDE;
}

// Exact inverse: G > 0 entrywise, so every step of the forward map is undone in turn.
mat DSDE_2_logU_box(const mat &DSDE, const mat &F, const mat &tau) {
    logU_spectral sp;
    sp.build(F, true);
    const mat That = sp.That_of_tau(tau);
    const arma::mat::fixed<6,6> DSDE_(DSDE);
    mat Lt(6, 6);
    vec ej(6);
    for (int j = 0; j < 6; ++j) {
        ej.zeros(); ej(j) = 1.;
        const mat dCp = (sp.N.t()*v2t_strain(ej)*sp.N)/sp.G;
        double y[6];
        box_apply(DSDE_, m33(0.5*sp.N*dCp*sp.N.t()), y);   // DSDE t2v_strain(.) by hand, as above
        const mat dS = v2t_stress_arr(y);
        const mat dTp = (0.5*sp.N.t()*dS*sp.N - sp.geometric(That, dCp))/sp.G;
        Lt.col(j) = t2v_stress(sp.N*dTp*sp.N.t());
    }
    return Lt;
}

}  // namespace

mat LogStressU2Kirchoff(const mat &T, const mat &F) {
    logU_spectral sp;
    sp.build(F, false);
    const mat RN = sp.R*sp.N;
    return RN*(sp.kappa % (sp.N.t()*T*sp.N))*RN.t();
}

mat Kirchoff2LogStressU(const mat &tau, const mat &F) {
    logU_spectral sp;
    sp.build(F, false);
    return sp.N*sp.That_of_tau(tau)*sp.N.t();
}

mat DtauDe_2_DSDE(const mat &Lt, const mat &F, const mat &tau){

    dk_spectral sp;
    sp.build(F, tau, false);
    const arma::mat::fixed<6,6> Lt_(Lt);
    const m33 invU_tau_hat = mm(sp.invU, sp.tau_hat);   // column-independent factor of dS

    mat DSDE(6, 6);
    double y[6];
    for (int j = 0; j < 6; ++j) {
        const m33 dC = 2. * unit_strain_tensor(j);     // engineering Voigt basis
        const m33 dCp = mm(mtm(sp.N, dC), sp.N);       // shared eigenbasis transform
        const m33 dh = mmt(mm(sp.N, m33(sp.G_h % dCp)), sp.N);   // d(1/2 ln C) = d(ln U)
        // The box acts in the spatial frame (d ln V = R d(ln U) R^T at fixed R):
        // conjugate in, apply Lt, conjugate back to the corotational frame.
        box_apply(Lt_, mmt(mm(sp.R, dh), sp.R), y);
        const m33 dtau_hat = mm(mtm(sp.R, v2t_stress_arr(y)), sp.R);
        const m33 dinvU = mmt(mm(sp.N, m33(sp.G_iU % dCp)), sp.N);   // d(C^{-1/2})
        const m33 dS = mm(mm(sp.invU, dtau_hat), sp.invU)
                     + mm(mm(dinvU, sp.tau_hat), sp.invU) + mm(invU_tau_hat, dinvU);
        set_stress_col(DSDE, j, dS);
    }
    return DSDE;
}

mat Dtau_LieDD_2_DSDE(const mat &Lt, const mat &F){
    
    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Dtau_LieDD_2_DSDE.");
    }

    return push_forward_voigt(mat::fixed<6,6>(Lt), mat::fixed<3,3>(invF));   // pull-back by F^-1
}

mat DtauDe_JaumannDD_2_DSDE(const mat &Lt, const mat &F, const mat &tau){
    
    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DtauDe_JaumannDD_2_DSDE.");
    }   
    if (is_exactly_symmetric(tau)) {
        // Dtau_LieDD = Dtau_JaumannDD - Jaumann correction, on the Voigt entries (both minor symmetries
        // hold exactly for a symmetric tau), pulled back by F^-1
        double R[6][3][3];
        voigt_expand(R, mat::fixed<6,6>(Lt));
        add_jaumann_correction(R, fastor_sym_copy(tau), false);
        return push_forward_voigt(voigt_collapse(R), mat::fixed<3,3>(invF));
    }

    auto invF_ = to_fastor2(invF, false);
    auto tau_ = to_fastor2(tau);
    auto Dtau_JaumannDD = to_fastor4(Lt);

    // Dtau_LieDD = Dtau_JaumannDD - Jaumann correction
    auto Dtau_LieDD = apply_jaumann_correction(Dtau_JaumannDD, tau_, false);
    // DSDE = pull-back of Dtau_LieDD
    auto DSDE = pullback_4(Dtau_LieDD, invF_);
    return fastor4_to_voigt(DSDE);
}

mat DsigmaDe_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DsigmaDe_2_DSDE.");
    }
    return DtauDe_2_DSDE(J*Lt, F, Cauchy2Kirchoff(sigma, F, J));
}

mat Dsigma_LieDD_2_DSDE(const mat &Lt, const mat &F){

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Dsigma_LieDD_2_DSDE.");
    }
    return Dtau_LieDD_2_DSDE(J*Lt, F);
}

mat DsigmaDe_JaumannDD_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DsigmaDe_JaumannDD_2_DSDE.");
    }
    return DtauDe_JaumannDD_2_DSDE(J*Lt, F, Cauchy2Kirchoff(sigma, F, J));
}

mat DtauDe_GreenNaghdiDD_2_DSDE(const mat &Lt, const mat &F, const mat &tau){

    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DtauDe_GreenNaghdiDD_2_DSDE.");
    }
    if (is_exactly_symmetric(tau)) {
        // Dtau_LieDD = Dtau_GN + B-correction, on the Voigt entries (both minor symmetries hold exactly
        // for a symmetric tau), pulled back by F^-1
        double R[6][3][3];
        add_B_correction(R, mat::fixed<6,6>(Lt), mat::fixed<6,6>(get_BBBB_GN(F)), fastor_sym_copy(tau), true);
        return push_forward_voigt(voigt_collapse(R), mat::fixed<3,3>(invF));
    }
    auto invF_ = to_fastor2(invF, false);
    auto tau_ = to_fastor2(tau);
    auto Dtau_GN = to_fastor4(Lt);
    auto B_ = to_fastor4(get_BBBB_GN(F));
    auto Dtau_LieDD = apply_B_correction(Dtau_GN, B_, tau_, true);
    auto DSDE = pullback_4(Dtau_LieDD, invF_);
    return fastor4_to_voigt(DSDE);
}

mat DsigmaDe_GreenNaghdiDD_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DsigmaDe_GreenNaghdiDD_2_DSDE.");
    }
    return DtauDe_GreenNaghdiDD_2_DSDE(J*Lt, F, Cauchy2Kirchoff(sigma, F, J));
}

mat DtauDe_2_DsigmaDe(const mat &Lt, const double &J) {
    
    assert(J > simcoon::iota);
    return (1./J)*Lt;
}

mat DsigmaDe_2_DtauDe(const mat &Lt, const double &J) {
    
    assert(J > simcoon::iota);    
    return Lt*J;
}

mat DSDE_2_DtauDe(const mat &DSDE, const mat &F, const mat &tau) {

    dk_spectral sp;
    sp.build(F, tau);
    const m33 G_h_inv = 1.0 / sp.G_h;      // elementwise; G_h > 0 (1/2 ln strictly increasing)
    const m33 G_comp = sp.G_iU % G_h_inv;  // composed d(C^{-1/2}) o [d(1/2 ln C)]^{-1}
    const arma::mat::fixed<6,6> DSDE_(0.5 * DSDE);
    const m33 invU_tau_hat = mm(sp.invU, sp.tau_hat);   // column-independent factor

    mat Lt(6, 6);
    double y[6];
    for (int j = 0; j < 6; ++j) {
        // box argument is ln V: pull dh back by R (inverse of the forward R conjugation)
        const m33 dh = mm(mtm(sp.R, unit_strain_tensor(j)), sp.R);
        const m33 dhp = mm(mtm(sp.N, dh), sp.N);                     // eigenframe
        const m33 dC = mmt(mm(sp.N, m33(G_h_inv % dhp)), sp.N);      // inverse of d(1/2 ln C)
        box_apply(DSDE_, dC, y);
        const m33 dS = v2t_stress_arr(y);
        const m33 dinvU = mmt(mm(sp.N, m33(G_comp % dhp)), sp.N);
        const m33 inner = dS - mm(mm(dinvU, sp.tau_hat), sp.invU) - mm(invU_tau_hat, dinvU);
        const m33 dtau_hat = mm(mm(sp.U, inner), sp.U);
        set_stress_col(Lt, j, mmt(mm(sp.R, dtau_hat), sp.R));
    }
    return Lt;
}

mat DSDE_2_DsigmaDe(const mat &DSDE, const mat &F, const mat &sigma) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_DsigmaDe.");
    }
    return (1./J)*DSDE_2_DtauDe(DSDE, F, Cauchy2Kirchoff(sigma, F, J));
}

mat DSDE_2_Dtau_LieDD(const mat &DSDE, const mat &F) {

    return push_forward_voigt(mat::fixed<6,6>(DSDE), mat::fixed<3,3>(F));
}

mat DSDE_2_Dsigma_LieDD(const mat &DSDE, const mat &F) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_DsigmaDe_LieDD.");
    }
    return (1./J)*push_forward_voigt(mat::fixed<6,6>(DSDE), mat::fixed<3,3>(F));   // DSDE_2_Dtau_LieDD / J
}

mat DSDE_2_Dtau_JaumannDD(const mat &DSDE, const mat &F, const mat &tau) {

    // Push-forward + Jaumann correction, on the Voigt entries only (the same sums as the
    // 3x3x3x3 Fastor route)
    double R[6][3][3];
    push_forward_voigt_core(R, mat::fixed<6,6>(DSDE), mat::fixed<3,3>(F));
    add_jaumann_correction(R, fastor_sym_copy(tau), true);
    return voigt_average(R);
}

mat DSDE_2_Dsigma_JaumannDD(const mat &DSDE, const mat &F, const mat &sigma) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_Dsigma_JaumannDD.");
    }
    return (1./J)*DSDE_2_Dtau_JaumannDD(DSDE, F, Cauchy2Kirchoff(sigma, F, J));
}

mat DSDE_2_Dtau_GreenNaghdiDD(const mat &DSDE, const mat &F, const mat &tau) {

    mat Dtau_LieDD = DSDE_2_Dtau_LieDD(DSDE, F);
    return Dtau_LieDD_Dtau_GreenNaghdiDD(Dtau_LieDD, F, tau);
}

mat DSDE_2_Dsigma_GreenNaghdiDD(const mat &DSDE, const mat &F, const mat &sigma) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_Dsigma_GreenNaghdiDD.");
    }
    return (1./J)*DSDE_2_Dtau_GreenNaghdiDD(DSDE, F, Cauchy2Kirchoff(sigma, F, J));
}

// Standard logarithmic reverse convenience functions
mat DSDE_2_Dtau_logarithmicDD(const mat &DSDE, const mat &F, const mat &tau) {

    // XBM log rate: the in-rate tangent is the box tangent — exact core
    // directly (the Lie round trip would cancel anyway).
    return DSDE_2_DtauDe(DSDE, F, tau);
}

mat DSDE_2_Dsigma_logarithmicDD(const mat &DSDE, const mat &F, const mat &sigma) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_Dsigma_logarithmicDD.");
    }
    return (1./J)*DSDE_2_Dtau_logarithmicDD(DSDE, F, Cauchy2Kirchoff(sigma, F, J));
}

// Corate-dispatched material<->box tangent maps (full doc in objective_rates.hpp): match the spin
// kernel to corate_type so the round-trip dS/dE <-> Lt is exact. DtauDe_corate_2_DSDE is the inverse.
mat DSDE_2_DtauDe_corate(const mat &DSDE, const int &corate_type, const mat &F, const mat &tau) {
    switch (corate_type) {
        case 0:  return DSDE_2_Dtau_JaumannDD(DSDE, F, tau);
        case 1:  return DSDE_2_Dtau_GreenNaghdiDD(DSDE, F, tau);    // plain GN: path-integral strain -> rate identity
        case 4:  return DSDE_2_Dtau_LieDD(DSDE, F);                 // Truesdell: the convected box IS the Lie tangent
        case 5:  return DSDE_2_Dtau_JaumannDD(DSDE, F, tau)*inv(A_F(F));   // log_F: De = A^F:D dt, stress carried by sym(DF X DF^-1)
        case corate_type_total:  return DSDE_2_logU_box(DSDE, F, tau);            // strain space: dT/dH in reference components
        case corate_type_material_tangent:  return DSDE;            // the material tangent itself
        case 2:                                                     // XBM (logarithmic)
        case 3:                                                     // log_R: A^R:D accumulates exactly ln U (a function of C) -> exact map
        default: return DSDE_2_DtauDe(DSDE, F, tau);
    }
}

mat dH_dE(const mat &C) {
    vec c;
    mat N;
    eig_sym(c, N, C);
    // g_ab: the divided difference of 1/2 ln, free of the cancellation of two close logs
    mat g(3, 3);
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            g(a, b) = logU_spectral::dd1(c(a), c(b));
    mat P(6, 6);
    for (int j = 0; j < 6; j++) {
        vec e = zeros(6);
        e(j) = 1.;
        P.col(j) = t2v_strain(spectral_congruence(N, 2.*g, v2t_strain(e)));
    }
    return P;
}

vec box_covector_2_dE(const vec &d_box, const mat &F, const int &corate_type) {
    mat P = dH_dE(F.t()*F);
    if (corate_type != corate_type_total) {   // rate form: the box strain is R H R^T
        mat R, U;
        RU_decomposition(R, U, F);
        for (int j = 0; j < 6; j++)
            P.col(j) = t2v_strain(R*v2t_strain(vec(P.col(j)))*R.t());
    }
    return P.t()*d_box;
}

mat DtauDe_corate_2_DSDE(const mat &Lt, const int &corate_type, const mat &F, const mat &tau) {
    switch (corate_type) {
        case 0:  return DtauDe_JaumannDD_2_DSDE(Lt, F, tau);
        case 1:  return DtauDe_GreenNaghdiDD_2_DSDE(Lt, F, tau);    // plain GN: rate identity (see forward map)
        case 4:  return Dtau_LieDD_2_DSDE(Lt, F);                   // Truesdell: convected box = Lie tangent
        case 5:  return DtauDe_JaumannDD_2_DSDE(Lt*A_F(F), F, tau); // log_F: inverse of the forward map
        case corate_type_total:  return logU_box_2_DSDE(Lt, F, tau);              // strain space: P^T Lt P + T:L
        case corate_type_material_tangent:  return Lt;              // already dS/dE
        case 2:                                                     // XBM and...
        case 3:                                                     // log_R: exact map (see DSDE_2_DtauDe_corate)
        default: return DtauDe_2_DSDE(Lt, F, tau);
    }
}

mat Dtau_LieDD_2_DtauDe_corate(const mat &Dtau_LieDD, const int &corate_type, const mat &F, const mat &tau) {
    switch (corate_type) {
        case 0:  return Dtau_LieDD_Dtau_JaumannDD(Dtau_LieDD, tau);
        case 1:  return Dtau_LieDD_Dtau_GreenNaghdiDD(Dtau_LieDD, F, tau);
        case 4:  return Dtau_LieDD;                                   // Truesdell: the convected box IS the Lie tangent
        case 5:  return Dtau_LieDD_Dtau_JaumannDD(Dtau_LieDD, tau)*inv(A_F(F));   // log_F: chain rule through De = A^F:D dt
        case corate_type_total:  return DSDE_2_logU_box(Dtau_LieDD_2_DSDE(Dtau_LieDD, F), F, tau);   // strain space, through dS/dE
        case corate_type_material_tangent:  return Dtau_LieDD_2_DSDE(Dtau_LieDD, F);   // the pull-back alone
        case 2:                                                       // XBM and...
        case 3:                                                       // log_R: exact spectral map
        default: return Dtau_LieDD_Dtau_logarithmicDD(Dtau_LieDD, F, tau);
    }
}

mat Dtau_LieDD_2_DtauDe_stretch(const mat &Dtau_LieDD, const vec::fixed<3> &e, const mat::fixed<3,3> &N, const mat &tau) {
    m33 H;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            const double d = e(a) - e(b);
            H(a, b) = (std::fabs(d) > 1.e-4) ? std::tanh(d)/d : 1. - d*d/3. + 2.*d*d*d*d/15.;   // Taylor of tanh(d)/d
        }
    }
    const arma::mat::fixed<6,6> Jc(Dtau_LieDD);
    const m33 tau_(tau);
    mat Lt(6, 6);
    double y[6];
    for (int j = 0; j < 6; ++j) {
        const m33 dt = spectral_congruence(N, H, unit_strain_tensor(j));
        box_apply(Jc, dt, y);
        const m33 dtau = v2t_stress_arr(y) + mm(dt, tau_) + mm(tau_, dt);
        set_stress_col(Lt, j, dtau);
    }
    return Lt;
}

// Canonical box tangent Lt = d(tau_hat)/d(De) (full doc in objective_rates.hpp).
// From the material tangent dS/dE:
mat box_DtauDe_from_dSdE(const mat &dSdE, const mat &F, const vec &sigma) {
    return DSDE_2_DtauDe(dSdE, F, det(F)*v2t_stress(sigma));
}

mat Dtau_LieDD_Dtau_JaumannDD(const mat &Dtau_LieDD, const mat &tau) {

    // the correction on the Voigt entries only (the same sums as the 3x3x3x3 Fastor route)
    double R[6][3][3];
    voigt_expand(R, mat::fixed<6,6>(Dtau_LieDD));
    add_jaumann_correction(R, fastor_sym_copy(tau), true);
    return voigt_average(R);
}

mat Dtau_LieDD_Dtau_objectiveDD(const mat &Dtau_LieDD, const mat &B, const mat &tau) {

    // Dtau_objectiveDD = Dtau_LieDD - B-correction, on the Voigt entries only (the same sums as
    // the 3x3x3x3 Fastor route)
    double R[6][3][3];
    add_B_correction(R, mat::fixed<6,6>(Dtau_LieDD), mat::fixed<6,6>(B), fastor_sym_copy(tau), false);
    return voigt_average(R);
}

mat Dtau_LieDD_Dtau_GreenNaghdiDD(const mat &Dtau_LieDD, const mat &F, const mat &tau) {

    mat B = get_BBBB_GN(F);
    return Dtau_LieDD_Dtau_objectiveDD(Dtau_LieDD, B, tau);
}

mat Dtau_LieDD_Dtau_logarithmicDD(const mat &Dtau_LieDD, const mat &F, const mat &tau) {

    // XBM in-rate tangent = box tangent: route through the exact spectral pair
    return DSDE_2_DtauDe(Dtau_LieDD_2_DSDE(Dtau_LieDD, F), F, tau);
}

mat Dsigma_LieDD_Dsigma_JaumannDD(const mat &Dsigma_LieDD, const mat &sigma) {

    // as Dtau_LieDD_Dtau_JaumannDD
    double R[6][3][3];
    voigt_expand(R, mat::fixed<6,6>(Dsigma_LieDD));
    add_jaumann_correction(R, fastor_sym_copy(sigma), true);
    return voigt_average(R);
}

mat Dsigma_LieDD_Dsigma_objectiveDD(const mat &Dsigma_LieDD, const mat &B, const mat &sigma) {

    // Dsigma_objectiveDD = Dsigma_LieDD - B-correction, as Dtau_LieDD_Dtau_objectiveDD
    double R[6][3][3];
    add_B_correction(R, mat::fixed<6,6>(Dsigma_LieDD), mat::fixed<6,6>(B), fastor_sym_copy(sigma), false);
    return voigt_average(R);
}

mat Dsigma_LieDD_Dsigma_GreenNaghdiDD(const mat &Dsigma_LieDD, const mat &F, const mat &sigma) {

    mat B = get_BBBB_GN(F);
    return Dsigma_LieDD_Dsigma_objectiveDD(Dsigma_LieDD, B, sigma);
}

mat Dsigma_LieDD_Dsigma_logarithmicDD(const mat &Dsigma_LieDD, const mat &F, const mat &sigma) {

    // Same exact corate-2 core as Dtau_LieDD_Dtau_logarithmicDD, on the
    // Cauchy (1/J) menu.
    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside Dsigma_LieDD_Dsigma_logarithmicDD.");
    }
    return (1./J)*DSDE_2_DtauDe(Dtau_LieDD_2_DSDE(J*Dsigma_LieDD, F), F, Cauchy2Kirchoff(sigma, F, J));
}

mat DSDE_DBiotStressDU(const mat &DSDE, const mat &U, const mat &S) {

    // C_sjpr = 0.5*δ_is*(U_lj*U_mp*δ_nr*DSDE_ilmn + U_lj*δ_mr*U_np*DSDE_ilmn)
    //        + 0.5*S_sp*δ_rj + 0.5*S_rj*δ_sp
    // Straight from the Voigt DSDE, for the Voigt pairs s <= j only, the (l, m, n) terms the deltas
    // keep summed in their original order (the same sums as the 3x3x3x3 contraction)
    const mat::fixed<6,6> D(DSDE);
    // U and S as the Fastor copies read them (a plain memcpy of the column-major storage)
    const mat::fixed<3,3> U_t = mat::fixed<3,3>(U).t(), S_t = mat::fixed<3,3>(S).t();
    double C[6][3][3];
    for (int I = 0; I < 6; ++I) {
        const int s = vi[I], j = vj[I];
        for (int p = 0; p < 3; ++p) {
            for (int r = 0; r < 3; ++r) {
                double sum = 0.0;
                for (int l = 0; l < 3; ++l) {
                    for (int m = 0; m < 3; ++m) {
                        // term1 (n = r) always; term2 for every n when m = r
                        sum += 0.5 * U_t(l,j) * U_t(m,p) * D(vm[s][l], vm[m][r]);
                        if (m != r) continue;
                        for (int n = 0; n < 3; ++n) {
                            sum += 0.5 * U_t(l,j) * U_t(n,p) * D(vm[s][l], vm[m][n]);
                        }
                    }
                }
                // Jaumann-like terms
                if (r == j) sum += 0.5 * S_t(s,p);
                if (s == p) sum += 0.5 * S_t(r,j);
                C[I][p][r] = sum;
            }
        }
    }
    mat L(6,6);
    for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J)
            L(I, J) = 0.5 * (C[I][vi[J]][vj[J]] + C[I][vj[J]][vi[J]]);
    return L;

}

} //namespace simcoon
