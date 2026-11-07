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

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <stdexcept>

using namespace std;
using namespace arma;

namespace simcoon{

namespace {

// Fixed-size kernels of the isochoric-invariant functions: the same products summed in the same
// order as the armadillo expressions they replace, temporaries on the stack (bit-identical under
// -ffp-contract=on with Armadillo 15.6; the tolerance-based tests remain the guarantee).
typedef mat::fixed<3,3> m33;
typedef mat::fixed<6,6> m66;

// Voigt index of the pair (i,j)
static const int voigt_id[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};

// trace(X), summed as armadillo's trace (15.6.0): two accumulators, (X00 + X22) + X11
inline double trace3(const m33 &X) {
    double val1 = 0.;
    double val2 = 0.;
    val1 += X(0,0);
    val2 += X(1,1);
    val1 += X(2,2);
    return val1 + val2;
}

// dev(X) = X - ((1/3) tr X) I
inline m33 dev3(const m33 &X) {
    const double s = (1./3.)*trace3(X);
    m33 D;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            const double e = (i == j) ? 1. : 0.;
            const double sph = e*s;
            D(i,j) = X(i,j) - sph;
        }
    }
    return D;
}

// powmat(X, 2): X*X as armadillo's tiny-matrix kernel (15.6.0); a diagonal X takes armadillo's own
// shortcut (a pow call per entry, not always the rounded square)
inline m33 square3(const m33 &X) {
    m33 P;
    if (X.is_diagmat()) {
        P = powmat(X, 2);
        return P;
    }
    for (int j = 0; j < 3; ++j) {
        const double x0 = X(0,j), x1 = X(1,j), x2 = X(2,j);
        for (int i = 0; i < 3; ++i) {
            P(i,j) = X(i,0)*x0 + X(i,1)*x1 + X(i,2)*x2;
        }
    }
    return P;
}

// t2v_sym(A)
inline void sym_voigt(const m33 &A, double v[6]) {
    v[0] = A(0,0);
    v[1] = A(1,1);
    v[2] = A(2,2);
    v[3] = 0.5*(A(0,1) + A(1,0));
    v[4] = 0.5*(A(0,2) + A(2,0));
    v[5] = 0.5*(A(1,2) + A(2,1));
}

// t2v_sym of the identity
const double I_voigt[6] = {1., 1., 1., 0., 0., 0.};

// Idev(), Ivol() and Ireal() entries (constitutive.cpp)
inline double Idev_at(const int i, const int j) {
    if (i < 3 && j < 3) return (i == j) ? (1.-1./3.) : (0.-1./3.);
    return (i == j) ? 0.5 : 0.;
}
inline double Ivol_at(const int i, const int j) {
    return (i < 3 && j < 3) ? (1./3.) : 0.;
}
inline double Ireal_at(const int i, const int j) {
    return (i == j) ? ((i < 3) ? 1. : 0.5) : 0.;
}

// b_bar = J^(-2/3) b
inline m33 b_bar_of(const mat &b, const double &J) {
    const double k = pow(J,-2./3.);
    m33 b_bar;
    for (int e = 0; e < 9; ++e) b_bar(e) = b(e)*k;
    return b_bar;
}

// The first two isochoric invariants of b_bar, as isochoric_invariants
inline void bar_invariants(const m33 &b_bar, double &I1, double &I2) {
    I1 = trace3(b_bar);
    I2 = 0.5*(pow(trace3(b_bar),2.) - trace3(square3(b_bar)));
}

double J_of_b(const mat &b, const double &mJ, const char *where) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det(std::string("Error in det function inside ") + where + ".");
        }
    }
    return J;
}

// 2 dW/dI1 dev(b_bar) + 2 dW/dI2 (tr(b_bar) dev(b_bar) - dev(b_bar^2)); the I_2 part is skipped
// when its derivative is zero (an exact zero term)
m33 tau_iso_invariants_fixed(const double &dWdI_1_bar, const double &dWdI_2_bar, const m33 &b_bar) {
    const m33 dev_b_bar = dev3(b_bar);
    const double c1 = 2.*dWdI_1_bar;
    m33 tau;
    for (int e = 0; e < 9; ++e) tau(e) = dev_b_bar(e)*c1;
    if (dWdI_2_bar != 0.) {
        const m33 dev_b_bar2 = dev3(square3(b_bar));
        const double tr = trace3(b_bar);
        const double c2 = 2.*dWdI_2_bar;
        for (int e = 0; e < 9; ++e) {
            const double p = dev_b_bar(e)*tr;
            const double d = p - dev_b_bar2(e);
            const double q = d*c2;
            tau(e) = tau(e) + q;
        }
    }
    return tau;
}

// gamma_linear: (4/3) (tr_A I_dev - (dev_A (x) I + I (x) dev_A)), scaled by d1, added to Lt
// (Lt = that product when first), with dev_A in Voigt form
void add_gamma_linear(m66 &Lt, const double a[6], const double &tr_A, const double &d1, const bool &first) {
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            const double s1 = a[i]*I_voigt[j];
            const double s2 = I_voigt[i]*a[j];
            const double s = s1 + s2;
            const double x = Idev_at(i,j)*tr_A;
            const double g = (x - s)*(4./3.);
            const double t = g*d1;
            Lt(i,j) = first ? t : Lt(i,j) + t;
        }
    }
}

// Lt += (4 a (x) a) d2, auto_sym_dyadic of the Voigt vector a
void add_gamma_quadratic(m66 &Lt, const double a[6], const double &d2) {
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            const double g = (a[i]*a[j])*4.;
            const double t = g*d2;
            Lt(i,j) = Lt(i,j) + t;
        }
    }
}

// L_iso_hyper_invariants without the 1/J, on the stack
m66 L_iso_invariants_fixed(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const m33 &b_bar) {
    const m33 dev_b_bar = dev3(b_bar);
    double I1, I2;
    bar_invariants(b_bar, I1, I2);
    double a[6];
    sym_voigt(dev_b_bar, a);

    // Only the terms of a non-zero potential derivative are built, summed in the order of the
    // full expression (a skipped term adds an exact zero): a Neo-Hookean potential needs gamma_1
    // alone, and the I_2 terms cost b_bar^2 and five more dyadic products.
    m66 Lt;
    add_gamma_linear(Lt, a, I1, dWdI_1_bar, true);
    // the I_2 work is split around the gamma_11 term to keep the summation order of the full expression
    const bool with_I2 = (dWdI_2_bar != 0.) || (dW2dI_12_bar != 0.) || (dW2dI_22_bar != 0.);
    double dd[6] = {0., 0., 0., 0., 0., 0.};   // Voigt form of I1 dev(b_bar) - dev(b_bar^2)
    if (with_I2) {
        const m33 dev_b_bar2 = dev3(square3(b_bar));
        m33 devdevbb2;
        for (int e = 0; e < 9; ++e) {
            const double p = dev_b_bar(e)*I1;
            devdevbb2(e) = p - dev_b_bar2(e);
        }
        sym_voigt(devdevbb2, dd);
        if (dWdI_2_bar != 0.) {
            double a2[6], bv[6];
            sym_voigt(dev_b_bar2, a2);
            sym_voigt(b_bar, bv);
            m66 H_bar;                        // auto_sym_dyadic_operator(b_bar)
            for (int i=0; i<3; i++) {
                for (int j=i; j<3; j++) {
                    for (int k=0; k<3; k++) {
                        for (int l=k; l<3; l++) {
                            H_bar(voigt_id[i][j],voigt_id[k][l]) = 0.5*(b_bar(i,k)*b_bar(j,l) + b_bar(i,l)*b_bar(j,k));
                        }
                    }
                }
            }
            // gamma_2 = (8/3) (I2 (I_real - 2 I_vol) - I1 (dev(b_bar) (x) I + I (x) dev(b_bar))
            //         + (dev(b_bar^2) (x) I + I (x) dev(b_bar^2))) + 4 (b_bar (x) b_bar - H_bar)
            for (int j = 0; j < 6; ++j) {
                for (int i = 0; i < 6; ++i) {
                    const double v2 = Ivol_at(i,j)*2.;
                    const double r = Ireal_at(i,j) - v2;
                    const double X = r*I2;
                    const double s1 = a[i]*I_voigt[j];
                    const double s2 = I_voigt[i]*a[j];
                    const double Y = (s1 + s2)*I1;
                    const double z1 = a2[i]*I_voigt[j];
                    const double z2 = I_voigt[i]*a2[j];
                    const double Z = z1 + z2;
                    const double T = ((X - Y) + Z)*(8./3.);
                    const double bb = bv[i]*bv[j];
                    const double U = (bb - H_bar(i,j))*4.;
                    const double g = T + U;
                    const double t = g*dWdI_2_bar;
                    Lt(i,j) = Lt(i,j) + t;
                }
            }
        }
    }
    if (dW2dI_11_bar != 0.) {
        add_gamma_quadratic(Lt, a, dW2dI_11_bar);
    }
    if (with_I2) {
        if (dW2dI_12_bar != 0.) {
            for (int j = 0; j < 6; ++j) {
                for (int i = 0; i < 6; ++i) {
                    const double s1 = a[i]*dd[j];
                    const double s2 = dd[i]*a[j];
                    const double g = (s1 + s2)*4.;
                    const double t = g*dW2dI_12_bar;
                    Lt(i,j) = Lt(i,j) + t;
                }
            }
        }
        if (dW2dI_22_bar != 0.) {
            add_gamma_quadratic(Lt, dd, dW2dI_22_bar);
        }
    }
    return Lt;
}

// L_vol_hyper on the stack: (dU/dJ + J d2U/dJ2) 3 I_vol - 2 dU/dJ I_real
m66 L_vol_fixed(const double &dUdJ, const double &dU2dJ2, const double &J) {
    const double k1 = (dUdJ+dU2dJ2*J)*3.;
    const double k2 = 2.*dUdJ;
    m66 L;
    for (int j = 0; j < 6; ++j) {
        for (int i = 0; i < 6; ++i) {
            const double p1 = Ivol_at(i,j)*k1;
            const double p2 = Ireal_at(i,j)*k2;
            L(i,j) = p1 - p2;
        }
    }
    return L;
}

// N_k = n_k n_k^T, the eigenprojectors of the columns of n (n_k(a) n_k(b), as their outer product)
void eigen_projectors(std::vector<mat> &N_projectors, const mat &n) {
    for (int k = 0; k < 3; k++) {
        mat P(3,3);
        for (int q = 0; q < 3; q++)
            for (int p = 0; p < 3; p++)
                P(p,q) = n(p,k)*n(q,k);
        N_projectors[k] = P;
    }
}

}  // namespace

vec isochoric_invariants(const mat &b, const double &mJ) {

    const double J = J_of_b(b, mJ, "isochoric_invariants");
    vec I = zeros(3);
    bar_invariants(b_bar_of(b, J), I(0), I(1));
    I(2) = 1.;
    return I;
}

vec isochoric_invariants(const vec &lambda, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = prod(lambda);
    }
    vec lambda_bar = pow(J,-1./3.)*lambda;
    vec I = zeros(3);    
    I(0) = pow(lambda_bar(0),2.) + pow(lambda_bar(1),2.) + pow(lambda_bar(2),2.);
    I(1) = pow(lambda_bar(0),-2.) + pow(lambda_bar(1),-2.) + pow(lambda_bar(2),-2.);
    I(2) = 1.;
    return I;
}

vec isochoric_pstretch_from_V(const mat &V, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = det(V);
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside isochoric_pstretch_from_V.");
        }    
    }
    vec lambda;
    try {
        lambda = eig_sym(V);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Failed to compute eigenvalues in isochoric_pstretch_from_V.");
    }
    vec lambda_bar = pow(J,-1./3.)*lambda;
    return lambda_bar;    
}

vec isochoric_pstretch_from_b(const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside isochoric_pstretch_from_b.");
        }    
    }
    vec lambda;
    try {
        lambda = eig_sym(b);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Failed to compute eigenvalues in isochoric_pstretch_from_b.");
    }
    lambda.transform( [](double val) { return (sqrt(val)); } );
    vec lambda_bar = pow(J,-1./3.)*lambda;
    return lambda_bar;
}

vec isochoric_pstretch(const mat &input, const string &input_tensor, const double &mJ) {
    if (input_tensor == "b") {
        return isochoric_pstretch_from_b(input, mJ);
    }
    if (input_tensor == "v" || input_tensor == "V") {
        return isochoric_pstretch_from_V(input, mJ);
    }
    throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
}

void pstretch(vec &lambda, mat &n_pvectors, const mat &input, const string &input_tensor, const double &mJ) {

    double J=mJ;
    if (input_tensor == "b") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside pstretch.");
            }    
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside pstretch.");
        }
        lambda.transform( [](double val) { return (sqrt(val)); } );
    }
    else if (input_tensor == "v" || input_tensor == "V") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside pstretch.");
            }   
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);        
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside pstretch.");
        }
    }
    else {
        throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
    }
}

void pstretch(vec &lambda, mat &n_pvectors, std::vector<mat> &N_projectors, const mat &input, const string &input_tensor, const double &mJ) {

    pstretch(lambda, n_pvectors, input, input_tensor, mJ);
    eigen_projectors(N_projectors, n_pvectors);
}

void isochoric_pstretch(vec &lambda_bar, mat &n_pvectors, const mat &input, const string &input_tensor, const double &mJ) {

    double J=mJ;
    vec lambda = zeros(3);

    if (input_tensor == "b") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside isochoric_pstretch.");
            }   
        }        
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside isochoric_pstretch.");
        }        
        lambda.transform( [](double val) { return (sqrt(val)); } );
    }
    else if (input_tensor == "v" || input_tensor == "V") {
        if (fabs(mJ) < simcoon::iota) {
            try {
                J = sqrt(det(input));
            } catch (const std::runtime_error &e) {
                cerr << "Error in det: " << e.what() << endl;
                throw simcoon::exception_det("Error in det function inside isochoric_pstretch.");
            }   
        }
        bool success_eig_sym = eig_sym(lambda, n_pvectors, input);
        if (!success_eig_sym) {
            throw simcoon::exception_eig_sym("Error in eig_sym function inside isochoric_pstretch.");
        }        
    }
    else {
        throw std::invalid_argument("Invalid input string to describe the input vector: it should be *b* for left Cauchy-Green tensor or *v* or *V* for Eulerian stretch tensor");
    }
    lambda_bar = pow(J,-1./3.)*lambda;
}

void isochoric_pstretch(vec &lambda_bar, mat &n_pvectors, std::vector<mat> &N_projectors, const mat &input, const string &input_tensor, const double &mJ) {

    isochoric_pstretch(lambda_bar, n_pvectors, input, input_tensor, mJ);
    eigen_projectors(N_projectors, n_pvectors);
}

vec beta_coefs(const vec &dWdlambda_bar, const vec &lambda_bar) {

    vec beta = lambda_bar%dWdlambda_bar - (1./3.)*sum(lambda_bar%dWdlambda_bar)*ones(3);
    return beta;
}

mat gamma_coefs(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const vec &lambda_bar) {

    mat gamma = zeros(3,3);
    // g_ab = lam_a lam_b W_ab + delta_ab lam_a W_a (W_a term diagonal); gamma = P.g.P = d beta_a / d ln lam_b
    mat factor = dW2dlambda_bar2%(lambda_bar*lambda_bar.t()) + diagmat(dWdlambda_bar%lambda_bar);
    vec factor_sum_col = ((-1./3.)*sum(factor,0)).t();  // sum along dim 0 returns row vec, transpose to col
    vec factor_sum_row = (-1./3.)*sum(factor,1);

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
                gamma(i,j) = factor_sum_col(i) + factor_sum_row(j);
        }            
    }    

    gamma += factor;
    gamma += (1./9.)*accu(factor)*ones(3,3);

    return gamma;
}

mat tau_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &b, const double &mJ) {
    vec lambda_bar = zeros(3);
    mat n_pvectors = zeros(3,3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", mJ);
    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);

    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return tau_iso;
}

mat tau_iso_hyper_pstretch(const vec &dWdlambda_bar, const vec &lambda_bar, const std::vector<mat> &N_projectors) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return tau_iso;
}

mat sigma_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside sigma_iso_hyper_pstretch.");
        }   
    }  
    vec lambda_bar = zeros(3);
    mat n_pvectors = zeros(3,3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", mJ);
    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);

    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return (1./J)*tau_iso;
}

mat sigma_iso_hyper_pstretch(const vec &dWdlambda_bar, const vec &lambda_bar, const std::vector<mat> &N_projectors, const double &J) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat tau_iso = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        tau_iso += beta(i)*N_projectors[i];
    }
    return (1./J)*tau_iso;
}

/*vec a_coefs(const double &dWdI_1_bar, const double &dWdI_2_bar, const vec &I_bar) {
    vec a = zeros(2);
    a(0) = 2.*(dWdI_1_bar + I_bar(0)*dWdI_2_bar);
    a(1) = 2.*(dWdI_2_bar);
    return a;
}

vec b_coefs(const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const vec &I_bar) {
    vec b = zeros(4);
    b(0) = 4.*(dW2dI_11_bar + dWdI_2_bar + pow(I_bar(0),2.)*dW2dI_22_bar + 2*I_bar(0)*dW2dI_12_bar);
    b(1) = 4.*(I_bar(0)*dW2dI_22_bar+dW2dI_12_bar);
    b(2) = 4.*dW2dI_22_bar;
    b(3) = 4.*dWdI_2_bar;
    return b;
}

vec delta_coefs(const vec &a_coefs, const vec &b_coefs, const mat &b) {

    vec I = Inv_X(b);    
    double trb2 = trace(powmat(b,2));
    vec delta = zeros(8);
    delta(0) = b_coefs(0);
    delta(1) = -b_coefs(1);
    delta(2) = (1./3.)*(2*a_coefs(0)-b_coefs(0)*I(0)+b_coefs(1)*trb2);
    delta(3) = b_coefs(2);
    delta(4) = 2.*a_coefs(1)+b_coefs(1)*I(0)-b_coefs(2)*trb2+b_coefs(3);
    delta(5) = (1./9.)*(2.*a_coefs(0)*I(0)-2.*a_coefs(1)*trb2+b_coefs(0)*pow(I(0),2.)) + 
               (1./9.)*(-2.*b_coefs(1)*I(0)*trb2+b_coefs(2)*pow(trb2,2.)-b_coefs(3)*trb2);
    delta(6) = (2./3.)*(a_coefs(0)*I(0)-a_coefs(1)*trb2);
    delta(7) = -b_coefs(3);
    return delta;
}*/

mat tau_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const mat &b, const double &mJ) {
    const double J = J_of_b(b, mJ, "tau_iso_hyper_invariants");
    return tau_iso_invariants_fixed(dWdI_1_bar, dWdI_2_bar, b_bar_of(b, J));
}

mat tau_vol_hyper(const double &dUdJ, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside tau_vol_hyper.");
        }   
    }    
    mat Id = eye(3,3);

    return J*dUdJ*eye(3,3);
}

mat sigma_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const mat &b, const double &mJ) {
    const double J = J_of_b(b, mJ, "sigma_iso_hyper_invariants");
    const m33 tau = tau_iso_invariants_fixed(dWdI_1_bar, dWdI_2_bar, b_bar_of(b, J));
    const double invJ = 1./J;
    mat sigma(3,3);
    for (int e = 0; e < 9; ++e) sigma(e) = tau(e)*invJ;
    return sigma;
}

mat sigma_vol_hyper(const double &dUdJ, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside sigma_vol_hyper.");
        }   
    }    
    mat Id = eye(3,3);
    return dUdJ*eye(3,3);
}

/*mat L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = sqrt(det(b));
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    vec I_bar = isochoric_invariants(b,J);
    mat Id = eye(3,3);

    mat H_1 = ;
    mat H_2 = ;

    mat gamma_1 = 2.*H_1 - (4./3.)*(I_bar(0)*Idev() - (sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar)));
    mat gamma_2 = (8./3.)*(I_bar(0)*(Ireal() - 2.*Ivol()) - I_bar(0)*(sym_dyadic(dev_b_bar,Id)+sym_dyadic(Id,dev_b_bar))
                    + (sym_dyadic(dev_b_bar2,Id)+sym_dyadic(Id,dev_b_bar2))) + 4*(auto_sym_dyadic(dev_b_bar)-H_bar);
    mat gamma_11 = 4.*auto_sym_dyadic(dev_b_bar);
    mat gamma_22 = 4.*auto_sym_dyadic(devdevbb2);
    mat gamma_12 = 4.*(sym_dyadic(dev_b_bar,devdevbb2)+sym_dyadic(devdevbb2,dev_b_bar));

    mat L_iso = zeros (6,6);
    L_iso = (1./J)*(gamma_1*dWdI_1_bar+gamma_2*dWdI_2_bar+gamma_11*dW2dI_11_bar+gamma_12*dW2dI_12_bar+gamma_22*dW2dI_22_bar);

    return L_iso;
}*/

/*mat L_iso_hyper_invariants(const vec &delta_coefs, const mat &b, const double &mJ) {
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        J = sqrt(det(b));
    }    
    mat b_bar = pow(J,-2./3.)*b;
    mat b_bar2 = powmat(b_bar,2);
    mat Id = eye(3,3);

    cout << "delta_coefs(0)*auto_sym_dyadic(b_bar)" << delta_coefs(0)*auto_sym_dyadic(b_bar) << endl;
    cout << "delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar))" << delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar)) << endl;
    cout << "delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar))" << delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar)) << endl;
    cout << "delta_coefs(3)*auto_sym_dyadic(b_bar2)" << delta_coefs(3)*auto_sym_dyadic(b_bar2) << endl;
    cout << "delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2))" << delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2)) << endl;
    cout << "delta_coefs(5)*auto_sym_dyadic(Id)" << delta_coefs(5)*auto_sym_dyadic(Id) << endl;    
    cout << "delta_coefs(6)*auto_sym_dyadic_operator(Id)" << delta_coefs(6)*auto_sym_dyadic_operator(Id) << endl;
    cout << "delta_coefs(7)*auto_sym_dyadic_operator(b_bar)" << delta_coefs(7)*auto_sym_dyadic_operator(b_bar) << endl;                

    mat L_iso = zeros (6,6);
    L_iso = delta_coefs(0)*auto_sym_dyadic(b_bar)+delta_coefs(1)*(sym_dyadic(b_bar,b_bar2)+sym_dyadic(b_bar2,b_bar))+ delta_coefs(2)*(sym_dyadic(b_bar,Id) + sym_dyadic(Id,b_bar)) + delta_coefs(3)*auto_sym_dyadic(b_bar2)
                + delta_coefs(4)*(sym_dyadic(b_bar2,Id) + sym_dyadic(Id,b_bar2)) + delta_coefs(5)*auto_sym_dyadic(Id)
                + delta_coefs(6)*auto_sym_dyadic_operator(Id) + delta_coefs(7)*auto_sym_dyadic_operator(b_bar);
    return L_iso;
}*/

mat L_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const vec &lambda_bar, const mat &n_pvectors, const double &J) {

    vec beta = beta_coefs(dWdlambda_bar, lambda_bar);
    mat gamma = gamma_coefs(dWdlambda_bar, dW2dlambda_bar2, lambda_bar);
    const m33 n(n_pvectors);
    double l2[3];   // lambda_bar^2
    for (int k = 0; k < 3; ++k) l2[k] = pow(lambda_bar(k), 2.);

    // c accumulated on the stack, the dyadic products built in place: the same products, rounded
    // and summed in the same order as dyadic_4vectors_sym and auto_sym_dyadic
    m66 c(fill::zeros);
    double factor = 0.;

    for(unsigned int i=0; i<3; i++) {
        for(unsigned int j=0; j<3; j++) {
            const double delta_ij = (i == j) ? 1. : 0.;
            const double coef = (gamma(i,j) - 2*delta_ij*beta(i));
            for (int p=0; p<3; p++) {
                for (int q=p; q<3; q++) {
                    for (int k=0; k<3; k++) {
                        for (int l=k; l<3; l++) {
                            const double D = n(p,i)*n(q,i)*n(k,j)*n(l,j);   // aabb dyadic of (n_i, n_j)
                            const double t = D*coef;
                            c(voigt_id[p][q],voigt_id[k][l]) = c(voigt_id[p][q],voigt_id[k][l]) + t;
                        }
                    }
                }
            }
        }

        // ordered pairs a!=b: the symmetrized dyadic collapses to s s^T, s = sym dyadic of (n_i n_j^T + n_j n_i^T)
        for(unsigned int j=i+1; j<3; j++) {
            if(l2[j]-l2[i] < 1.E-6) {
                factor = l2[i] * pow(lambda_bar(j), -2.) * (0.5*gamma(j,j) - beta(j)) - 0.5*gamma(i,j);
            }
            else {
                factor = (beta(j)*l2[i] - beta(i)*l2[j])/(l2[j]-l2[i]);
            }

            m33 M;
            for (int q=0; q<3; q++) {
                for (int p=0; p<3; p++) {
                    const double m1 = n(p,i)*n(q,j);
                    const double m2 = n(p,j)*n(q,i);
                    M(p,q) = m1 + m2;
                }
            }
            double s[6];
            sym_voigt(M, s);
            for (int q=0; q<6; q++) {
                for (int p=0; p<6; p++) {
                    const double t = (s[p]*s[q])*factor;
                    c(p,q) = c(p,q) + t;
                }
            }
        }
    }
    const double invJ = 1./J;
    mat L(6,6);
    for (int e = 0; e < 36; ++e) L(e) = c(e)*invJ;
    return L;
}

mat L_iso_hyper_pstretch(const vec &dWdlambda_bar, const mat &dW2dlambda_bar2, const mat &b, const double &mJ) {

    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = sqrt(det(b));
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside L_iso_hyper_pstretch.");
        }
    }

    vec lambda_bar = zeros(3);
    mat n_pvectors = zeros(3,3);
    isochoric_pstretch(lambda_bar, n_pvectors, b, "b", J);
    // the shear coefficient is invariant under the J^(-1/3) bar-scaling, so the
    // lambda_bar-based assembly is exact for b as well
    return L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, lambda_bar, n_pvectors, J);
}

// box operators of any symmetric A with convected rate l A + A l^T - 2/3 tr(d) A (b_bar, fibre
// structure tensors): the anisotropic tangent reuses the I_1 algebra verbatim
static mat gamma_linear(const mat &dev_A, const double &tr_A, const mat &I_dev, const mat &Id) {
    return (4./3.)*(tr_A*I_dev - (sym_dyadic(dev_A,Id)+sym_dyadic(Id,dev_A)));
}

static mat gamma_quadratic(const mat &dev_A) {
    return 4.*auto_sym_dyadic(dev_A);
}

mat L_iso_hyper_invariants(const double &dWdI_1_bar, const double &dWdI_2_bar, const double &dW2dI_11_bar, const double &dW2dI_12_bar, const double &dW2dI_22_bar, const mat &b, const double &mJ) {
    const double J = J_of_b(b, mJ, "L_iso_hyper_invariants");
    const m66 Lt = L_iso_invariants_fixed(dWdI_1_bar, dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar, b_bar_of(b, J));
    const double invJ = 1./J;
    mat L(6,6);
    for (int e = 0; e < 36; ++e) L(e) = Lt(e)*invJ;
    return L;
}

mat L_vol_hyper(const double &dUdJ, const double &dU2dJ2, const mat &b, const double &mJ) {
    return L_vol_fixed(dUdJ, dU2dJ2, J_of_b(b, mJ, "L_vol_hyper"));
}

namespace {

// props(i) is unchecked in release builds: a short props vector would be read
// out of bounds silently.
void require_props(const vec &props, const uword n, const char *name) {
    if (props.n_elem < n) {
        throw std::invalid_argument(std::string("hyper_potential_derivatives: ") + name + " needs "
                                    + std::to_string(n) + " parameters, got "
                                    + std::to_string(props.n_elem));
    }
}

}  // namespace

VolumetricPotential volumetric_potential_of(const vec &props, const uword n_used) {
    if (props.n_elem <= n_used) {
        return VolumetricPotential::LOG_J;
    }
    const double code = props(n_used);
    if (code == 0.) return VolumetricPotential::LOG_J;
    if (code == 1.) return VolumetricPotential::QUADRATIC;
    throw std::invalid_argument("volumetric potential: props(" + std::to_string(n_used) + ") = "
                                + std::to_string(code) + " is neither 0 (kappa (J ln J - J + 1)) nor 1 (kappa/2 (J - 1)^2)");
}

void volumetric_derivatives(const VolumetricPotential &vol, const double &kappa, const double &J, double &dUdJ, double &dU2dJ2) {
    switch (vol) {
        case VolumetricPotential::LOG_J:
            dUdJ = kappa*log(J);
            dU2dJ2 = kappa/J;
            break;
        case VolumetricPotential::QUADRATIC:
            dUdJ = kappa*(J-1.);
            dU2dJ2 = kappa;
            break;
    }
}

hyper_anisotropy hyper_potential_anisotropy(const HyperPotential &potential, const vec &props) {

    hyper_anisotropy an;
    if (potential != HyperPotential::HOLZA) {
        return an;      // isotropic potential: no fibres, no dispersion
    }
    require_props(props, 5, "HOLZA");
    an.kappa_d = props(3);
    if (an.kappa_d < 0. || an.kappa_d > 1./3.) {
        throw std::invalid_argument("HOLZA: the dispersion kappa_d must lie in [0, 1/3], got "
                                    + std::to_string(an.kappa_d));
    }
    if (props(4) < 1.) {
        throw std::invalid_argument("HOLZA: at least one fibre family is required, got "
                                    + std::to_string(props(4)));
    }
    const uword n_fam = uword(props(4));
    require_props(props, 6 + 3*n_fam, "HOLZA");     // the a0 triplets, plus the trailing kappa
    an.a0 = zeros(3, n_fam);
    for (uword i = 0; i < n_fam; i++) {
        vec a = props.subvec(5 + 3*i, 7 + 3*i);
        const double a_norm = norm(a, 2);
        if (a_norm < simcoon::iota) {
            throw std::invalid_argument("HOLZA: fibre direction " + std::to_string(i)
                                        + " has a zero norm");
        }
        an.a0.col(i) = a/a_norm;    // the API already sends unit cosines; normalise anyway
    }
    return an;
}

std::vector<mat> structure_tensors_push_forward(const mat &F, const mat &a0, const double &kappa_d, const double &mJ) {

    std::vector<mat> A;
    if (a0.n_elem == 0) {
        return A;
    }
    double J=mJ;
    if (fabs(mJ) < simcoon::iota) {
        try {
            J = det(F);
        } catch (const std::runtime_error &e) {
            cerr << "Error in det: " << e.what() << endl;
            throw simcoon::exception_det("Error in det function inside structure_tensors_push_forward.");
        }
    }
    const double Jm13 = pow(J,-1./3.);

    // kappa_d = 0 (perfectly aligned fibres, the HGO 2000 model) is the common case and
    // needs no b_bar at all -- building it there is a 3x3 product multiplied by zero.
    mat kd_b_bar;
    if (kappa_d > 0.) {
        kd_b_bar = (kappa_d*Jm13*Jm13)*(F*F.t());
    }

    A.reserve(a0.n_cols);
    for (uword i = 0; i < a0.n_cols; i++) {
        vec a_bar = Jm13*(F*a0.col(i));
        mat A_i = (1.-3.*kappa_d)*(a_bar*a_bar.t());
        if (kappa_d > 0.) {
            A_i += kd_b_bar;
        }
        A.push_back(A_i);
    }
    return A;
}

hyper_invariants_dW hyper_potential_derivatives(const HyperPotential &potential, const vec &props, const vec &I_bar, const double &J, const std::vector<mat> &A) {

    hyper_invariants_dW dW;
    double kappa = 0.;  // the volumetric term U(J) is shared: see volumetric_derivatives
    uword n_used = 0;   // props consumed by the isochoric potential; an optional props(n_used) selects U(J)

    switch (potential) {
        case HyperPotential::NEOHC: {
            // \f$ W = \frac{\mu}{2}*\left(\bar{I}_1 -3 \right) + U(J) \f$
            require_props(props, 2, "NEOHC");
            n_used = 2;
            double mu = props(0);
            kappa = props(1);
            dW.dWdI_1_bar = 0.5*mu;
            break;
        }
        case HyperPotential::MOORI: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{01} left(\bar{I}_2 -3\right) + U(J) \f$
            require_props(props, 3, "MOORI");
            n_used = 3;
            double C_10 = props(0);
            double C_01 = props(1);
            kappa = props(2);
            dW.dWdI_1_bar = C_10;
            dW.dWdI_2_bar = C_01;
            break;
        }
        case HyperPotential::YEOHH: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{20} left(\bar{I}_1 -3\right)^2 + C_{30} left(\bar{I}_1 -3\right)^3 + U(J) \f$
            require_props(props, 4, "YEOHH");
            n_used = 4;
            double C_10 = props(0);
            double C_20 = props(1);
            double C_30 = props(2);
            kappa = props(3);
            dW.dWdI_1_bar = C_10 + 2.*C_20*(I_bar(0)-3.) + 3.*C_30*pow((I_bar(0)-3.),2.);
            dW.dW2dI_11_bar = 2.*C_20 + 6.*C_30*(I_bar(0)-3.);
            break;
        }
        case HyperPotential::ISHAH: {
            // Isihara model (1951)
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{20} left(\bar{I}_1 -3\right)^2 + C_{01} left(\bar{I}_2 -3\right) + U(J) \f$
            require_props(props, 4, "ISHAH");
            n_used = 4;
            double C_10 = props(0);
            double C_20 = props(1);
            double C_01 = props(2);
            kappa = props(3);
            dW.dWdI_1_bar = C_10 + 2.*C_20*(I_bar(0)-3.);
            dW.dW2dI_11_bar = 2.*C_20;
            dW.dWdI_2_bar = C_01;
            break;
        }
        case HyperPotential::GETHH: {
            // Gent-Thomas model (1958)
            // \f$ W = c_1 left(\bar{I}_1 -3\right) + c_2 \textrm{ln} left( \frac{\bar{I}_2}{3}\right) + U(J) \f$
            require_props(props, 3, "GETHH");
            n_used = 3;
            double c_1 = props(0);
            double c_2 = props(1);
            kappa = props(2);
            dW.dWdI_1_bar = c_1;
            if(fabs(I_bar(1)) > simcoon::iota) {
                dW.dWdI_2_bar = c_2/I_bar(1);
                dW.dW2dI_22_bar = -1.*c_2/pow(I_bar(1),2.);
            }
            break;
        }
        case HyperPotential::SWANH: {
            // Swanson model (1985)
            // \f$ W = \frac{3}{2} \sum_{i=1}^n \frac{A_i}{1+\alpha_i} left(\frac{\bar{I}_1}{3}\right)^{1+\alpha_i} + \frac{3}{2} \sum_{i=1}^n \frac{B_i}{1+\beta_i} left(\frac{\bar{I}_2}{3}\right)^{1+\beta_i} + U(J) \f$
            require_props(props, 2, "SWANH");
            int N_Swanson = int(props(0));
            require_props(props, 2 + 4*std::max(N_Swanson, 0), "SWANH");
            n_used = 2 + 4*std::max(N_Swanson, 0);
            kappa = props(1);
            for (int i=0; i<N_Swanson; i++) {
                const double A = props(2+i*4);
                const double B = props(2+i*4+1);
                const double alpha = props(2+i*4+2);
                const double beta = props(2+i*4+3);
                dW.dWdI_1_bar += 1./2.*A*pow((I_bar(0)/3.),alpha);
                dW.dW2dI_11_bar += A*alpha/6.*pow((I_bar(0)/3.),alpha-1.);
                dW.dWdI_2_bar += 1./2.*B*pow((I_bar(1)/3.),beta);
                dW.dW2dI_22_bar += B*beta/6.*pow((I_bar(1)/3.),beta-1.);
            }
            break;
        }
        case HyperPotential::HOLZA: {
            // Gasser-Ogden-Holzapfel (2006); kappa_d = 0 is Holzapfel-Gasser-Ogden (2000)
            // \f$ W = C_{10} \left(\bar{I}_1 - 3\right) + \sum_i \frac{k_1}{2 k_2} \left[ \textrm{exp}\left(k_2 \left(\bar{I}^*_{4,i} - 1\right)^2\right) - 1 \right] + U(J) \f$
            require_props(props, 5, "HOLZA");
            const double C_10 = props(0);
            const double k_1 = props(1);
            const double k_2 = props(2);
            const uword n_fam = uword(std::lround(props(4)));
            require_props(props, 6 + 3*n_fam, "HOLZA");
            n_used = 6 + 3*n_fam;
            kappa = props(5 + 3*n_fam);
            dW.dWdI_1_bar = C_10;
            if (A.size() != n_fam) {
                throw std::invalid_argument("HOLZA: " + std::to_string(A.size())
                                            + " structure tensors for " + std::to_string(n_fam)
                                            + " fibre families (see structure_tensors_push_forward)");
            }
            dW.dWdI_a_bar = zeros(n_fam);
            dW.dW2dI_aa_bar = zeros(n_fam);
            for (uword i = 0; i < n_fam; i++) {
                // tr(A_i) IS the fibre pseudo-invariant: no separate channel for it
                const double E_bar = trace(A[i]) - 1.;
                // switch on I*_4 (not on fibre compression), threshold iota not 0: psi'' jumps 0 -> k1 across it,
                // so round-off at I*_4 = 1 must not select a finite tangent
                if (E_bar <= simcoon::iota) {
                    continue;
                }
                // finite-arithmetic guard: inf here becomes a NaN state the solver cannot step-cut out of
                const double e = exp(std::min(k_2*E_bar*E_bar, 350.));
                dW.dWdI_a_bar(i) = k_1*E_bar*e;
                dW.dW2dI_aa_bar(i) = k_1*(1. + 2.*k_2*E_bar*E_bar)*e;
            }
            break;
        }
        default:
            throw std::invalid_argument("hyper_potential_derivatives: unknown potential "
                                        + std::to_string(static_cast<int>(potential)));
    }

    volumetric_derivatives(volumetric_potential_of(props, n_used), kappa, J, dW.dUdJ, dW.dU2dJ2);
    return dW;
}

void hyper_invariants_response(const hyper_invariants_dW &dW, const mat &b, const double &J, const mat &F, const int &corate_type, vec &tau, mat &Lt_box, const std::vector<mat> &A) {

    if (A.size() != dW.dWdI_a_bar.n_elem || A.size() != dW.dW2dI_aa_bar.n_elem) {
        throw std::invalid_argument("hyper_invariants_response: " + std::to_string(A.size())
                                    + " structure tensors for "
                                    + std::to_string(dW.dWdI_a_bar.n_elem) + " fibre derivatives");
    }

    // Kirchhoff throughout: Cauchy = tau/J only at the output boundary
    const double Jb = J_of_b(b, J, "hyper_invariants_response");
    const m33 b_bar = b_bar_of(b, Jb);
    mat m_tau_iso = tau_iso_invariants_fixed(dW.dWdI_1_bar, dW.dWdI_2_bar, b_bar);
    mat m_tau_vol = tau_vol_hyper(dW.dUdJ, b, J);
    mat m_tau = m_tau_iso + m_tau_vol;

    // builders return c = (1/J) d(L_v tau)/dD (not d(L_v sigma)/dD), with J written asymmetrically
    // (1/J in L_iso, none in L_vol): scale them together as a sum, never one at a time
    const m66 Lt_iso = L_iso_invariants_fixed(dW.dWdI_1_bar, dW.dWdI_2_bar, dW.dW2dI_11_bar, dW.dW2dI_12_bar, dW.dW2dI_22_bar, b_bar);
    const m66 Lt_vol = L_vol_fixed(dW.dUdJ, dW.dU2dJ2, Jb);
    const double invJb = 1./Jb;
    m66 Lt_spatial;
    for (int e = 0; e < 36; ++e) {
        const double l_iso = Lt_iso(e)*invJb;
        Lt_spatial(e) = l_iso + Lt_vol(e);
    }

    // fibre terms: I_1 algebra with b_bar -> structure tensor (kappa_d folded into A); zero below the switch
    if (!A.empty()) {
        mat Id = eye(3,3);
        mat I_dev = Idev();
        for (uword i = 0; i < A.size(); i++) {
            const double d1 = dW.dWdI_a_bar(i);
            const double d2 = dW.dW2dI_aa_bar(i);
            if (d1 == 0. && d2 == 0.) {
                continue;
            }
            const mat dev_A = dev(A[i]);
            m_tau += 2.*d1*dev_A;                                   // Kirchhoff: no 1/J
            Lt_spatial += (1./J)*(gamma_linear(dev_A, trace(A[i]), I_dev, Id)*d1
                                  + gamma_quadratic(dev_A)*d2);     // in c, like the two above
        }
    }
    tau = t2v_stress(m_tau);

    // Box tangent Lt = d(tau_hat)/d(De) in corate_type (Kirchhoff, no J). J*Lt_spatial is the
    // Kirchhoff-Lie tangent d(L_v tau)/dD that the map consumes.
    for (int e = 0; e < 36; ++e) Lt_spatial(e) = Lt_spatial(e)*J;
    Lt_box = Dtau_LieDD_2_DtauDe_corate(Lt_spatial, corate_type, F, m_tau);
}

} //namespace simcoon
