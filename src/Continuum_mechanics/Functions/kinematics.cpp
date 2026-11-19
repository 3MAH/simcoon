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
#include <simcoon/Simulation/Maths/eig_sym3.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

mat ER_to_F(const mat &E, const mat &R) {

    assert(E.n_cols == 3);
    assert(E.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //C = 2E + I; symmatu() so the symmetry check of sqrtmat_sympd accepts it
    mat C = symmatu(2.*E+eye(3,3));

    mat::fixed<3,3> U;
    if (!sqrtmat_sympd_3x3(U, C)) {
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside ER_to_F.");
    }
    return R*U;
}

mat eR_to_F(const mat &e, const mat &R) {

    assert(e.n_cols == 3);
    assert(e.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //From e we compute V : e = 1/2 (C-I) --> C = U^2 = 2E+I
    mat::fixed<3,3> V;
    if (!expmat_sym_3x3(V, e)) {
        throw simcoon::exception_expmat_sym("Error in expmat_sym function inside eR_to_F.");
    }
    return (V*R);
}
    
mat G_UdX(const mat &F) {
    return F - eye(3,3);
}

mat G_Udx(const mat &F) {

    try {
        return eye(3,3) - inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv : " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside G_Udx.");
    }        
}
    
mat R_Cauchy_Green(const mat &F) {
    return F.t()*F;
}
    
mat L_Cauchy_Green(const mat &F) {
    return F*F.t();
}

//Right polar decomposition F = R U from one eigen-solve of C = F^T F = N diag(c) N^T
void RU_spectral(mat::fixed<3,3> &R, vec::fixed<3> &c, mat::fixed<3,3> &N, const mat &F) {
    const mat::fixed<3,3> F_(F);
    // same failure modes as sqrtmat_sympd(C) then inv(U): a negative (or non-finite) spectrum,
    // then a singular stretch
    if (!eig_sym_3x3(c, N, F_.t()*F_) || c(0) < 0.) {
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside RU_decomposition.");
    }
    if (!(c(0) > 0.)) {
        throw simcoon::exception_inv("Error in inv function inside RU_decomposition.");
    }
    vec::fixed<3> inv_s;
    for (int k = 0; k < 3; k++) inv_s(k) = 1./std::sqrt(c(k));
    R = F_*spectral_compose_3x3(N, inv_s);
}

void RU_decomposition(mat &R, mat &U, const mat &F) {
    mat::fixed<3,3> R_, N;
    vec::fixed<3> c;
    RU_spectral(R_, c, N, F);
    R = R_;
    U = spectral_compose_3x3(N, sqrt(c));
}

void VR_decomposition(mat &V, mat &R, const mat &F) {
    const mat::fixed<3,3> F_(F);
    vec::fixed<3> b;
    mat::fixed<3,3> N;
    if (!eig_sym_3x3(b, N, F_*F_.t()) || b(0) < 0.) {
        throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside VR_decomposition.");
    }
    if (!(b(0) > 0.)) {
        throw simcoon::exception_inv("Error in inv function inside VR_decomposition.");
    }
    vec::fixed<3> inv_s;
    for (int k = 0; k < 3; k++) inv_s(k) = 1./std::sqrt(b(k));
    V = spectral_compose_3x3(N, sqrt(b));
    R = spectral_compose_3x3(N, inv_s)*F_;
}
    
vec Inv_X(const mat &X) {
    
    vec I = zeros(3);
    I(0) = trace(X);
    I(1) = 0.5*(pow(trace(X),2.) - trace(X*X));
    I(2) = det(X);    
    return I;
}

mat Cauchy(const mat &F) {
    try {
        return inv(L_Cauchy_Green(F));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Cauchy.");
    }     
}
    
mat Green_Lagrange(const mat &F) {
    return 0.5*(R_Cauchy_Green(F) - eye(3,3));
}

mat Euler_Almansi(const mat &F) {
    return 0.5*(eye(3,3) - Cauchy(F));
}

mat Log_strain(const mat &F) {
    const mat::fixed<3,3> F_(F);
    vec::fixed<3> b;
    mat::fixed<3,3> N;
    if (!eig_sym_3x3(b, N, F_*F_.t()) || !(b(0) > 0.)) {
        throw simcoon::exception_logmat_sympd("Error in logmat_sympd function inside Log_strain.");
    }
    for (int k = 0; k < 3; k++) b(k) = 0.5*std::log(b(k));
    return spectral_compose_3x3(N, b);
}

mat finite_L(const mat &F0, const mat &F1, const double &DTime) {
    if (DTime <= simcoon::iota) {
        return zeros(3,3);   // a rate is undefined at DTime = 0 (see finite_rotation)
    }
    
    //Definition of L = dot(F)*F^-1
    try {
        return (1./DTime)*(F1-F0)*inv(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_L.");
    }   
}

mat finite_D(const mat &F0, const mat &F1, const double &DTime) {
    const mat L = finite_L(F0, F1, DTime);
    return 0.5*(L+L.t());
}

mat finite_W(const mat &F0, const mat &F1, const double &DTime) {

    const mat L = finite_L(F0, F1, DTime);
    return 0.5*(L-L.t());

}
    
void finite_rotation(const mat &F0, const mat &F1, const double &DTime, mat &DR, mat &Omega) {
    mat::fixed<3,3> R0, R1, N;
    vec::fixed<3> c;
    RU_spectral(R0, c, N, F0);
    RU_spectral(R1, c, N, F1);
    // Exact and DTime-independent: DR stays the true relative rotation even
    // when the time increment degenerates (the spin, a rate, does not).
    DR = R1*R0.t();

    if (DTime <= simcoon::iota) {
        Omega = zeros(3,3);
        return;
    }
    try {
        // singular at a 180-degree increment or det F < 0 -> exception_inv (solver step-cut signal)
        Omega = (2./DTime)*(DR-eye(3,3))*inv(DR+eye(3,3));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_rotation.");
    }
}

mat finite_Omega(const mat &F0, const mat &F1, const double &DTime) {
    mat DR, Omega;
    finite_rotation(F0, F1, DTime, DR, Omega);
    return Omega;
}

mat finite_DQ(const mat &Omega0, const mat &Omega1, const double &DTime) {
    
    try {
        return (eye(3,3)+0.5*DTime*Omega0)*(inv(eye(3,3)-0.5*DTime*Omega1));
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside finite_DQ.");
    }   
    
}
    
    
} //namespace simcoon
