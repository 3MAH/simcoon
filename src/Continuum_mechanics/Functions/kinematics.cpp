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

///@file kinematics.cpp
///@brief A set of function that allows various strain transformation (in Finite strains)
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//This function returns F (in a vectorized), from E (Green-Lagrange strain) and R (Rotation matrix), according to a RU decomposition
mat ER_to_F(const mat &E, const mat &R) {

    assert(E.n_cols == 3);
    assert(E.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //From E we compute C : E = 1/2 (C-I) --> C = U^2 = 2E+I
    mat C = 2.*E+eye(3,3);

    vec lambda2_alpha;
    vec lambda_alpha = zeros(3);
    mat N_alpha;

    //Since C=U^2, an eigenvalue decomposition allows to find \lambda_alpha^2 (eigenvalues for U^2), therefore finding \lambda_alpha (eigenvalues for U) is straightforward.
    /*eig_sym(lambda2_alpha, N_alpha, C);
    mat U = zeros(3,3);
    for(unsigned int i=0; i<3; i++) {
        lambda_alpha(i) = sqrt(lambda2_alpha(i));
        vec N = N_alpha.col(i);
        U = U + (lambda_alpha(i)*(N*N.t()));
    }*/
    return (R*sqrtmat_sympd(C));
}

mat eR_to_F(const mat &e, const mat &R) {

    assert(e.n_cols == 3);
    assert(e.n_rows == 3);
    assert(R.n_cols == 3);
    assert(R.n_rows == 3);

    //From e we compute V : e = 1/2 (C-I) --> C = U^2 = 2E+I
    mat V = expmat_sym(e);

    return (V*R);
}
    
//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
mat G_UdX(const mat &F) {
    return F - eye(3,3);
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
mat G_Udx(const mat &F) {
    return eye(3,3) - inv(F);
}
    
//This function computes the Right Cauchy-Green C
mat R_Cauchy_Green(const mat &F) {
    return F.t()*F;
}
    
//This function computes the Left Cauchy-Green B
mat L_Cauchy_Green(const mat &F) {
    return F*F.t();
}

void RU_decomposition(mat &R, mat &U, const mat &F) {
    mat U2 = F.t()*F;
    U = sqrtmat_sympd(U2);
    R = F*inv(U);
}

void VR_decomposition(mat &R, mat &V, const mat &F) {
    mat V2 = F*F.t();
    V = sqrtmat_sympd(V2);
    R = inv(V)*F;
}
    
//This function computes the common Right (or Left) Cauchy-Green invariants
vec Inv_X(const mat &X) {
    
    vec lambda = eig_sym(X);
    vec I = zeros(3);
    I(0) = trace(X);
    I(1) = 0.5*(pow(trace(X),2.) - trace(X*X));
    I(2) = det(X);
    return I;
}

//This function computes the Cauchy deformation tensor b
mat Cauchy(const mat &F) {
    return inv(L_Cauchy_Green(F));
}
    
//This function computes the Green-Lagrange finite strain tensor E
mat Green_Lagrange(const mat &F) {
    return 0.5*(R_Cauchy_Green(F) - eye(3,3));
}

//This function computes the Euler-Almansi finite strain tensor e
mat Euler_Almansi(const mat &F) {
    return 0.5*(eye(3,3) - Cauchy(F));
}

//This function computes the Euler-Almansi finite strain tensor h
mat Log_strain(const mat &F) {
    return  0.5*logmat_sympd(L_Cauchy_Green(F));
}

//This function computes the velocity difference
mat finite_L(const mat &F0, const mat &F1, const double &DTime) {
    
    //Definition of L = dot(F)*F^-1
    return (1./DTime)*(F1-F0)*inv(F1);
}

//This function computes the deformation rate D
mat finite_D(const mat &F0, const mat &F1, const double &DTime) {
    
    //Definition of L = dot(F)*F^-1
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //Definition of the deformation rate D
    return 0.5*(L+L.t());
    
}

//This function computes the spin tensor W (correspond to Jaumann rate)
mat finite_W(const mat &F0, const mat &F1, const double &DTime) {

    //Definition of L = dot(F)*F^-1
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //Definition of the rotation matrix Q
    return 0.5*(L-L.t());

}
    
//This function computes the spin tensor Omega (correspond to Green-Naghdi rate)
// Note : here R is the rigid body rotation in the polar decomposition of the deformation gradient F
mat finite_Omega(const mat &F0, const mat &F1, const double &DTime) {
    
    //Definition of Omega = dot(R)*R^-1 (or R.t() since R is a rotation matrix)
    mat R0 = zeros(3,3);
    mat U0 = zeros(3,3);
    mat R1 = zeros(3,3);
    mat U1 = zeros(3,3);
    RU_decomposition(mat &R0, mat &U0, const mat &F0);
    RU_decomposition(mat &R1, mat &U1, const mat &F1);
    return (1./DTime)*(R1-R0)*R1.t();
}

//This function computes the increment of finite rotation
mat finite_DQ(const mat &Omega0, const mat &Omega1, const double &DTime) {
    
    return (eye(3,3)+0.5*DTime*Omega0)*(inv(eye(3,3)-0.5*DTime*Omega1));
}
    
    
} //namespace simcoon
