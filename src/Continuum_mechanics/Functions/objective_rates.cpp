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

///@file objective_rate.cpp
///@brief A set of function that help to define different quantities, depending on a selected objective rate
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <FTensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;
using namespace FTensor;

namespace simcoon{

void Jaumann(mat &DR, mat &R, mat &D, mat &W, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    mat L0 = (1./DTime)*(F1-F0)*inv(F0);
    mat L1 = (1./DTime)*(F1-F0)*inv(F1);
    mat L = 0.5*L0+0.5*L1;
    
    //decomposition of L
    D = 0.5*(L+L.t());
    W = 0.5*(L-L.t());
    
    //Jaumann
    DR = (inv(I-0.5*DTime*W))*(I+0.5*DTime*W);
    R += DR;
}
    
void Green_Naghdi(mat &DR, mat &R, mat &D, mat &W, const double &DTime, const mat &F0, const mat &F1) {
    //Green-Naghdi
    mat U0;
    mat R0;
    mat U1;
    mat R1;
    RU_decomposition(R0,U0,F0);
    RU_decomposition(R1,U1,F1);
    
    mat L0 = (1./DTime)*(F1-F0)*inv(F0);
    mat L1 = (1./DTime)*(F1-F0)*inv(F1);
    mat L = 0.5*L0+0.5*L1;
    
    //decomposition of L
    D = 0.5*(L+L.t());
    W = 0.5*(L-L.t());

    //alternative ... to test
    //mat U = 0.5*(U0 + U1);
    //R = 0.5*(R0 + R1);
    //DR = (F1-F0)*inv(U)-R*(U1-U0)*inv(U);
    R = R0;
    DR = (F1-F0)*inv(U1)-R0*(U1-U0)*inv(U1);
}

//This function computes the tangent modulus that links the Lie derivative of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E) and the transformation gradient F
mat DSDE_2_Dtau_LieDD(const mat &DSDE, const mat &F) {

    Tensor2<double,3,3> F_ = mat_FTensor2(F);
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(DSDE);
    Tensor4<double,3,3,3,3> C_;
    
    Index<'i', 3> i;
    Index<'s', 3> s;
    Index<'r', 3> r;
    Index<'p', 3> p;
    
    Index<'L', 3> L;
    Index<'J', 3> J;
    Index<'M', 3> M;
    Index<'N', 3> N;
    
    C_(i,s,r,p) = F_(i,L)*(F_(s,J)*(F_(r,M)*(F_(p,N)*DSDE_(L,J,M,N))));
    return FTensor4_mat(C_);
}
    
//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E), the transformation gradient F and the Kirchoff stress tau
mat DSDE_2_Dtau_JaumannDD(const mat &DSDE, const mat &F, const mat &tau) {
    
    Tensor2<double,3,3> F_ = mat_FTensor2(F);
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(DSDE);
    Tensor4<double,3,3,3,3> C_;
    
    Index<'i', 3> i;
    Index<'s', 3> s;
    Index<'r', 3> r;
    Index<'p', 3> p;
    
    Index<'L', 3> L;
    Index<'J', 3> J;
    Index<'M', 3> M;
    Index<'N', 3> N;
    
    C_(i,s,r,p) = F_(i,L)*(F_(s,J)*(F_(r,M)*(F_(p,N)*DSDE_(L,J,M,N)))) + 0.5*tau_(p,s)*delta_(i,r) + 0.5*tau_(r,s)*delta_(i,p) + 0.5*tau_(i,r)*delta_(s,p) + 0.5*tau_(i,p)*delta_(s,r);
    return FTensor4_mat(C_);
}

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D and the Kirchoff stress tau
mat Dtau_LieDD_Dtau_JaumannDD(const mat &Dtau_LieDD, const mat &tau) {

    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> Dtau_LieDD_ = mat_FTensor4(Dtau_LieDD);
    Tensor4<double,3,3,3,3> Dtau_JaumannDD_;
    
    Index<'i', 3> i;
    Index<'s', 3> s;
    Index<'r', 3> r;
    Index<'p', 3> p;
    
    Dtau_JaumannDD_(i,s,p,r) = Dtau_LieDD_(i,s,p,r) - 0.5*tau_(p,s)*delta_(i,r) - 0.5*tau_(r,s)*delta_(i,p) - 0.5*tau_(i,r)*delta_(s,p) - 0.5*tau_(i,p)*delta_(s,r);
    return FTensor4_mat(Dtau_JaumannDD_);
}

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D and the Kirchoff stress tau
mat Dtau_JaumannDD_Dtau_LieDD(const mat &Dtau_JaumannDD, const mat &tau) {
    
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> Dtau_JaumannDD_ = mat_FTensor4(Dtau_JaumannDD);
    Tensor4<double,3,3,3,3> Dtau_LieDD_;
    
    Index<'i', 3> i;
    Index<'j', 3> s;
    Index<'p', 3> p;
    Index<'r', 3> r;
    
    Dtau_LieDD_(i,s,p,r) = Dtau_JaumannDD_(i,s,p,r) - 0.5*tau_(p,s)*delta_(i,r) - 0.5*tau_(r,s)*delta_(i,p) - 0.5*tau_(i,r)*delta_(s,p) - 0.5*tau_(i,p)*delta_(s,r);
    return FTensor4_mat(Dtau_LieDD_);
}
    
    
//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, weighted by J (the determinant of F), from the tangent modulus that links the Lie derivative of the Kirchoff stress tau to the rate of deformation D
mat Dtau_LieDD_Lt_Aba(const mat &Dtau_LieDD, const mat &tau, const double &J) {

    mat Dtau_JaumannDD = Dtau_LieDD_Dtau_JaumannDD(Dtau_LieDD, tau);
    return (1./J)*Dtau_JaumannDD;
}

//This function computes the tangent modulus that computes the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, weighted by J (the determinant of F), from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D
//Note : This is the tangent modulus utilized by Abaqus
mat Dtau_JaumannDD_Lt_Aba(const mat &Dtau_JaumannDD, const double &J) {
    
    return (1./J)*Dtau_JaumannDD;
}
 
} //namespace simcoon
