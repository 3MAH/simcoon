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
#include <simcoon/FTensor.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>

using namespace std;
using namespace arma;
using namespace FTensor;

namespace simcoon{

void Jaumann(mat &DR, mat &D, mat &W, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //decomposition of L
    D = 0.5*(L+L.t());
    W = 0.5*(L-L.t());
    
    //Jaumann
    DR = (inv(I-0.5*DTime*W))*(I+0.5*DTime*W);
}
    
void Green_Naghdi(mat &DR, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    //Green-Naghdi
    mat I = eye(3,3);
    mat U0;
    mat R0;
    mat U1;
    mat R1;
    RU_decomposition(R0,U0,F0);
    RU_decomposition(R1,U1,F1);
    
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    Omega = (1./DTime)*(R1-R0)*R1.t();

    DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    //alternative ... to test
    //    DR = (F1-F0)*inv(U1)-R0*(U1-U0)*inv(U1);
}

void logarithmic_R(mat &DR, mat &N_1, mat &N_2, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    //Green-Naghdi
    mat I = eye(3,3);
    mat U0;
    mat R0;
    mat U1;
    mat R1;
    RU_decomposition(R0,U0,F0);
    RU_decomposition(R1,U1,F1);
    
    mat L = (1./DTime)*(F1-F0)*inv(F1);
    
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    Omega = (1./DTime)*(R1-R0)*R1.t();

    DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    //alternative ... to test
    //    DR = (F1-F0)*inv(U1)-R0*(U1-U0)*inv(U1);
    
    //Logarithmic
    mat B = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    bool success_eig_sym = eig_sym(bi, Bi, B);
    if (!success_eig_sym) {
        throw simcoon::exception_eig_sym("Error in eig_sym function inside logarithmic_R.");
    }

    std::vector<mat> Bi_proj(3);
    Bi_proj[0] = Bi.col(0)*(Bi.col(0)).t();
    Bi_proj[1] = Bi.col(1)*(Bi.col(1)).t();
    Bi_proj[2] = Bi.col(2)*(Bi.col(2)).t();
    
    N_1 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                N_1+=((1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j)))*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }

    N_2 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                N_2+=((1.-(pow(bi(i)/bi(j),0.5)))/(1.+(pow(bi(i)/bi(j),0.5))))*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }
}

void logarithmic_F(mat &DF, mat &N_1, mat &N_2, mat &D, mat &L, const double &DTime, const mat &F0, const mat &F1) {
    //Green-Naghdi
    mat I = eye(3,3);
    mat U0;
    mat R0;
    mat U1;
    mat R1;
    RU_decomposition(R0,U0,F0);
    RU_decomposition(R1,U1,F1);
    
    L = (1./DTime)*(F1-F0)*inv(F1);
    
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());

    //alternative ... to test
    
    //Logarithmic
    mat B = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    bool success_eig_sym = eig_sym(bi, Bi, B);
    if (!success_eig_sym) {
        throw simcoon::exception_eig_sym("Error in eig_sym function inside logarithmic_R.");
    }
    std::vector<mat> Bi_proj(3);
    Bi_proj[0] = Bi.col(0)*(Bi.col(0)).t();
    Bi_proj[1] = Bi.col(1)*(Bi.col(1)).t();
    Bi_proj[2] = Bi.col(2)*(Bi.col(2)).t();
    
    N_1 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                N_1+=((1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j)))*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }

    N_2 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                N_2+=((1.-(pow(bi(i)/bi(j),0.5)))/(1.+(pow(bi(i)/bi(j),0.5))))*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }
    
    DF = (inv(I-0.5*DTime*L))*(I+0.5*DTime*L);
    
}

void Truesdell(mat &DF, mat &D, mat &L, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    L = (1./DTime)*(F1-F0)*inv(F1);
    //Note that The "spin" is actually L (spin for rigid frames of reference, "flot" for Truesdell)    
    D = 0.5*(L+L.t());
    
    //Truesdell
    DF = (inv(I-0.5*DTime*L))*(I+0.5*DTime*L);
}

mat get_BBBB(const mat &F1) {
    mat B = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    bool success_eig_sym = eig_sym(bi, Bi, B);
    if (!success_eig_sym) {
        throw simcoon::exception_eig_sym("Error in eig_sym function inside logarithmic_R.");
    }
    mat BBBB = zeros(6,6);
    
    double f_z = 0.;
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                BBBB = BBBB + f_z*B_klmn(Bi.col(i),Bi.col(j));
            }
        }
    }
    return BBBB;
}

mat get_BBBB_GN(const mat &F1) {
    mat B = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    bool success_eig_sym = eig_sym(bi, Bi, B);
    if (!success_eig_sym) {
        throw simcoon::exception_eig_sym("Error in eig_sym function inside logarithmic_R.");
    }
    mat BBBB = zeros(6,6);
    
    double f_z = 0.;
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                f_z = (sqrt(bi(j)) - sqrt(bi(i)))/(sqrt(bi(j)) + sqrt(bi(i)));
                BBBB = BBBB + f_z*B_klmn(Bi.col(i),Bi.col(j));
            }
        }
    }
    return BBBB;
}

void logarithmic(mat &DR, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    mat L = zeros(3,3);
    
    if(DTime > simcoon::iota) {
        L = (1./DTime)*(F1-F0)*inv(F1);
    }
        
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    
    //Logarithmic
    mat B = L_Cauchy_Green(F1);
    
    vec bi = zeros(3);
    mat Bi;
    bool success_eig_sym = eig_sym(bi, Bi, B);
    if (!success_eig_sym) {
        throw simcoon::exception_eig_sym("Error in eig_sym function inside logarithmic_R.");
    }
    std::vector<mat> Bi_proj(3);
    Bi_proj[0] = Bi.col(0)*(Bi.col(0)).t();
    Bi_proj[1] = Bi.col(1)*(Bi.col(1)).t();
    Bi_proj[2] = Bi.col(2)*(Bi.col(2)).t();
    
    mat N = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if ((i!=j)&&(fabs(bi(i)-bi(j))>simcoon::iota)) {
                N+=((1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j)))*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }
    Omega = W + N;
    DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
}

mat Delta_log_strain(const mat &D, const mat &Omega, const double &DTime) {
    mat I = eye(3,3);
    mat DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    return 0.5*(D+(DR*D*DR.t()))*DTime;
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DtauDe_2_DSDE(const mat &Lt, const mat &B, const mat &F, const mat &tau){
    
    mat invF = inv(F);
    Tensor2<double,3,3> invF_ = mat_FTensor2(invF);
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> Dtau_logarithmicDD_ = mat_FTensor4(Lt);
    Tensor4<double,3,3,3,3> Dtau_LieDD_ = mat_FTensor4(zeros(6,6));
    Tensor4<double,3,3,3,3> B_ = mat_FTensor4(B);
    Tensor4<double,3,3,3,3> I_ = mat_FTensor4(zeros(6,6));
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(zeros(6,6));
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'p', 3> p;
    
    Index<'L', 3> L;
    Index<'J', 3> H;
    Index<'M', 3> M;
    Index<'N', 3> N;
    
    I_(i,j,k,l) = 0.5*delta_(i,k)*delta_(j,l) + 0.5*delta_(i,l)*delta_(j,k);
    Dtau_LieDD_(i,j,k,l) = Dtau_logarithmicDD_(i,j,k,l) + (B_(i,p,k,l)-I_(i,p,k,l))*tau_(p,j) + tau_(i,p)*(B_(j,p,k,l)-I_(j,p,k,l));
    DSDE_(L,H,M,N) = invF_(l,N)*(invF_(k,M)*(invF_(j,H)*(invF_(i,L)*Dtau_LieDD_(i,j,k,l))));
    return FTensor4_mat(DSDE_);
}

mat Dtau_LieDD_2_DSDE(const mat &Lt, const mat &F){
    
    mat invF = inv(F);
    Tensor2<double,3,3> invF_ = mat_FTensor2(invF);
    Tensor4<double,3,3,3,3> Dtau_LieDD_ = mat_FTensor4(Lt);
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(zeros(6,6));
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    
    Index<'L', 3> L;
    Index<'J', 3> H;
    Index<'M', 3> M;
    Index<'N', 3> N;

    DSDE_(L,H,M,N) = invF_(l,N)*(invF_(k,M)*(invF_(j,H)*(invF_(i,L)*Dtau_LieDD_(i,j,k,l))));
    return FTensor4_mat(DSDE_);
}

mat DtauDe_JaumannDD_2_DSDE(const mat &Lt, const mat &F, const mat &tau){
    
    mat invF = inv(F);
    Tensor2<double,3,3> invF_ = mat_FTensor2(invF);
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> Dtau_JaumannDD_ = mat_FTensor4(Lt);
    Tensor4<double,3,3,3,3> Dtau_LieDD_ = mat_FTensor4(zeros(6,6));
    Tensor4<double,3,3,3,3> I_ = mat_FTensor4(zeros(6,6));
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(zeros(6,6));
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'p', 3> p;
    
    Index<'L', 3> L;
    Index<'J', 3> H;
    Index<'M', 3> M;
    Index<'N', 3> N;
    
    Dtau_LieDD_(i,j,k,l) = Dtau_JaumannDD_(i,j,k,l) - 0.5*tau_(k,j)*delta_(i,l) - 0.5*tau_(l,j)*delta_(i,k) - 0.5*tau_(i,l)*delta_(j,k) - 0.5*tau_(i,k)*delta_(j,l);
    DSDE_(L,H,M,N) = invF_(l,N)*(invF_(k,M)*(invF_(j,H)*(invF_(i,L)*Dtau_LieDD_(i,j,k,l))));
    return FTensor4_mat(DSDE_);
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DsigmaDe_2_DSDE(const mat &Lt, const mat &B, const mat &F, const mat &sigma){
    
    double J = det(F);
    return J*DtauDe_2_DSDE(Lt, B, F, Cauchy2Kirchoff(sigma, F, J));
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DsigmaDe_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){

    double J = det(F);
    mat B = get_BBBB(F);
    return J*DtauDe_2_DSDE(Lt, B, F, Cauchy2Kirchoff(sigma, F, J));
}

mat Dsigma_LieDD_2_DSDE(const mat &Lt, const mat &F){
    
    double J = det(F);
    return J*Dtau_LieDD_2_DSDE(Lt, F);
}

mat DsigmaDe_JaumannDD_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){
    
    double J = det(F);
    return J*DtauDe_JaumannDD_2_DSDE(Lt, F, Cauchy2Kirchoff(sigma, F, J));
}


mat DtauDe_2_DsigmaDe(const mat &Lt, const double &J) {
    
    return (1./J)*Lt;
}

mat DsigmaDe_2_DtauDe(const mat &Lt, const double &J) {
    
    return Lt*J;
}

mat DSDE_2_DtauDe(const mat &DSDE, const mat &B, const mat &F, const mat &tau) {
    
    Tensor2<double,3,3> F_ = mat_FTensor2(F);
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(DSDE);
    Tensor4<double,3,3,3,3> B_ = mat_FTensor4(B);
    Tensor4<double,3,3,3,3> I_;
    Tensor4<double,3,3,3,3> C_;
    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'p', 3> p;
    
    Index<'L', 3> L;
    Index<'J', 3> J;
    Index<'M', 3> M;
    Index<'N', 3> N;
    
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    I_(i,j,k,l) = 0.5*delta_(i,k)*delta_(j,l) + 0.5*delta_(i,l)*delta_(j,k);
    C_(i,j,k,l) = F_(i,L)*(F_(j,J)*(F_(k,M)*(F_(l,N)*DSDE_(L,J,M,N)))) - (B_(i,p,k,l)-I_(i,p,k,l))*tau_(p,j)-tau_(i,p)*(B_(j,p,k,l)-I_(j,p,k,l));
    return FTensor4_mat(C_);
}

mat DSDE_2_DsigmaDe(const mat &DSDE, const mat &B, const mat &F, const mat &sigma) {

    double J = det(F);
    return (1./J)*DSDE_2_DtauDe(DSDE, B, F, Cauchy2Kirchoff(sigma, F, J));
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

mat DSDE_2_DsigmaDe_LieDD(const mat &DSDE, const mat &F) {

    double J = det(F);
    return (1./J)*DSDE_2_Dtau_LieDD(DSDE, F);
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

mat DSDE_2_Dsigma_JaumannDD(const mat &DSDE, const mat &F, const mat &sigma) {
    
    double J = det(F);
    return (1./J)*DSDE_2_Dtau_JaumannDD(DSDE, F, Cauchy2Kirchoff(sigma, F, J));
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
    
    Dtau_JaumannDD_(i,s,p,r) = Dtau_LieDD_(i,s,p,r) + 0.5*tau_(p,s)*delta_(i,r) + 0.5*tau_(r,s)*delta_(i,p) + 0.5*tau_(i,r)*delta_(s,p) + 0.5*tau_(i,p)*delta_(s,r);
    return FTensor4_mat(Dtau_JaumannDD_);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dtau_LieDD_Dtau_objectiveDD(const mat &Dtau_LieDD, const mat &B, const mat &tau) {

    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> tau_ = mat_FTensor2(tau);
    Tensor4<double,3,3,3,3> Dtau_LieDD_ = mat_FTensor4(Dtau_LieDD);
    Tensor4<double,3,3,3,3> B_ = mat_FTensor4(B);
    Tensor4<double,3,3,3,3> I_;
    
    Tensor4<double,3,3,3,3> Dtau_logarithmicDD_;

    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'p', 3> p;
    
    I_(i,j,k,l) = 0.5*delta_(i,k)*delta_(j,l) + 0.5*delta_(i,l)*delta_(j,k);
    Dtau_logarithmicDD_(i,j,k,l) = Dtau_LieDD_(i,j,k,l) - (B_(i,p,k,l)-I_(i,p,k,l))*tau_(p,j)-tau_(i,p)*(B_(j,p,k,l)-I_(j,p,k,l));
    return FTensor4_mat(Dtau_logarithmicDD_);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dtau_LieDD_Dtau_GreenNaghdiDD(const mat &Dtau_LieDD, const mat &F, const mat &tau) {

    mat B = get_BBBB_GN(F);
    return Dtau_LieDD_Dtau_objectiveDD(Dtau_LieDD, B, tau);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dtau_LieDD_Dtau_logarithmicDD(const mat &Dtau_LieDD, const mat &F, const mat &tau) {

    mat B = get_BBBB(F);
    return Dtau_LieDD_Dtau_objectiveDD(Dtau_LieDD, B, tau);
}

//This function computes the tangent modulus that links the Jaumann rate of the Cauchy stress tau to the rate of deformation D, from the tangent modulus that links the Lie derivative of the Cauchy stress tau to the rate of deformation D
mat Dsigma_LieDD_Dsigma_JaumannDD(const mat &Dsigma_LieDD, const mat &sigma) {

    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> sigma_ = mat_FTensor2(sigma);
    Tensor4<double,3,3,3,3> Dsigma_LieDD_ = mat_FTensor4(Dsigma_LieDD);
    Tensor4<double,3,3,3,3> Dsigma_JaumannDD_;
    
    Index<'i', 3> i;
    Index<'s', 3> s;
    Index<'r', 3> r;
    Index<'p', 3> p;
    
    Dsigma_JaumannDD_(i,s,p,r) = Dsigma_LieDD_(i,s,p,r) + 0.5*sigma_(p,s)*delta_(i,r) + 0.5*sigma_(r,s)*delta_(i,p) + 0.5*sigma_(i,r)*delta_(s,p) + 0.5*sigma_(i,p)*delta_(s,r);
    return FTensor4_mat(Dsigma_JaumannDD_);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dsigma_LieDD_Dsigma_objectiveDD(const mat &Dsigma_LieDD, const mat &B, const mat &sigma) {

    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> sigma_ = mat_FTensor2(sigma);
    Tensor4<double,3,3,3,3> Dsigma_LieDD_ = mat_FTensor4(Dsigma_LieDD);
    Tensor4<double,3,3,3,3> B_ = mat_FTensor4(B);
    Tensor4<double,3,3,3,3> I_;
    
    Tensor4<double,3,3,3,3> Dsigma_logarithmicDD_;

    
    Index<'i', 3> i;
    Index<'j', 3> j;
    Index<'k', 3> k;
    Index<'l', 3> l;
    Index<'p', 3> p;
    
    I_(i,j,k,l) = 0.5*delta_(i,k)*delta_(j,l) + 0.5*delta_(i,l)*delta_(j,k);
    Dsigma_logarithmicDD_(i,j,k,l) = Dsigma_LieDD_(i,j,k,l) - (B_(i,p,k,l)-I_(i,p,k,l))*sigma_(p,j)-sigma_(i,p)*(B_(j,p,k,l)-I_(j,p,k,l));
    return FTensor4_mat(Dsigma_logarithmicDD_);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dsigma_LieDD_Dsigma_GreenNaghdiDD(const mat &Dsigma_LieDD, const mat &F, const mat &sigma) {

    mat B = get_BBBB_GN(F);
    return Dsigma_LieDD_Dsigma_objectiveDD(Dsigma_LieDD, B, sigma);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dsigma_LieDD_Dsigma_logarithmicDD(const mat &Dsigma_LieDD, const mat &F, const mat &sigma) {

    mat B = get_BBBB(F);
    return Dsigma_LieDD_Dsigma_objectiveDD(Dsigma_LieDD, B, sigma);
}
 
mat DSDE_DBiotStressDU(const mat &DSDE, const mat &U, const mat &S) {

    Tensor2<double,3,3> U_ = mat_FTensor2(U);
    Tensor2<double,3,3> delta_ = mat_FTensor2(eye(3,3));
    Tensor2<double,3,3> S_ = mat_FTensor2(S);
    Tensor4<double,3,3,3,3> DSDE_ = mat_FTensor4(DSDE);
    Tensor4<double,3,3,3,3> C_;
    
    Index<'s', 3> s;
    Index<'j', 3> j;
    Index<'p', 3> p;
    Index<'r', 3> r;
    
    Index<'i', 3> i;
    Index<'l', 3> l;
    Index<'m', 3> m;
    Index<'n', 3> n;
    
    C_(s,j,p,r) = 0.5*delta_(i,s)*(U_(l,j)*(U_(m,p)*(delta_(n,r)*DSDE_(i,l,m,n)))) 
                + 0.5*delta_(i,s)*(U_(l,j)*(delta_(m,r)*(U_(n,p)*DSDE_(i,l,m,n))))
                + 0.5*S_(s,p)*delta_(r,j) + 0.5*S_(r,j)*delta_(s,p);
    return FTensor4_mat(C_);

}

} //namespace simcoon
