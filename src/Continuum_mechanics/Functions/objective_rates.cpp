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
#include <Fastor/Fastor.h>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

// Helper: copy arma 3x3 to Fastor tensor (both column-major, 72 bytes)
static inline Fastor::Tensor<double,3,3> to_fastor2(const mat &m) {
    Fastor::Tensor<double,3,3> T;
    std::memcpy(T.data(), mat::fixed<3,3>(m).memptr(), 9 * sizeof(double));
    return T;
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

// Helper: 4th-order push-forward DSDE_ijkl = F_iL F_jH F_kM F_lN C_LHMN
static inline Fastor::Tensor<double,3,3,3,3> pushforward_4(
    const Fastor::Tensor<double,3,3,3,3> &C,
    const Fastor::Tensor<double,3,3> &F)
{
    return push_forward_4(C, F);
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

void Jaumann(mat &DR, mat &D, mat &W, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    
    mat L;
    if(DTime > simcoon::iota) {    
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Jaumann (L).");
        }          
    }   

    //decomposition of L
    D = 0.5*(L+L.t());
    W = 0.5*(L-L.t());
    
    //Jaumann
    try {
        DR = (inv(I-0.5*DTime*W))*(I+0.5*DTime*W);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Jaumann (DR).");
    }     
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
    
    mat L;
    if(DTime > simcoon::iota) {    
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Green_Naghdi (L).");
        }          
    }   
    
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    Omega = (1./DTime)*(R1-R0)*R1.t();


    try {
        DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Green_Naghdi (DR).");
    }         
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
    
    mat L;
    if(DTime > simcoon::iota) {    
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic_R (L).");
        }          
    }   
    
    //decomposition of L
    D = 0.5*(L+L.t());
    mat W = 0.5*(L-L.t());
    Omega = (1./DTime)*(R1-R0)*R1.t();

    try {
        DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside logarithmic_R (DR).");
    }    
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
            if (i!=j) {
                double f_z = 0.;
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                } else {
                    // Taylor expansion: f(r) = -s/6 + s^2/12 + O(s^3) where s = r - 1
                    f_z = s*(-1./6. + s/12.);
                }
                N_1 += f_z*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }

    N_2 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if (i!=j) {
                double g_z = 0.;
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    g_z = (1.-(pow(bi(i)/bi(j),0.5)))/(1.+(pow(bi(i)/bi(j),0.5)));
                } else {
                    // Taylor expansion: g(r) = -s/4 + 3s^2/16 + O(s^3) where s = r - 1
                    g_z = s*(-1./4. + 3.*s/16.);
                }
                N_2 += g_z*Bi_proj[i]*D*Bi_proj[j];
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

    if(DTime > simcoon::iota) {
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic_F (L).");
        }
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

    N_1 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if (i!=j) {
                double f_z = 0.;
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                } else {
                    // Taylor expansion: f(r) = -s/6 + s^2/12 + O(s^3) where s = r - 1
                    f_z = s*(-1./6. + s/12.);
                }
                N_1 += f_z*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }

    N_2 = zeros(3,3);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            if (i!=j) {
                double g_z = 0.;
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    g_z = (1.-(pow(bi(i)/bi(j),0.5)))/(1.+(pow(bi(i)/bi(j),0.5)));
                } else {
                    // Taylor expansion: g(r) = -s/4 + 3s^2/16 + O(s^3) where s = r - 1
                    g_z = s*(-1./4. + 3.*s/16.);
                }
                N_2 += g_z*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }
    
    try {
        DF = (inv(I-0.5*DTime*L))*(I+0.5*DTime*L);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside logarithmic_F (DF).");
    }         
}

void Truesdell(mat &DF, mat &D, mat &L, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    if(DTime > simcoon::iota) {    
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside Truesdell (L).");
        }          
    }      

    //Note that The "spin" is actually L (spin for rigid frames of reference, "flot" for Truesdell)    
    D = 0.5*(L+L.t());
    
    //Truesdell
    try {
        DF = (inv(I-0.5*DTime*L))*(I+0.5*DTime*L);
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
            if (i!=j) {
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                } else {
                    // Taylor expansion: f(r) = -s/6 + s^2/12 + O(s^3) where s = r - 1
                    f_z = s*(-1./6. + s/12.);
                }
                BBBB = BBBB + f_z*linearop_eigsym(Bi.col(i),Bi.col(j));
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
            if (i!=j) {
                f_z = (sqrt(bi(j)) - sqrt(bi(i)))/(sqrt(bi(j)) + sqrt(bi(i)));
                BBBB = BBBB + f_z*linearop_eigsym(Bi.col(i),Bi.col(j));
            }
        }
    }
    return BBBB;
}

void logarithmic(mat &DR, mat &D, mat &Omega, const double &DTime, const mat &F0, const mat &F1) {
    mat I = eye(3,3);
    mat L = zeros(3,3);

    if(DTime > simcoon::iota) {    
        try {
            L = (1./DTime)*(F1-F0)*inv(F1);
        } catch (const std::runtime_error &e) {
            cerr << "Error in inv: " << e.what() << endl;
            throw simcoon::exception_inv("Error in inv function inside logarithmic (L).");
        }          
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
            if (i!=j) {
                double f_z = 0.;
                double s = bi(i)/bi(j) - 1.;
                if (fabs(s) > 1.e-4) {
                    f_z = (1.+(bi(i)/bi(j)))/(1.-(bi(i)/bi(j)))+2./log(bi(i)/bi(j));
                } else {
                    // Taylor expansion: h(r) = -s/6 + s^2/12 + O(s^3) where s = r - 1
                    f_z = s*(-1./6. + s/12.);
                }
                N += f_z*Bi_proj[i]*D*Bi_proj[j];
            }
        }
    }
    Omega = W + N;

    try {
        DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside logarithmic (DR).");
    }       
}

mat Delta_log_strain(const mat &D, const mat &Omega, const double &DTime) {
    mat I = eye(3,3);
    mat DR;
    try {
        DR = (inv(I-0.5*DTime*Omega))*(I+0.5*DTime*Omega);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside logarithmic (DR).");
    }           
    return 0.5*(D+(DR*D*DR.t()))*DTime;
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DtauDe_2_DSDE(const mat &Lt, const mat &B, const mat &F, const mat &tau){
    
    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DtauDe_2_DSDE.");
    }   
    auto invF_ = to_fastor2(invF);
    auto tau_ = to_fastor2(tau);
    auto Dtau_logDD = to_fastor4(Lt);
    auto B_ = to_fastor4(B);

    // Dtau_LieDD = Dtau_logDD + B-correction
    auto Dtau_LieDD = apply_B_correction(Dtau_logDD, B_, tau_, true);
    // DSDE = pull-back of Dtau_LieDD
    auto DSDE = pullback_4(Dtau_LieDD, invF_);
    return fastor4_to_voigt(DSDE);
}

mat Dtau_LieDD_2_DSDE(const mat &Lt, const mat &F){
    
    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside Dtau_LieDD_2_DSDE.");
    }   

    auto invF_ = to_fastor2(invF);
    auto Dtau_LieDD = to_fastor4(Lt);
    auto DSDE = pullback_4(Dtau_LieDD, invF_);
    return fastor4_to_voigt(DSDE);
}

mat DtauDe_JaumannDD_2_DSDE(const mat &Lt, const mat &F, const mat &tau){
    
    mat invF;
    try {
        invF = inv(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside DtauDe_JaumannDD_2_DSDE.");
    }   
    auto invF_ = to_fastor2(invF);
    auto tau_ = to_fastor2(tau);
    auto Dtau_JaumannDD = to_fastor4(Lt);

    // Dtau_LieDD = Dtau_JaumannDD - Jaumann correction
    auto Dtau_LieDD = apply_jaumann_correction(Dtau_JaumannDD, tau_, false);
    // DSDE = pull-back of Dtau_LieDD
    auto DSDE = pullback_4(Dtau_LieDD, invF_);
    return fastor4_to_voigt(DSDE);
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DsigmaDe_2_DSDE(const mat &Lt, const mat &B, const mat &F, const mat &sigma){
    
    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DsigmaDe_2_DSDE.");
    }     
    return DtauDe_2_DSDE(J*Lt, B, F, Cauchy2Kirchoff(sigma, F, J));
}

//This function computes the tangent modulus that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E to the tangent modulus that links the Kirchoff elastic tensor and logarithmic strain, through the log rate and the and the transformation gradient F
mat DsigmaDe_2_DSDE(const mat &Lt, const mat &F, const mat &sigma){

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DsigmaDe_2_DSDE.");
    }
    mat B = get_BBBB(F);
    return DtauDe_2_DSDE(J*Lt, B, F, Cauchy2Kirchoff(sigma, F, J));
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

    mat B = get_BBBB_GN(F);
    return DtauDe_2_DSDE(Lt, B, F, tau);
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

mat DSDE_2_DtauDe(const mat &DSDE, const mat &B, const mat &F, const mat &tau) {
    
    auto F_ = to_fastor2(F);
    auto tau_ = to_fastor2(tau);
    auto DSDE_ = to_fastor4(DSDE);
    auto B_ = to_fastor4(B);

    // Push-forward DSDE to get Dtau_LieDD, then subtract B-correction
    auto Dtau_LieDD = pushforward_4(DSDE_, F_);
    auto C = apply_B_correction(Dtau_LieDD, B_, tau_, false);
    return fastor4_to_voigt(C);
}

mat DSDE_2_DsigmaDe(const mat &DSDE, const mat &B, const mat &F, const mat &sigma) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_DsigmaDe.");
    }   
    return (1./J)*DSDE_2_DtauDe(DSDE, B, F, Cauchy2Kirchoff(sigma, F, J));
}

//This function computes the tangent modulus that links the Lie derivative of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E) and the transformation gradient F
mat DSDE_2_Dtau_LieDD(const mat &DSDE, const mat &F) {

    auto F_ = to_fastor2(F);
    auto DSDE_ = to_fastor4(DSDE);
    auto C = pushforward_4(DSDE_, F_);
    return fastor4_to_voigt(C);
}

mat DSDE_2_Dsigma_LieDD(const mat &DSDE, const mat &F) {

    double J;
    try {
        J = det(F);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside DSDE_2_DsigmaDe_LieDD.");
    }
    return (1./J)*DSDE_2_Dtau_LieDD(DSDE, F);
}

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the Saint-Venant Kirchoff elastic tensor (that links the Piola-Kirchoff II stress S to the Green-Lagrange stress E), the transformation gradient F and the Kirchoff stress tau
mat DSDE_2_Dtau_JaumannDD(const mat &DSDE, const mat &F, const mat &tau) {

    auto F_ = to_fastor2(F);
    auto tau_ = to_fastor2(tau);
    auto DSDE_ = to_fastor4(DSDE);

    // Push-forward + Jaumann correction
    auto Dtau_LieDD = pushforward_4(DSDE_, F_);
    auto C = apply_jaumann_correction(Dtau_LieDD, tau_, true);
    return fastor4_to_voigt(C);
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

    mat Dtau_LieDD = DSDE_2_Dtau_LieDD(DSDE, F);
    return Dtau_LieDD_Dtau_logarithmicDD(Dtau_LieDD, F, tau);
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

//This function computes the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D, from the tangent modulus that links the Jaumann rate of the Kirchoff stress tau to the rate of deformation D and the Kirchoff stress tau
mat Dtau_LieDD_Dtau_JaumannDD(const mat &Dtau_LieDD, const mat &tau) {

    auto tau_ = to_fastor2(tau);
    auto Dtau_LieDD_ = to_fastor4(Dtau_LieDD);
    auto result = apply_jaumann_correction(Dtau_LieDD_, tau_, true);
    return fastor4_to_voigt(result);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dtau_LieDD_Dtau_objectiveDD(const mat &Dtau_LieDD, const mat &B, const mat &tau) {

    auto tau_ = to_fastor2(tau);
    auto Dtau_LieDD_ = to_fastor4(Dtau_LieDD);
    auto B_ = to_fastor4(B);

    // Dtau_objectiveDD = Dtau_LieDD - B-correction
    auto result = apply_B_correction(Dtau_LieDD_, B_, tau_, false);
    return fastor4_to_voigt(result);
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

    auto sigma_ = to_fastor2(sigma);
    auto Dsigma_LieDD_ = to_fastor4(Dsigma_LieDD);
    auto result = apply_jaumann_correction(Dsigma_LieDD_, sigma_, true);
    return fastor4_to_voigt(result);
}

//This function computes the tangent modulus that links the Lie rate of the Kirchoff stress tau to the rate of deformation D to the logarithmic rate of the Kirchoff stress and the rate of deformation D
mat Dsigma_LieDD_Dsigma_objectiveDD(const mat &Dsigma_LieDD, const mat &B, const mat &sigma) {

    auto sigma_ = to_fastor2(sigma);
    auto Dsigma_LieDD_ = to_fastor4(Dsigma_LieDD);
    auto B_ = to_fastor4(B);

    // Dsigma_objectiveDD = Dsigma_LieDD - B-correction
    auto result = apply_B_correction(Dsigma_LieDD_, B_, sigma_, false);
    return fastor4_to_voigt(result);
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

    auto U_ = to_fastor2(U);
    auto S_ = to_fastor2(S);
    auto DSDE_ = to_fastor4(DSDE);

    // C_sjpr = 0.5*δ_is*(U_lj*U_mp*δ_nr*DSDE_ilmn + U_lj*δ_mr*U_np*DSDE_ilmn)
    //        + 0.5*S_sp*δ_rj + 0.5*S_rj*δ_sp
    Fastor::Tensor<double,3,3,3,3> C;
    C.zeros();
    for (int s = 0; s < 3; ++s)
    for (int j = 0; j < 3; ++j)
    for (int p = 0; p < 3; ++p)
    for (int r = 0; r < 3; ++r) {
        double sum = 0.0;
        // i=s (delta_is), contract over l,m,n
        for (int l = 0; l < 3; ++l)
        for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n) {
            double dsde_val = DSDE_(s,l,m,n);
            // term1: U_lj * U_mp * δ_nr
            if (n == r) sum += 0.5 * U_(l,j) * U_(m,p) * dsde_val;
            // term2: U_lj * δ_mr * U_np
            if (m == r) sum += 0.5 * U_(l,j) * U_(n,p) * dsde_val;
        }
        // Jaumann-like terms
        if (r == j) sum += 0.5 * S_(s,p);
        if (s == p) sum += 0.5 * S_(r,j);
        C(s,j,p,r) = sum;
    }
    return fastor4_to_voigt(C);

}

} //namespace simcoon
