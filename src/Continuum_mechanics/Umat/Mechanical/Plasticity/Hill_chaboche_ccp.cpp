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

///@file plastic_kin_iso_ccp.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Linear Kinematical hardening coupled with a power-law hardenig is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Hill_chaboche_ccp.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

using namespace std;
using namespace arma;

namespace simcoon {
    
///@brief The mechanical elastic-plastic UMAT with kinematic + isotropic hardening requires 7 constants:

///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE
///@brief props[3] : J2 equivalent yield stress limit : sigmaY
///@brief props[4] : hardening parameter k
///@brief props[5] : exponent m
///@brief props[6] : linear kinematical hardening h

///@brief The elastic-plastic UMAT with kinematic + isotropic hardening requires 14 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)
///@brief statev[8] : Backstress 11: X(0,0)
///@brief statev[9] : Backstress 11: X(1,1)
///@brief statev[10] : Backstress 11: X(2,2)
///@brief statev[11] : Backstress 11: X(0,1)
///@brief statev[12] : Backstress 11: X(0,2)
///@brief statev[13] : Backstress 11: X(1,2)


void umat_hill_chaboche_CCP(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double G = props(2);    
    double alpha_iso = props(3);
    double sigmaY = props(4);
    double Q=props(5);
    double b=props(6);
    double C_1 = props(7);
    double D_1 = props(8);
    double C_2 = props(9);
    double D_2 = props(10);

    double F_hill = props(11);
    double G_hill = props(12);
    double H_hill = props(13);
    double L_hill = props(14);
    double M_hill = props(15);
    double N_hill = props(16);
    
    vec Hill_params = {F_hill,G_hill,H_hill,L_hill,M_hill,N_hill};    
    //double X_0 = props(10);
    //double ep_eq0 = props(11);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
       
    ///@brief Temperature initialization
    double T_init = statev(0);
    //From the statev to the internal variables
    double p = statev(1);
    vec EP(6);
    EP(0) = statev(2);
    EP(1) = statev(3);
    EP(2) = statev(4);
    EP(3) = statev(5);
    EP(4) = statev(6);
    EP(5) = statev(7);
    
    ///@brief a is the internal variable associated with kinematical hardening
    vec a_1 = zeros(6);
    a_1(0) = statev(8);
    a_1(1) = statev(9);
    a_1(2) = statev(10);
    a_1(3) = statev(11);
    a_1(4) = statev(12);
    a_1(5) = statev(13);

    vec a_2 = zeros(6);
    a_2(0) = statev(14);
    a_2(1) = statev(15);
    a_2(2) = statev(16);
    a_2(3) = statev(17);
    a_2(4) = statev(18);
    a_2(5) = statev(19);
    
    vec X_1 = zeros(6);
    X_1(0) = statev(20);
    X_1(1) = statev(21);
    X_1(2) = statev(22);
    X_1(3) = statev(23);
    X_1(4) = statev(24);
    X_1(5) = statev(25);

    vec X_2 = zeros(6);
    X_2(0) = statev(26);
    X_2(1) = statev(27);
    X_2(2) = statev(28);
    X_2(3) = statev(29);
    X_2(4) = statev(30);
    X_2(5) = statev(31);
    
    double Hp = statev(32);
    
    vec X = X_1 + X_2;
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    a_1 = rotate_strain(a_1, DR);
    a_2 = rotate_strain(a_2, DR);
    
    //Elstic stiffness tensor
    L = L_cubic(E, nu, G, "EnuG");
        
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        stress = vide;
        EP = vide;
        a_1 = vide;
        a_2 = vide;
        p = 0.;
        Hp = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
    //Additional parameters and variables
    double dHpdp=0.;
    
    if (p > simcoon::iota)	{
        dHpdp = b*(Q-Hp);
    }
    else {
        dHpdp = 0.;
    }
    
    //Variables values at the start of the increment
    vec stress_start = stress;
    vec EP_start = EP;
    vec a_1start = a_1;
    vec a_2start = a_2;
    vec X_1start = X_1;
    vec X_2start = X_2;
    
    double A_p_start = -Hp;
    vec A_a1_start = -X_1start;
    vec A_a2_start = -X_2start;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    stress = el_pred(L, Eel, ndi);
    
    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
    double dPhidp=0.;
    vec dPhida_1 = zeros(6);
    vec dPhida_2 = zeros(6);
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec Lambdap = dHill_stress(stress-X,Hill_params);
    vec Lambdaa_1 = dHill_stress(stress-X,Hill_params) - D_1*a_1;
    vec Lambdaa_2 = dHill_stress(stress-X,Hill_params) - D_2*a_2;

    std::vector<vec> kappa_j(1);
    kappa_j[0] = L*Lambdap;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    // tangent_mode 2: closest-point projection replaces the CCP loop (doc §cpp_return_mapping).
    // Backward-Euler closed forms resolve the state at every iterate:
    //   a_i = (a_{n,i} + Dl*eta(xi)) / (1 + Dl*D_i)   (Armstrong-Frederick recall)
    //   Hp  = (Hp_n + b*Q*Dl) / (1 + b*Dl)            (Voce isotropic hardening)
    // iterated on xi = sigma - X; dLambda/dsigma is the TOTAL derivative of the
    // inner-consistent map (central FD), carrying the dLambda/dX chain.
    ReturnMappingResult rm;
    const bool use_cpp = (tangent_mode == 2) && (ndi == 3);
    if (use_cpp) {
        const double p_n = s_j(0);
        const double Hp_n = Hp;
        const vec a_1n = a_1, a_2n = a_2;
        const vec X_1n = X_1, X_2n = X_2;
        // Backward-Euler state at (sig, Dl): fixed point on X; returns nothing, updates captures.
        auto solve_state = [&](const vec &sig, double Dl) {
            for (int it_fp = 0; it_fp < 20; it_fp++) {
                const vec eta_xi = dHill_stress(sig - X, Hill_params);
                const vec a1l = (a_1n + Dl*eta_xi)/(1. + Dl*D_1);
                const vec a2l = (a_2n + Dl*eta_xi)/(1. + Dl*D_2);
                const vec X1l = X_1n + (2./3.)*C_1*((a1l - a_1n) % Ir05());
                const vec X2l = X_2n + (2./3.)*C_2*((a2l - a_2n) % Ir05());
                const vec Xn = X1l + X2l;
                const double dX = norm(Xn - X, 2);
                a_1 = a1l; a_2 = a2l; X_1 = X1l; X_2 = X2l; X = Xn;
                if (dX < 1e-12*(norm(Xn, 2) + 1.)) break;
            }
            Hp = (Hp_n + b*Q*Dl)/(1. + b*Dl);
            dHpdp = b*(Q - Hp);
        };
        // Guarded probe of Phi at (sig, Dl) with the inner state solved, captures restored.
        auto probe_Phi = [&](const vec &sig, double Dl) {
            const vec Xs1 = X_1, Xs2 = X_2, Xs = X, as1 = a_1, as2 = a_2;
            const double Hs = Hp, dHs = dHpdp;
            solve_state(sig, Dl);
            const double v = Hill_stress(sig - X, Hill_params) - Hp - sigmaY;
            X_1 = Xs1; X_2 = Xs2; X = Xs; a_1 = as1; a_2 = as2;
            Hp = Hs; dHpdp = dHs;
            return v;
        };
        // With the inner-consistent state map (X, Hp) = V_hat(sigma, Dl), ALL derivative
        // callbacks must be TOTAL derivatives of the composed map (partial-only inputs
        // leave an O(dX/dsigma) error in the consistent tangent).
        ReturnMechanism mech;
        mech.Phi            = [&](const vec &sig) { return Hill_stress(sig - X, Hill_params) - Hp - sigmaY; };
        mech.dPhi_dsigma    = [&](const vec &sig) {
            const double hfd = 1.e-5*(norm(sig, 2) + 1.);
            const double Dl = p - p_n;
            vec g(6);
            for (int c6 = 0; c6 < 6; c6++) {
                vec sp = sig, sm = sig;
                sp(c6) += hfd;
                sm(c6) -= hfd;
                g(c6) = (probe_Phi(sp, Dl) - probe_Phi(sm, Dl))/(2.*hfd);
            }
            return g;
        };
        mech.Lambda         = [&](const vec &sig) { return dHill_stress(sig - X, Hill_params); };
        mech.dLambda_dsigma = [&](const vec &sig) {
            const double hfd = 1.e-5*(norm(sig, 2) + 1.);
            const double Dl = p - p_n;
            const vec Xsave1 = X_1, Xsave2 = X_2, Xsave = X, asave1 = a_1, asave2 = a_2;
            const double Hpsave = Hp, dHsave = dHpdp;
            mat D(6, 6);
            for (int c6 = 0; c6 < 6; c6++) {
                vec sp = sig, sm = sig;
                sp(c6) += hfd;
                sm(c6) -= hfd;
                solve_state(sp, Dl);
                const vec ep = dHill_stress(sp - X, Hill_params);
                solve_state(sm, Dl);
                const vec em = dHill_stress(sm - X, Hill_params);
                D.col(c6) = (ep - em)/(2.*hfd);
            }
            X_1 = Xsave1; X_2 = Xsave2; X = Xsave; a_1 = asave1; a_2 = asave2;
            Hp = Hpsave; dHpdp = dHsave;
            return D;
        };
        rm = closest_point_return_mapping(stress, L, mech,
                 [&](const vec &sig, double Dl) {
                     p = p_n + Dl;
                     solve_state(sig, Dl);
                     return true;
                 },
                 // TOTAL dPhi/dDlambda at fixed sigma (Voce + both backstress chains).
                 [&](const vec &sig, double Dl) {
                     const double hDl = 1.e-6*(fabs(Dl) + 1.e-8);
                     const double Dlp = Dl + hDl;
                     const double Dlm = std::max(Dl - hDl, 0.);
                     return (probe_Phi(sig, Dlp) - probe_Phi(sig, Dlm))/(Dlp - Dlm);
                 },
                 sigmaY,
                 // multiplier-side state chain c = Dl*L*dLambda/dDl (FD of the consistent map):
                 // required for the exact consistent tangent with backstress (doc eq:Lambda_tilde_state).
                 [&](const vec &sig, double Dl) -> vec {
                     const vec Xs1 = X_1, Xs2 = X_2, Xs = X, as1 = a_1, as2 = a_2;
                     const double Hs = Hp, dHs = dHpdp;
                     const double hDl = 1.e-6*(fabs(Dl) + 1.e-8);
                     const double Dlp = Dl + hDl;
                     const double Dlm = std::max(Dl - hDl, 0.);
                     solve_state(sig, Dlp);
                     const vec ep = dHill_stress(sig - X, Hill_params);
                     solve_state(sig, Dlm);
                     const vec em = dHill_stress(sig - X, Hill_params);
                     X_1 = Xs1; X_2 = Xs2; X = Xs; a_1 = as1; a_2 = as2;
                     Hp = Hs; dHpdp = dHs;
                     return vec(Dl*(L*((ep - em)/(Dlp - Dlm))));
                 });
        if (rm.converged) {
            stress = rm.sigma;
            Ds_j(0) = rm.Dlambda(0);
            s_j(0) = p_n + Ds_j(0);
            p = s_j(0);
            dPhidsigma = rm.dPhidsigma_l[0];
            Lambdap = dHill_stress(stress - X, Hill_params);       // flow direction (NOT the total dPhi/dsigma)
            Lambdaa_1 = Lambdap - D_1*a_1;
            Lambdaa_2 = Lambdap - D_2*a_2;
            kappa_j = rm.kappa_j;
            dPhida_1 = -1.*(2./3.)*C_1*(Lambdap % Ir05());
            dPhida_2 = -1.*(2./3.)*C_2*(Lambdap % Ir05());
            K(0,0) = -dHpdp + sum(dPhida_1 % Lambdaa_1) + sum(dPhida_2 % Lambdaa_2);
            EP = EP_start + Ds_j(0)*Lambdap;   // CPP: flow at the CONVERGED stress
            // a_i, X_i, Hp already inner-consistent from the last state solve.
        }
        else {
            tnew_dt = 0.5;                     // step-cut; never a silent CCP fallback
            stress = stress_start;
            EP = EP_start;
            a_1 = a_1n; a_2 = a_2n; X_1 = X_1n; X_2 = X_2n; X = X_1 + X_2;
            Hp = Hp_n;
            dHpdp = (p_n > simcoon::iota) ? b*(Q - Hp) : 0.;
            s_j(0) = p_n;
            p = p_n;
            Ds_j(0) = 0.;
            dPhidsigma = dHill_stress(stress - X, Hill_params);
            Lambdap = dPhidsigma;
            Lambdaa_1 = dPhidsigma - D_1*a_1;
            Lambdaa_2 = dPhidsigma - D_2*a_2;
            kappa_j[0] = L*Lambdap;
            dPhida_1 = -1.*(2./3.)*C_1*(dPhidsigma % Ir05());
            dPhida_2 = -1.*(2./3.)*C_2*(dPhidsigma % Ir05());
            K(0,0) = -dHpdp + sum(dPhida_1 % Lambdaa_1) + sum(dPhida_2 % Lambdaa_2);
        }
    }
    else {
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        p = s_j(0);
        if (p > simcoon::iota)	{
            dHpdp = b*(Q-Hp);
        }
        else {
            dHpdp = 0.;
        }
        dPhidsigma = dHill_stress(stress-X,Hill_params);
        dPhidp = -1.*dHpdp;
//        dPhidp = -1.*dHpdp;
        dPhida_1 = -1.*(2./3.)*C_1*(dHill_stress(stress-X,Hill_params)%Ir05());
        dPhida_2 = -1.*(2./3.)*C_2*(dHill_stress(stress-X,Hill_params)%Ir05());
//        dPhida = 0.*(eta_stress(stress - X)%Ir05());
        
        //compute Phi and the derivatives
        Phi(0) = Hill_stress(stress-X,Hill_params) - Hp - sigmaY;
        
        Lambdap = dHill_stress(stress-X,Hill_params);
        Lambdaa_1 = dHill_stress(stress-X,Hill_params) - D_1*a_1;
        Lambdaa_2 = dHill_stress(stress-X,Hill_params) - D_2*a_2;
        kappa_j[0] = L*Lambdap;
        
        K(0,0) = dPhidp + sum(dPhida_1%Lambdaa_1) + sum(dPhida_2%Lambdaa_2);
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        Hp += b*(Q-Hp)*ds_j(0);
        EP = EP + ds_j(0)*Lambdap;
        a_1 = a_1 + ds_j(0)*Lambdaa_1;
        a_2 = a_2 + ds_j(0)*Lambdaa_2;
        X_1 += ds_j(0)*(2./3.)*C_1*(Lambdaa_1%Ir05());
        X_2 += ds_j(0)*(2./3.)*C_2*(Lambdaa_2%Ir05());
        X = X_1 + X_2;
        
        //the stress is now computed using the relationship stress = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        stress = el_pred(L, Eel, ndi);
    }
    }

    //Computation of the increments of variables
    vec Dsigma = stress - stress_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];
    vec Da_1 = a_1 - a_1start;
    vec Da_2 = a_2 - a_2start;
        
    //Computation of the tangent modulus — continuum elastic-plastic operator
    //assembled via the shared leading-mechanism helper (doc §7.4).
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const std::vector<vec> dPhidsigma_l = { dPhidsigma };
    ContinuumTangent ct;
    if (use_cpp) {
        // Exact consistent tangent of the converged CPP map (doc §cpp_return_mapping).
        ct = cpp_consistent_tangent(rm, L);
    } else if (tangent_mode >= 1) {
        // Simo-Hughes algorithmic tangent (closest-point). Hill flow on effective stress (sigma-X):
        // dLambda_eps/dsigma = ddHill_stress(stress-X). Backstress state-coupling deferred (CPP, future).
        const std::vector<mat> dLambda_dsigma_l = { ddHill_stress(stress - X, Hill_params) };
        ct = assemble_algorithmic_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L, dLambda_dsigma_l);
    } else {
        ct = assemble_continuum_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L);
    }
    Lt = ct.Lt;
    const std::vector<vec>& P_epsilon = ct.P_epsilon;

    std::vector<double> P_theta(1);
    P_theta[0] = dPhidtheta - sum(dPhidsigma%(L*alpha));
    
    double A_p = -Hp;
    vec A_a1 = -X_1;
    vec A_a2 = -X_2;
    
    double Dgamma_loc = 0.5*sum((stress_start+stress)%DEP) + 0.5*(A_p_start + A_p)*Dp + 0.5*sum((A_a1_start + A_a1)%Da_1)+0.5*sum((A_a2_start + A_a2)%Da_2);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%(DEtot-DEP)) - 0.5*sum((A_a1_start + A_a1)%Da_1) - 0.5*sum((A_a2_start + A_a2)%Da_2);
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;
            
    ///@brief statev evolving variables
    //statev
    statev(0) = T_init;
    statev(1) = p;
    
    statev(2) = EP(0);
    statev(3) = EP(1);
    statev(4) = EP(2);
    statev(5) = EP(3);
    statev(6) = EP(4);
    statev(7) = EP(5);
    
    statev(8) = a_1(0);
    statev(9) = a_1(1);
    statev(10) = a_1(2);
    statev(11) = a_1(3);
    statev(12) = a_1(4);
    statev(13) = a_1(5);
    
    statev(14) = a_2(0);
    statev(15) = a_2(1);
    statev(16) = a_2(2);
    statev(17) = a_2(3);
    statev(18) = a_2(4);
    statev(19) = a_2(5);
    
    statev(20) = X_1(0);
    statev(21) = X_1(1);
    statev(22) = X_1(2);
    statev(23) = X_1(3);
    statev(24) = X_1(4);
    statev(25) = X_1(5);
    
    statev(26) = X_2(0);
    statev(27) = X_2(1);
    statev(28) = X_2(2);
    statev(29) = X_2(3);
    statev(30) = X_2(4);
    statev(31) = X_2(5);

    statev(32) = Hp;
}

    
} //namespace simcoon
