///@file Zener_fast.cpp
///@brief User subroutine for Zener viscoelastic model in 3D case, with thermoelastic effect(Poynting-Thomson model)
///@brief This implementation uses a single scalar internal variable for the evaluation of the viscoelastic strain increment
///@author Chemisky, Chatzigeorgiou
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
using namespace std;
using namespace arma;

///@brief The viscoelastic Zener model requires 10 constants:
//      -------------------
///@brief      props(0) = rho           - density
///@brief      props(1) = c_p           - specific heat capacity
///@brief      props(2) = E0            - Thermoelastic Young's modulus
///@brief      props(3) = nu0           - Thermoelastic Poisson's ratio
///@brief      props(4) = alpha_iso     - Thermoelastic CTE
///@brief      props(5) = E1            - Viscoelastic Young modulus of Zener branch
///@brief      props(6) = nu0           - Viscoelastic Poisson ratio of Zener branch
///@brief      props(7) = etaB1         - Viscoelastic Bulk viscosity of Zener branch
///@brief      props(8) = etaS1         - Viscoelastic shear viscosity of Zener branch

///@brief Number of statev required for thermoelastic constitutive law

namespace simcoon {
    
void umat_zener_fast_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
    {
        
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    vec sigma_start = sigma;
    
    ///@brief Temperature initialization
    double T_init = statev(0);
    double v = statev(1);
    
    //From the statev to the internal variables
    vec EV1 = zeros(6);
    EV1(0) = statev(2);
    EV1(1) = statev(3);
    EV1(2) = statev(4);
    EV1(3) = statev(5);
    EV1(4) = statev(6);
    EV1(5) = statev(7);
    
    //Rotation of internal variables (tensors)
    EV1 = rotate_strain(EV1, DR);
    
    //From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double E0 = props(2);
    double nu0 = props(3);
    double alpha_iso = props(4);
    double E1 = props(5);
    double nu1 = props(6);
    double etaB1 = props(7);
    double etaS1 = props(8);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Define the viscoelastic stiffness
    mat L0 = L_iso(E0, nu0, "Enu");
    mat L1 = L_iso(E1, nu1, "Enu");
    
    mat H1 = H_iso(etaB1, etaS1);                  //dimension of stiffness tensor
    mat invH1 = inv(H1);
    
    if(start) { //Initialization
        T_init = T;
        EV1 = zeros(6);
        sigma = zeros(6);
        sigma_start = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
        
        Wt = 0.;
        Wt_r = 0.;
        Wt_ir = 0.;
    }
    
    //Additional parameters and variables
    double c_0 = rho*c_p;
        
    //Variables at the start of the increment
    vec DEV1 = zeros(6);
    vec EV1_start = EV1;
    vec A_v_start = L1*EV1;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = v;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);        
    
    //Determination of the initial, predicted stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EV1;
    vec DEel = DEtot - alpha*(DT);
    if (ndi == 1) {
        sigma(0) = sigma_start(0) + E0*DEel(0);
    }
    else if (ndi == 2) {
        sigma(0) = sigma_start(0) + E0/(1. - (nu0*nu0))*(DEel(0)) + nu0*(DEel(1));
        sigma(1) = sigma_start(1) + E0/(1. - (nu0*nu0))*(DEel(1)) + nu0*(DEel(0));
        sigma(3) = sigma_start(3) + E0/(1.+nu0)*0.5*DEel(3);
    }
    else
    sigma = sigma_start + (L0*DEel);
    

    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
    double dPhidv=0.;
    vec dPhidEv = zeros(6);
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec sigma_tildeV1 = invH1*(sigma-L1*EV1);
    vec Lambdav = eta_norm_strain(sigma_tildeV1);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = L0*Lambdav;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        v = s_j(0);

        sigma_tildeV1 = (invH1*(sigma-L1*EV1))%Ir05();
        Lambdav = eta_norm_stress(sigma_tildeV1);
        dPhidsigma = invH1*(eta_norm_stress(sigma_tildeV1)%Ir05()); //Dimension of strain (similar to Lambda in general)
        
        if (DTime > simcoon::iota) {
            Phi(0) = norm_stress(sigma_tildeV1) - Ds_j(0)/DTime;
            dPhidv = -1.*sum((dPhidsigma%Ir2())%(L1*Lambdav))-1./DTime;
        }
        else {
            //No time, no flow: the branch is INACTIVE. The stationary condition
            //Phi = ||flow|| has root EV = eps, committed as relaxed by the zero-time probe.
            Phi(0) = 0.;
            dPhidv = -1.;
        }
        kappa_j[0] = L0*Lambdav;
        
        K(0,0) = dPhidv;
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = norm_stress(sigma_tildeV1);
        if (Y_crit(0) < simcoon::precision_umat) {
            Y_crit(0) = simcoon::precision_umat;
        }
        
        Newton_Raphon(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        EV1 = EV1 + ds_j(0)*Lambdav;
        DEV1 = DEV1 + ds_j(0)*Lambdav;
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV1;
        DEel = DEtot - alpha*(DT) - DEV1;
        if (ndi == 1) {
            sigma(0) = sigma_start(0) + E0*DEel(0);
        }
        else if (ndi == 2) {
            sigma(0) = sigma_start(0) + E0/(1. - (nu0*nu0))*(DEel(0)) + nu0*(DEel(1));
            sigma(1) = sigma_start(1) + E0/(1. - (nu0*nu0))*(DEel(1)) + nu0*(DEel(0));
            sigma(3) = sigma_start(3) + E0/(1.+nu0)*0.5*DEel(3);
        }
        else
        sigma = sigma_start + (L0*DEel);
    }
    
    //Computation of the tangent modulus — continuum operator via shared helper (doc §7.4).
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const std::vector<vec> dPhidsigma_l = { dPhidsigma };
    const ContinuumTangent ct = assemble_continuum_tangent(Bhat, kappa_j, dPhidsigma_l, Ds_j, L0);
    dSdE = ct.Lt;
    const std::vector<vec>& P_epsilon = ct.P_epsilon;

    std::vector<double> P_theta(1);
    P_theta[0] = dPhidtheta - sum(dPhidsigma%(L0*alpha));

    dSdT = -1.*L0*alpha - (kappa_j[0]*P_theta[0]);
        
    //computation of the internal energy production
    double eta_r = c_0*log((T+DT)/T_init) + sum(alpha%sigma);
    double eta_r_start = c_0*log(T/T_init) + sum(alpha%sigma_start);

    double eta_ir = 0.;
    double eta_ir_start = 0.;
    
    double eta = eta_r + eta_ir;
    double eta_start = eta_r_start + eta_ir_start;
    
    double Deta = eta - eta_start;
    double Deta_r = eta_r - eta_r_start;
    double Deta_ir = eta_ir - eta_ir_start;
        
    // r = (Dgamma - Tm alpha:(sigma - sigma_start) - rho c_p DT)/DTime, midpoint Tm as Wt: the
    // heat source of the actual increments (a linearisation in (DEtot, DT) would vanish during a
    // strain hold while the branches relax). drdE/drdT differentiate it with the flow directions
    // frozen and the continuum dDs_i/dE = P_epsilon[i], dDs_i/dT = P_theta[i]: approximate, like
    // dSdE (no algorithmic tangent here); the converged r is exact.
    // Kelvin branch force sigma - L1 EV1 (the viscous stress), as in the mechanical twin
    const vec A_mid2 = (sigma_start - L1*EV1_start) + (sigma - L1*EV1);
    const vec LLambda = L1*Lambdav;
    const double Dgamma_loc = 0.5*sum(A_mid2%DEV1);
    const vec dDgamma_dE = 0.5*(dSdE.t()*DEV1 - sum(DEV1%LLambda)*P_epsilon[0] + sum(A_mid2%Lambdav)*P_epsilon[0]);
    const double dDgamma_dT = 0.5*(sum(dSdT%DEV1) - sum(DEV1%LLambda)*P_theta[0] + sum(A_mid2%Lambdav)*P_theta[0]);
    if (DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        const double Tm = T + 0.5*DT;
        r = (Dgamma_loc - Tm*sum(alpha%(sigma - sigma_start)) - rho*c_p*DT)/DTime;
        drdE = (dDgamma_dE - Tm*(dSdE.t()*alpha))/DTime;
        drdT = (dDgamma_dT - 0.5*sum(alpha%(sigma - sigma_start)) - Tm*sum(alpha%dSdT) - rho*c_p)/DTime;
    }
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot) - Dgamma_loc;
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir += (T+0.5*DT)*Deta_ir;
        
    //Return the statev;
    statev(0) = T_init;
    statev(1) = v;
    statev(2) = EV1(0);
    statev(3) = EV1(1);
    statev(4) = EV1(2);
    statev(5) = EV1(3);
    statev(6) = EV1(4);
    statev(7) = EV1(5);
}
    
} //namespace smart

