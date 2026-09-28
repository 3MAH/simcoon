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
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>
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
    const vec EV1_start = EV1;

    // Implicit (backward-Euler) step of the Kelvin branch in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const LinearViscoStep st = kelvin_series_step(L0, {L1}, {H1}, {EV1_start},
                                                  Etot + DEtot - alpha*(T + DT - T_init), alpha, DTime);
    EV1 = st.EV_i[0];
    const vec DEV1 = EV1 - EV1_start;
    v += norm_strain(DEV1);
    const vec Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV1;
    sigma = el_pred(L0, Eel, ndi);
    if (tangent_mode == tangent_none) {
        dSdE = L0;
        dSdT = -L0*alpha;
    }
    else {
        dSdE = st.dSdE;
        dSdT = st.dSdT;
    }

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

    // Kelvin branch force sigma - L1 EV1 (the viscous stress), as in the mechanical twin
    const vec A_mid2 = (sigma_start - L1*EV1_start) + (sigma - L1*EV1);
    const mat &dEVdE = st.dEVdE_i[0];
    const vec &dEVdT = st.dEVdT_i[0];
    const double Dgamma_loc = 0.5*sum(A_mid2%DEV1);
    const vec dDgamma_dE = 0.5*((st.dSdE - L1*dEVdE).t()*DEV1 + dEVdE.t()*A_mid2);
    const double dDgamma_dT = 0.5*(sum((st.dSdT - L1*dEVdT)%DEV1) + sum(dEVdT%A_mid2));

    // r = (Dgamma - Tm alpha:(sigma - sigma_start) - rho c_p DT)/DTime, midpoint Tm as Wt: the heat
    // source of the actual increments (it keeps flowing while the branches relax under a strain
    // hold). With the closed-form step, drdE/drdT are its exact derivatives.
    if (DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        const double Tm = T + 0.5*DT;
        r = (Dgamma_loc - Tm*sum(alpha%(sigma - sigma_start)) - rho*c_p*DT)/DTime;
        drdE = (dDgamma_dE - Tm*(st.dSdE.t()*alpha))/DTime;
        drdT = (dDgamma_dT - 0.5*sum(alpha%(sigma - sigma_start)) - Tm*sum(alpha%st.dSdT) - rho*c_p)/DTime;
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

