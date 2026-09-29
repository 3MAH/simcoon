///@file Zener_fast.hpp
///@brief User subroutine for Zener viscoelastic model with  N viscoelastic Kelvin branches in series in 3D case, with thermoelastic effect
///@brief This implementation uses a single scalar internal variable for the evaluation of the viscoelastic strain increment
///@author Chemisky, Chatzigeorgiou
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>
using namespace std;
using namespace arma;

// Model, props and statev layout: see the Doxygen block in
// simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Zener_Nfast.hpp

namespace simcoon {
    
void umat_zener_Nfast_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double E0 = props(2);
    double nu0 = props(3);
    double alpha_iso = props(4);
    int N_kelvin = int(props(5));
    
    vec E_visco = zeros(N_kelvin);
    vec nu_visco = zeros(N_kelvin);
    vec etaB_visco = zeros(N_kelvin);
    vec etaS_visco = zeros(N_kelvin);
    
    for (int i=0; i<N_kelvin; i++) {
        E_visco(i) = props(6+i*4);
        nu_visco(i) = props(6+i*4+1);
        etaB_visco(i) = props(6+i*4+2);
        etaS_visco(i) = props(6+i*4+3);
    }
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Define the viscoelastic stiffness
    mat L0 = L_iso(E0, nu0, "Enu");
    ///@brief Temperature initialization
    double T_init = statev(0);
    
    //From the statev to the internal variables
    vec EV = zeros(6);
    
    std::vector<vec> EV_i(N_kelvin);
    vec v = zeros(N_kelvin);
    
    for (int i=0; i<N_kelvin; i++) {
        v(i) = statev(i*7+7);
        EV_i[i] = zeros(6);
        EV_i[i](0) = statev(i*7+7+1);
        EV_i[i](1) = statev(i*7+7+2);
        EV_i[i](2) = statev(i*7+7+3);
        EV_i[i](3) = statev(i*7+7+4);
        EV_i[i](4) = statev(i*7+7+5);
        EV_i[i](5) = statev(i*7+7+6);
    }
    
    //Rotation of internal variables (tensors)
    for (int i=0; i<N_kelvin; i++) {
        EV_i[i] = rotate_strain(EV_i[i], DR);
    }
    
    vec sigma_start = sigma;
    std::vector<vec> DEV_i(N_kelvin);

    std::vector<mat> L_i(N_kelvin);
    std::vector<mat> H_i(N_kelvin);
    
    for (int i=0; i<N_kelvin; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
    }
    
    
    if(start) { //Initialization
        T_init = T;
        for (int i=0; i<N_kelvin; i++) {
            EV_i[i] = zeros(6);
        }
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
    const std::vector<vec> EV_i_start = EV_i;

    // Implicit (backward-Euler) step of the Kelvin branches in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const vec eps_e = Etot + DEtot - alpha*(T + DT - T_init);
    const LinearViscoStep st = kelvin_series_step(L0, L_i, H_i, EV_i_start, eps_e, alpha, DTime);
    EV = zeros(6);
    for (int i=0; i<N_kelvin; i++) {
        EV_i[i] = st.EV_i[i];
        DEV_i[i] = EV_i[i] - EV_i_start[i];
        v(i) += norm_strain(DEV_i[i]);
        EV += EV_i[i];
    }
    sigma = el_pred(L0, eps_e - EV, ndi);
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

    // Kelvin branch forces sigma - L_i EV_i (the viscous stresses), as in the mechanical twin
    double Dgamma_loc = 0.;
    vec dDgamma_dE = zeros(6);
    double dDgamma_dT = 0.;
    for (int i=0; i<N_kelvin; i++) {
        const vec A_mid2 = (sigma_start - L_i[i]*EV_i_start[i]) + (sigma - L_i[i]*EV_i[i]);
        const mat &dEVdE = st.dEVdE_i[i];
        const vec &dEVdT = st.dEVdT_i[i];
        Dgamma_loc += 0.5*sum(A_mid2%DEV_i[i]);
        dDgamma_dE += 0.5*((st.dSdE - L_i[i]*dEVdE).t()*DEV_i[i] + dEVdE.t()*A_mid2);
        dDgamma_dT += 0.5*(sum((st.dSdT - L_i[i]*dEVdT)%DEV_i[i]) + sum(dEVdT%A_mid2));
    }

    // heat source of the actual increments and its exact derivatives (linear_viscoelastic.hpp)
    viscous_heat_source(Dgamma_loc, dDgamma_dE, dDgamma_dT, st, alpha, sigma, sigma_start, T, DT,
                        rho*c_p, DTime, r, drdE, drdT);
    
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
    //From the statev to the internal variables
    statev(1) = EV(0);
    statev(2) = EV(1);
    statev(3) = EV(2);
    statev(4) = EV(3);
    statev(5) = EV(4);
    statev(6) = EV(5);
    
    for (int i=0; i<N_kelvin; i++) {
        statev(i*7+7) = v(i);
        statev(i*7+7+1) = EV_i[i](0);
        statev(i*7+7+2) = EV_i[i](1);
        statev(i*7+7+3) = EV_i[i](2);
        statev(i*7+7+4) = EV_i[i](3);
        statev(i*7+7+5) = EV_i[i](4);
        statev(i*7+7+6) = EV_i[i](5);
    }
}

} //namespace smart

