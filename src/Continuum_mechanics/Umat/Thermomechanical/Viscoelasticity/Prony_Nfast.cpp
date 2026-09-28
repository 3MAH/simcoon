///@file Prony.hpp
///@brief User subroutine for Prony series viscoelastic model in 3D case
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

// Model, props and statev layout: see the Doxygen block in
// simcoon/Continuum_mechanics/Umat/Thermomechanical/Viscoelasticity/Prony_Nfast.hpp

namespace simcoon {
    
void umat_prony_Nfast_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
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
    int N_prony = int(props(5));
    
    vec E_visco = zeros(N_prony);
    vec nu_visco = zeros(N_prony);
    vec etaB_visco = zeros(N_prony);
    vec etaS_visco = zeros(N_prony);
    
    for (int i=0; i<N_prony; i++) {
        E_visco(i) = props(6+i*4);
        nu_visco(i) = props(6+i*4+1);
        etaB_visco(i) = props(6+i*4+2);
        etaS_visco(i) = props(6+i*4+3);
    }
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Define the viscoelastic stiffness
    mat L0 = L_iso(E0, nu0, "Enu");
    mat M0 = M_iso(E0, nu0, "Enu");
    ///@brief Temperature initialization
    double T_init = statev(0);
    
    //From the statev to the internal variables
    vec EV_tilde = zeros(6);
    EV_tilde(0) = statev(1);
    EV_tilde(1) = statev(2);
    EV_tilde(2) = statev(3);
    EV_tilde(3) = statev(4);
    EV_tilde(4) = statev(5);
    EV_tilde(5) = statev(6);
    
    std::vector<vec> EV_i(N_prony);
    vec v = zeros(N_prony);
    
    for (int i=0; i<N_prony; i++) {
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
    EV_tilde = rotate_strain(EV_tilde, DR);
    for (int i=0; i<N_prony; i++) {
        EV_i[i] = rotate_strain(EV_i[i], DR);
    }
    
    vec sigma_start = sigma;
    std::vector<vec> DEV_i(N_prony);
    std::vector<vec> A_v(N_prony);
    std::vector<vec> A_v_start(N_prony);

    std::vector<mat> L_i(N_prony);
    std::vector<mat> H_i(N_prony);
    
    for (int i=0; i<N_prony; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));

        //Unconditionally, as the mechanical Prony_Nfast does: a default-constructed
        //arma::vec has size 0, so the `A_v_start[i] +=` below threw
        //"addition: incompatible matrix dimensions: 0x1 and 6x1" on the first call.
        DEV_i[i] = zeros(6);
        A_v[i] = zeros(6);
        A_v_start[i] = zeros(6);
    }
    
    
    if(start) { //Initialization
        T_init = T;
        EV_tilde = zeros(6);
        for (int i=0; i<N_prony; i++) {
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
    for (int i=0; i<N_prony; i++) {
        A_v_start[i] += L_i[i]*(Etot - alpha*(T-T_init) - EV_i[i]);   // start state: T, EV_i before the update
    }

    // Implicit (backward-Euler) step of the Maxwell branches in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const LinearViscoStep st = maxwell_parallel_step(L0, L_i, H_i, EV_i_start,
                                                     Etot + DEtot - alpha*(T + DT - T_init), alpha, DTime);
    EV_tilde = zeros(6);
    for (int i=0; i<N_prony; i++) {
        EV_i[i] = st.EV_i[i];
        DEV_i[i] = EV_i[i] - EV_i_start[i];
        v(i) += norm_strain(DEV_i[i]);
        EV_tilde += (M0*L_i[i])*EV_i[i];
    }
    const vec Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV_tilde;
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

    // branch forces L_i (E - alpha dT - EV_i): dA_i/dE = L_i (I - dEV_i/dE),
    // dA_i/dT = -L_i (alpha + dEV_i/dT)
    double Dgamma_loc = 0.;
    vec dDgamma_dE = zeros(6);
    double dDgamma_dT = 0.;
    for (int i=0; i<N_prony; i++) {
        A_v[i] = L_i[i]*(Etot + DEtot - alpha*(T+DT-T_init) - EV_i[i]);
        const vec A_mid2 = A_v_start[i] + A_v[i];
        const mat &dEVdE = st.dEVdE_i[i];
        const vec &dEVdT = st.dEVdT_i[i];
        Dgamma_loc += 0.5*sum(A_mid2%DEV_i[i]);
        dDgamma_dE += 0.5*((L_i[i]*(eye(6,6) - dEVdE)).t()*DEV_i[i] + dEVdE.t()*A_mid2);
        dDgamma_dT += 0.5*(-sum((L_i[i]*(alpha + dEVdT))%DEV_i[i]) + sum(dEVdT%A_mid2));
    }

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
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    for (int i=0; i<N_prony; i++) {
        Wm_r += -0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
    }
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir += (T+0.5*DT)*Deta_ir;
    
    //Return the statev;
    statev(0) = T_init;
    //From the statev to the internal variables
    statev(1) = EV_tilde(0);
    statev(2) = EV_tilde(1);
    statev(3) = EV_tilde(2);
    statev(4) = EV_tilde(3);
    statev(5) = EV_tilde(4);
    statev(6) = EV_tilde(5);
    
    for (int i=0; i<N_prony; i++) {
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

