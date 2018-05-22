///@file Zener_fast.hpp
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

using namespace std;
using namespace arma;

///@brief The viscoelastic Zener model requires 10 constants:
//      -------------------
///@brief      props(0) = E0            - Thermoelastic Young's modulus
///@brief      props(1) = nu0           - Thermoelastic Poisson's ratio
///@brief      props(2) = alpha_iso     - Thermoelastic CTE
///@brief      props(3) = E1            - Viscoelastic Young modulus of Zener branch i
///@brief      props(4) = nu1           - Viscoelastic Poisson ratio of Zener branch i
///@brief      props(5) = etaB1         - Viscoelastic Bulk viscosity of Zener branch i
///@brief      props(6) = etaS1         - Viscoelastic Bulk viscosity of Zener branch i

///@brief Number of statev required for thermoelastic constitutive law

namespace simcoon {
    
void umat_zener_fast(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
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
    double E0 = props(0);
    double nu0 = props(1);
    double alpha_iso = props(2);
    double E1 = props(3);
    double nu1 = props(4);
    double etaB1 = props(5);
    double etaS1 = props(6);

    
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
    }
    
    //Variables at the start of the increment
    vec DEV1 = zeros(6);
    vec EV1_start = EV1;
    vec A_v_start = sigma_start - L1*EV1;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = v;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);        
    
    //Determination of the initial, predicted stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EV1;
    sigma = el_pred(L0, Eel, ndi);

    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
    double dPhidv=0.;
    vec dPhidEv = zeros(6);
    vec dPhidsigma = zeros(6);
    
    //Compute the explicit flow direction
    vec flow_V1 = invH1*(sigma-L1*EV1);
    vec Lambdav = eta_norm_strain(flow_V1);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = L0*Lambdav;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        v = s_j(0);

        flow_V1 = (invH1*(sigma-L1*EV1));
        Lambdav = eta_norm_strain(flow_V1);
        dPhidsigma = invH1*(eta_norm_strain(flow_V1)%Ir05()); //Dimension of strain (similar to Lambda in general)
        
        if (DTime > iota) {
            Phi(0) = norm_strain(flow_V1) - Ds_j(0)/DTime;
            dPhidv = -1.*sum((dPhidsigma)%(L1*Lambdav))-1./DTime;
        }
        else {
            Phi(0) = norm_strain(flow_V1);
            dPhidv = -1.*sum((dPhidsigma)%(L1*Lambdav));
        }
        kappa_j[0] = L0*Lambdav;
        
        K(0,0) = dPhidv;
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = norm_stress(flow_V1);
        if (Y_crit(0) < precision_umat) {
            Y_crit(0) = precision_umat;
        }
        
        Newton_Raphon(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        EV1 = EV1 + ds_j(0)*Lambdav;
        DEV1 = DEV1 + ds_j(0)*Lambdav;
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV1;
        sigma = el_pred(L0, Eel, ndi);
    }
    
    //Computation of the tangent modulus
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);
    
    vec op = zeros(1);
    mat delta = eye(1,1);
    
    for (int i=0; i<1; i++) {
        if(Ds_j[i] > iota)
            op(i) = 1.;
    }
    
    mat Bbar = zeros(1,1);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }
    
    mat invBbar = zeros(1, 1);
    mat invBhat = zeros(1, 1);
    invBbar = inv(Bbar);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    std::vector<vec> P_epsilon(1);
    P_epsilon[0] = invBhat(0, 0)*(L0*dPhidsigma);
    Lt = L0 - (kappa_j[0]*P_epsilon[0].t());
    
    vec A_v = sigma-L1*EV1;
    double Dgamma_loc = 0.5*sum((A_v_start + A_v)%DEV1);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot) - 0.5*sum((A_v_start + A_v)%DEV1);
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
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
    
} //namespace simcoon

