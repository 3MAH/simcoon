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
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
using namespace std;
using namespace arma;

///@brief The viscoelastic burger model requires 4+N*4 constants:
//      -------------------
//
///@brief
///@brief      props(0) = E0                   - Thermoelastic Young's modulus
///@brief      props(1) = nu0                  - Thermoelastic Poisson's ratio
///@brief      props(2) = alpha_iso            - Thermoelastic CTE
///@brief      props(3) = N_kelvin             - Number of Kelvin branches
///@brief      props(4+i*4) = E_visco(i)       - Viscoelastic Young modulus of Zener branch i
///@brief      props(4+i*4+1) = nu_visco(i)    - Viscoelastic Poisson ratio of Zener branch i
///@brief      props(4+i*4+2) = etaB_visco     - Viscoelastic Bulk viscosity of Zener branch i
///@brief      props(4+i*4+3) = etaS_visco     - Viscoelastic Bulk viscosity of Zener branch i

///@brief Number of statev required for thermoelastic constitutive law : 7+N*7

namespace simcoon {
    
void umat_zener_Nfast_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
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
    EV(0) = statev(1);
    EV(1) = statev(2);
    EV(2) = statev(3);
    EV(3) = statev(4);
    EV(4) = statev(5);
    EV(5) = statev(6);
    
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
    EV = rotate_strain(EV, DR);
    for (int i=0; i<N_kelvin; i++) {
        EV_i[i] = rotate_strain(EV_i[i], DR);
    }
    
    std::vector<mat> L_i(N_kelvin);
    std::vector<mat> H_i(N_kelvin);
    std::vector<mat> invH_i(N_kelvin);
    
    for (int i=0; i<N_kelvin; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
        invH_i[i] = inv(H_i[i]);
    }
    
    vec sigma_start = sigma;
    std::vector<vec> DEV_i(N_kelvin);
    std::vector<vec> A_v(N_kelvin);
    std::vector<vec> A_v_start(N_kelvin);
    
    if(start) { //Initialization
        T_init = T;
        EV = zeros(6);
        for (int i=0; i<N_kelvin; i++) {
            EV_i[i] = zeros(6);
            DEV_i[i] = zeros(6);
            A_v[i] = zeros(6);
            A_v_start[i] = zeros(6);
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
    
    //Variables at the start of the increment
    vec DEV = zeros(6);
    vec EV_start = EV;
    for (int i=0; i<N_kelvin; i++) {
        A_v_start[i] += L_i[i]*EV_i[i];
    }
    
    //Variables required for the loop
    vec s_j = zeros(N_kelvin);
    s_j = v;
    vec Ds_j = zeros(N_kelvin);
    vec ds_j = zeros(N_kelvin);
    
    //Determination of the initial, predicted stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EV;
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
    vec Phi = zeros(N_kelvin);
    mat B = zeros(N_kelvin,N_kelvin);
    vec Y_crit = zeros(N_kelvin);
    
    vec dPhidv = zeros(N_kelvin);
    std::vector<vec> dPhidEv(N_kelvin);
    std::vector<vec> dPhi_idsigma(N_kelvin);
    for (int i=0; i<N_kelvin; i++) {
        dPhidEv[i] = zeros(6);
        dPhi_idsigma[i] = zeros(6);
    }
    
    //Compute the explicit flow direction
    std::vector<vec> flow_visco(N_kelvin);
    std::vector<vec> Lambdav(N_kelvin);
    std::vector<vec> kappa_j(N_kelvin);
    for (int i=0; i<N_kelvin; i++) {
        flow_visco[i] = invH_i[i]*(sigma-L_i[i]*EV_i[i]);
        Lambdav[i] = eta_norm_strain(flow_visco[i]);
        kappa_j[i] = L0*Lambdav[i];
    }
    
    mat K = zeros(N_kelvin,N_kelvin);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        v = s_j;

        for (int i=0; i<N_kelvin; i++) {
            flow_visco[i] = invH_i[i]*(sigma-L_i[i]*EV_i[i]);
            Lambdav[i] = eta_norm_strain(flow_visco[i]);
            dPhi_idsigma[i] = invH_i[i]*(eta_norm_strain(flow_visco[i])%Ir05()); //Dimension of strain (The flow is of stress type here)
            kappa_j[i] = L0*Lambdav[i];
            
            if (DTime > iota) {
                Phi(i) = norm_strain(flow_visco[i]) - Ds_j(i)/DTime;
                dPhidv[i] = -1.*sum((dPhi_idsigma[i])%(L_i[i]*Lambdav[i]))-1./DTime;
            }
            else {
                Phi(i) = norm_strain(flow_visco[i]);
                dPhidv[i] = -1.*sum((dPhi_idsigma[i])%(L_i[i]*Lambdav[i]));
            }
            kappa_j[i] = L0*Lambdav[i];
            K(i,i) = dPhidv[i];
        }
        
        B = zeros(N_kelvin,N_kelvin);
        for (int i=0; i<N_kelvin; i++) {
            B(i, i) = K(i,i);
            Y_crit(i) = norm_strain(flow_visco[i]);
            if (Y_crit(i) < precision_umat) {
                Y_crit(i) = precision_umat;
            }
            for (int j=0; j<N_kelvin; j++) {
                B(i,j) += -1.*sum(dPhi_idsigma[i]%kappa_j[j]);
            }
        }
        
        Newton_Raphon(Phi, Y_crit, B, Ds_j, ds_j, error);

        EV = zeros(6);
        for (int i=0; i<N_kelvin; i++) {
            s_j(i) += ds_j(i);
            DEV += ds_j(i)*Lambdav[i];
            EV_i[i] += ds_j(i)*Lambdav[i];
            DEV_i[i] += ds_j(i)*Lambdav[i];
            EV += EV_i[i];
        }
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV;
        DEel = DEtot - alpha*(DT) - DEV;
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
    
    mat Bhat = zeros(N_kelvin, N_kelvin);
    
    vec op = zeros(N_kelvin);
    mat delta = eye(N_kelvin,N_kelvin);
    mat Bbar = zeros(N_kelvin,N_kelvin);
    mat invBbar = zeros(N_kelvin, N_kelvin);
    mat invBhat = zeros(N_kelvin, N_kelvin);
    std::vector<vec> P_epsilon(N_kelvin);
    std::vector<double> P_theta(N_kelvin);
    
    dSdE = L0;
    dSdT = -1.*L0*alpha;    
    
    for (int i=0; i<N_kelvin; i++) {
        P_epsilon[i] = zeros(6);
        P_theta[i] = 0.;        
    }
    
    for (int i=0; i<N_kelvin; i++) {
        
        if(Ds_j(i) > iota)
            op(i) = 1.;
        
        for (int j = 0; j <N_kelvin; j++) {
            Bhat(i, j) = sum(dPhi_idsigma[i]%kappa_j[j]) - K(i,j);
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }

    invBbar = inv(Bbar);
    
    for (int i=0; i<N_kelvin; i++) {
        for (int j = 0; j <N_kelvin; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    for (int i=0; i<N_kelvin; i++) {
        for (int j = 0; j <N_kelvin; j++) {
            P_epsilon[i] += invBhat(j, i)*(L0*dPhi_idsigma[j]);
            P_theta[i] += invBhat(j, i)*(dPhi_idsigma[j]*(L0*alpha));
            
        }
        dSdE += -1.*(kappa_j[i]*P_epsilon[i].t());
        dSdT +=  -1.*(kappa_j[i]*P_theta[i]);        
    }
    
    for (int i=0; i<N_kelvin; i++) {
        A_v[i] += L_i[i]*EV_i[i];
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
    
    vec Gamma_epsilon = zeros(6);
    double Gamma_theta = 0.;
    
    vec N_epsilon = zeros(6);
    double N_theta = 0.;
    
    double A_p = -Hp;
    double dA_pdp = -dHpdp;
    //    double A_theta = 0;
    
    if(DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        for (int i=0; i<N_prony; i++) {
            Gamma_epsilon += (dSdE*DEV)*(1./DTime) + sum(sigma%Lambdav)*P_epsilon[0]/DTime;
            Gamma_theta = sum(dSdT%DEV)*(1./DTime) + sum(sigma%Lambdav)*P_theta[0]/DTime;
        }
        
        N_epsilon = -1./DTime*(T + DT)*(dSdE*alpha);
        N_theta = -1./DTime*(T + DT)*sum(dSdT%alpha) -1.*Deta/DTime - rho*c_p*(1./DTime);
        
        drdE = N_epsilon + Gamma_epsilon;
        drdT = N_theta + Gamma_theta;
        
        r = sum(N_epsilon%DEtot) + N_theta*DT + sum(Gamma_epsilon%DEtot) + Gamma_theta*DT;
    }
    
    double Dgamma_loc = 0.;
    for (int i=0; i<N_kelvin; i++) {
        Dgamma_loc += 0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
    }
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEV));
    for (int i=0; i<N_kelvin; i++) {
        Wm_r += -0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
    }
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir = (T+0.5*DT)*Deta_ir;
        
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

