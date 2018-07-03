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
using namespace std;
using namespace arma;

///@brief The viscoelastic Prony series model requires 4+N*4 constants:
//      -------------------
///@brief      props(0) = E0                - Thermoelastic Young's modulus
///@brief      props(1) = nu0               - Thermoelastic Poisson's ratio
///@brief      props(2) = alpha_iso         - Thermoelastic CTE
///@brief      props(3) = N_prony           - Number of Prony series
///@brief      props(4+i*4) = E_visco(i)    - Viscoelastic Young modulus of Prony branch i
///@brief      props(4+i*4+1) = nu_visco(i) - Viscoelastic Poisson ratio of Prony branch i
///@brief      props(4+i*4+2) = etaB_visco  - Viscoelastic Bulk viscosity of Prony branch i
///@brief      props(4+i*4+3) = etaS_visco  - Viscoelastic Bulk viscosity of Prony branch i

///@brief Number of statev required for thermoelastic constitutive law : 7+N*7

namespace simcoon {
    
void umat_prony_Nfast_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
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
    
    std::vector<mat> L_i(N_prony);
    std::vector<mat> H_i(N_prony);
    std::vector<mat> invH_i(N_prony);
    
    for (int i=0; i<N_prony; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
        invH_i[i] = inv(H_i[i]);
    }
    
    vec sigma_start = sigma;
    std::vector<vec> DEV_i(N_prony);
    std::vector<vec> A_v(N_prony);
    std::vector<vec> A_v_start(N_prony);
    
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
    vec DEV_tilde = zeros(6);
    vec EV_tilde_start = EV_tilde;
    for (int i=0; i<N_prony; i++) {
        A_v_start[i] += L_i[i]*(Etot - alpha*(T+DT-T_init) - EV_i[i]);
    }
    
    //Variables required for the loop
    vec s_j = v;
    vec Ds_j = zeros(N_prony);
    vec ds_j = zeros(N_prony);
    
    //Determination of the initial, predicted stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EV_tilde;
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
    vec Phi = zeros(N_prony);
    mat B = zeros(N_prony,N_prony);
    vec Y_crit = zeros(N_prony);
    
    vec dPhidv = zeros(N_prony);
    std::vector<vec> dPhidEv(N_prony);
    std::vector<vec> dPhi_idv_temp(N_prony);
    for (int i=0; i<N_prony; i++) {
        dPhidEv[i] = zeros(6);
        dPhi_idv_temp[i] = zeros(6);
    }
    
    //Compute the explicit flow direction
    std::vector<vec> flow_visco(N_prony);
    std::vector<vec> Lambdav(N_prony);
    std::vector<vec> kappa_j(N_prony);
    for (int i=0; i<N_prony; i++) {
        flow_visco[i] = invH_i[i]*(L_i[i]*(Etot+DEtot-EV_i[i]));
        Lambdav[i] = eta_norm_strain(flow_visco[i]);
        kappa_j[i] = L_i[i]*Lambdav[i];
    }
    
    mat K = zeros(N_prony,N_prony);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        v = s_j;

        for (int i=0; i<N_prony; i++) {
            flow_visco[i] = invH_i[i]*(L_i[i]*(Etot+DEtot)-L_i[i]*EV_i[i]);
            Lambdav[i] = eta_norm_strain(flow_visco[i]);
            dPhi_idv_temp[i] = invH_i[i]*(eta_norm_strain(flow_visco[i])%Ir05()); //Dimension of strain (The flow is of stress type here)
            kappa_j[i] = L_i[i]*Lambdav[i];
            
            if (DTime > iota) {
                Phi(i) = norm_strain(flow_visco[i]) - Ds_j(i)/DTime;
                dPhidv[i] = -1.*sum((dPhi_idv_temp[i])%(L_i[i]*Lambdav[i]))-1./DTime;
            }
            else {
                Phi(i) = norm_strain(flow_visco[i]);
                dPhidv[i] = -1.*sum((dPhi_idv_temp[i])%(L_i[i]*Lambdav[i]));
            }
            kappa_j[i] = L_i[i]*Lambdav[i];
            K(i,i) = dPhidv[i];
        }
        
        B = zeros(N_prony,N_prony);
        for (int i=0; i<N_prony; i++) {
            B(i, i) = K(i,i);
            Y_crit(i) = norm_strain(flow_visco[i]);
            if (Y_crit(i) < precision_umat) {
                Y_crit(i) = precision_umat;
            }
        }
        
        Newton_Raphon(Phi, Y_crit, B, Ds_j, ds_j, error);

        EV_tilde = zeros(6);
        for (int i=0; i<N_prony; i++) {
            s_j(i) += ds_j(i);
            DEV_tilde += (M0*L_i[i])*(ds_j(i)*Lambdav[i]);
            EV_i[i] += ds_j(i)*Lambdav[i];
            DEV_i[i] += ds_j(i)*Lambdav[i];
            EV_tilde += (M0*L_i[i])*EV_i[i];
        }
        
        //the stress is now computed using the relationship sigma = L0 E-sum LpEp
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV_tilde;
        DEel = DEtot - alpha*(DT) - DEV_tilde;
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
    
    // Tangent modulus for prony series
    // L0 - summation( (L_i[i]*Lambdav[i]) \dyad (dPhi_idv_temp[i]*Lambdav[i])/A[i] )
    // where A[i]= K(i,i)
                          
    mat Bhat = zeros(N_prony, N_prony);
    
    vec op = zeros(N_prony);
    mat delta = eye(N_prony,N_prony);
    mat Bbar = zeros(N_prony,N_prony);
    mat invBbar = zeros(N_prony, N_prony);
    mat invBhat = zeros(N_prony, N_prony);
    std::vector<vec> P_epsilon(N_prony);
    std::vector<double> P_theta(N_prony);
    dSdE = L0;
    dSdT = -1.*L0*alpha;
    
    for (int i=0; i<N_prony; i++) {
        P_epsilon[i] = zeros(6);
        P_theta[i] = 0.;
    }
    
    for (int i=0; i<N_prony; i++) {
        
        if(Ds_j(i) > iota)
            op(i) = 1.;
        
        for (int j = 0; j <N_prony; j++) {
            Bhat(i, j) = - K(i,j);
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }

    invBbar = inv(Bbar);
    
    for (int i=0; i<N_prony; i++) {
        for (int j = 0; j <N_prony; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    for (int i=0; i<N_prony; i++) {
        for (int j = 0; j <N_prony; j++) {
            P_epsilon[i] += invBhat(j, i)*(L_i[j]*dPhi_idv_temp[j]);
            P_theta[i] += invBhat(j, i)*sum(dPhi_idv_temp[j]%(L_i[j]*alpha));
        }
        dSdE += -1.*(kappa_j[i]*P_epsilon[i].t());
        dSdT +=  -1.*(kappa_j[i]*P_theta[i]);
    }
    
    for (int i=0; i<N_prony; i++) {
        A_v[i] += L_i[i]*(Etot - alpha*(T+DT-T_init) - EV_i[i]);
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
    for (int i=0; i<N_prony; i++) {
        Dgamma_loc += 0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
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
    Wt_ir = (T+0.5*DT)*Deta_ir;
    
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

