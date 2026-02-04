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
    
void umat_zener_Nfast(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double E0 = props(0);
    double nu0 = props(1);
    double alpha_iso = props(2);
    int N_kelvin = int(props(3));
    
    vec E_visco = zeros(N_kelvin);
    vec nu_visco = zeros(N_kelvin);
    vec etaB_visco = zeros(N_kelvin);
    vec etaS_visco = zeros(N_kelvin);
    
    for (int i=0; i<N_kelvin; i++) {
        E_visco(i) = props(4+i*4);
        nu_visco(i) = props(4+i*4+1);
        etaB_visco(i) = props(4+i*4+2);
        etaS_visco(i) = props(4+i*4+3);
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
    
    vec sigma_start = sigma;
    std::vector<vec> DEV_i(N_kelvin);
    std::vector<vec> A_v(N_kelvin);
    std::vector<vec> A_v_start(N_kelvin);
    
    if(start) { //Initialization
        T_init = T;
        EV = zeros(6);
        for (int i=0; i<N_kelvin; i++) {
            EV_i[i] = zeros(6);
        }
        sigma = zeros(6);
        sigma_start = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
    for (int i=0; i<N_kelvin; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
        invH_i[i] = inv(H_i[i]);
        
        DEV_i[i] = zeros(6);
        A_v[i] = zeros(6);
        A_v_start[i] = zeros(6);
    }
    
    //Variables at the start of the increment
    vec DEV = zeros(6);
    vec EV_start = EV;
    for (int i=0; i<N_kelvin; i++) {
        A_v_start[i] = sigma_start - L_i[i]*EV_i[i];
    }
    
    //Variables required for the loop
    vec s_j = zeros(N_kelvin);
    s_j = v;
    vec Ds_j = zeros(N_kelvin);
    vec ds_j = zeros(N_kelvin);
    
    //Determination of the initial, predicted stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EV;
    sigma = el_pred(L0, Eel, ndi);

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
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        v = s_j;

        for (int i=0; i<N_kelvin; i++) {
            flow_visco[i] = invH_i[i]*(sigma-L_i[i]*EV_i[i]);
            Lambdav[i] = eta_norm_strain(flow_visco[i]);
            dPhi_idsigma[i] = invH_i[i]*(eta_norm_strain(flow_visco[i])%Ir05()); //Dimension of strain (The flow is of stress type here)
            
            if (DTime > simcoon::iota) {
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
            if (Y_crit(i) < simcoon::precision_umat) {
                Y_crit(i) = simcoon::precision_umat;
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
        sigma = el_pred(L0, Eel, ndi);
    }
    
    mat Bhat = zeros(N_kelvin, N_kelvin);
    
    vec op = zeros(N_kelvin);
    mat delta = eye(N_kelvin,N_kelvin);
    mat Bbar = zeros(N_kelvin,N_kelvin);
    mat invBbar = zeros(N_kelvin, N_kelvin);
    mat invBhat = zeros(N_kelvin, N_kelvin);
    std::vector<vec> P_epsilon(N_kelvin);
    
    Lt = L0;    
    for (int i=0; i<N_kelvin; i++) {
        P_epsilon[i] = zeros(6);
    }
    
    for (int i=0; i<N_kelvin; i++) {
        
        if(Ds_j(i) > simcoon::iota)
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
        }
        Lt += -1.*(kappa_j[i]*P_epsilon[i].t());
    }
    
    for (int i=0; i<N_kelvin; i++) {
        A_v[i] = sigma - L_i[i]*EV_i[i];
    }
    double Dgamma_loc = 0.;
    for (int i=0; i<N_kelvin; i++) {
        Dgamma_loc += 0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
    }
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    for (int i=0; i<N_kelvin; i++) {
        Wm_r += -0.5*sum((A_v_start[i] + A_v[i])%DEV_i[i]);
    }
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
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

} //namespace simcoon

