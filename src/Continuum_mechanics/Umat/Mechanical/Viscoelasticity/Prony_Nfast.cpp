///@file Prony.hpp
///@brief User subroutine for Prony series viscoelastic model in 3D case
///@author Chemisky, Chatzigeorgiou
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
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
    
void umat_prony_Nfast(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double E0 = props(0);
    double nu0 = props(1);
    double alpha_iso = props(2);
    int N_prony = int(props(3));
    
    vec E_visco = zeros(N_prony);
    vec nu_visco = zeros(N_prony);
    vec etaB_visco = zeros(N_prony);
    vec etaS_visco = zeros(N_prony);
    
    for (int i=0; i<N_prony; i++) {
        E_visco(i) = props(4+i*4);
        nu_visco(i) = props(4+i*4+1);
        etaB_visco(i) = props(4+i*4+2);
        etaS_visco(i) = props(4+i*4+3);
    }
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Define the viscoelastic stiffness
    mat L0 = L_iso(E0, nu0, "Enu");
    L = L0;
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
    }
    
    for (int i=0; i<N_prony; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
        invH_i[i] = inv(H_i[i]);
        
        DEV_i[i] = zeros(6);
        A_v[i] = zeros(6);
        A_v_start[i] = zeros(6);
    }
    
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
    sigma = el_pred(L0, Eel, ndi);

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
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        v = s_j;

        for (int i=0; i<N_prony; i++) {
            flow_visco[i] = invH_i[i]*(L_i[i]*(Etot+DEtot)-L_i[i]*EV_i[i]);
            Lambdav[i] = eta_norm_strain(flow_visco[i]);
            dPhi_idv_temp[i] = invH_i[i]*(eta_norm_strain(flow_visco[i])%Ir05()); //Dimension of strain (The flow is of stress type here)
            
            if (DTime > simcoon::iota) {
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
            if (Y_crit(i) < simcoon::precision_umat) {
                Y_crit(i) = simcoon::precision_umat;
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
        sigma = el_pred(L0, Eel, ndi);
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
    
    Lt = L0;
    for (int i=0; i<N_prony; i++) {
        P_epsilon[i] = zeros(6);
    }
    
    for (int i=0; i<N_prony; i++) {
        
        if(Ds_j(i) > simcoon::iota)
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
        }
        Lt += -1.*(kappa_j[i]*P_epsilon[i].t());                     
    }
    
    for (int i=0; i<N_prony; i++) {
        A_v[i] += L_i[i]*(Etot - alpha*(T+DT-T_init) - EV_i[i]);
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

} //namespace simcoon

