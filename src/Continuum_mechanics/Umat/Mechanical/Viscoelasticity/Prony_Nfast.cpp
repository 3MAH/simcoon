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
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>

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
    
void umat_prony_Nfast(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
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
    
    //From the statev to the internal variables (EV_tilde, statev(1..6), is an output only)
    std::vector<vec> EV_i(N_prony);
    vec v = zeros(N_prony);
    for (int i=0; i<N_prony; i++) {
        v(i) = statev(i*7+7);
        EV_i[i] = rotate_strain(vec(statev.subvec(i*7+8, i*7+13)), DR);
    }

    std::vector<mat> L_i(N_prony);
    std::vector<mat> H_i(N_prony);
    for (int i=0; i<N_prony; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
    }

    vec stress_start = stress;
    if(start) { //Initialization
        T_init = T;
        for (int i=0; i<N_prony; i++) {
            EV_i[i] = zeros(6);
        }
        stress = zeros(6);
        stress_start = zeros(6);

        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }

    // Implicit (backward-Euler) step of the Maxwell branches in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const std::vector<vec> EV_i_start = EV_i;
    const vec eps_e_start = Etot - alpha*(T - T_init);
    const vec eps_e = Etot + DEtot - alpha*(T + DT - T_init);
    const LinearViscoStep st = maxwell_parallel_step(L0, L_i, H_i, EV_i_start, eps_e, alpha, DTime);
    vec EV_tilde = zeros(6);
    double Dgamma_loc = 0.;
    for (int i=0; i<N_prony; i++) {
        EV_i[i] = st.EV_i[i];
        const vec DEV_i = EV_i[i] - EV_i_start[i];
        v(i) += norm_strain(DEV_i);
        EV_tilde += (M0*L_i[i])*EV_i[i];
        // dashpot dissipation on the branch forces L_i (eps_e - EV_i), trapezoidal
        Dgamma_loc += 0.5*sum((L_i[i]*(eps_e_start - EV_i_start[i]) + L_i[i]*(eps_e - EV_i[i]))%DEV_i);
    }
    stress = el_pred(L0, eps_e - EV_tilde, ndi);
    Lt = (tangent_mode == tangent_none) ? L0 : st.dSdE;

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%DEtot) - Dgamma_loc;
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

