///@file Zener_fast.hpp
///@brief User subroutine for Zener viscoelastic model with  N viscoelastic Kelvin branches in series in 3D case, with thermoelastic effect
///@brief This implementation uses a single scalar internal variable for the evaluation of the viscoelastic strain increment
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
    
void umat_zener_Nfast(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
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
    L = L0;
    ///@brief Temperature initialization
    double T_init = statev(0);
    
    //From the statev to the internal variables (EV, statev(1..6), is an output only)
    std::vector<vec> EV_i(N_kelvin);
    vec v = zeros(N_kelvin);
    for (int i=0; i<N_kelvin; i++) {
        v(i) = statev(i*7+7);
        EV_i[i] = rotate_strain(vec(statev.subvec(i*7+8, i*7+13)), DR);
    }

    std::vector<mat> L_i(N_kelvin);
    std::vector<mat> H_i(N_kelvin);
    for (int i=0; i<N_kelvin; i++) {
        L_i[i] = L_iso(E_visco(i), nu_visco(i), "Enu");
        H_i[i] = H_iso(etaB_visco(i), etaS_visco(i));
    }

    vec stress_start = stress;
    if(start) { //Initialization
        T_init = T;
        for (int i=0; i<N_kelvin; i++) {
            EV_i[i] = zeros(6);
        }
        stress = zeros(6);
        stress_start = zeros(6);

        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }

    // Implicit (backward-Euler) step of the Kelvin branches in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const std::vector<vec> EV_i_start = EV_i;
    const vec eps_e = Etot + DEtot - alpha*(T + DT - T_init);
    const LinearViscoStep st = kelvin_series_step(L0, L_i, H_i, EV_i_start, eps_e, alpha, DTime);
    vec EV = zeros(6);
    for (int i=0; i<N_kelvin; i++) {
        EV_i[i] = st.EV_i[i];
        v(i) += norm_strain(EV_i[i] - EV_i_start[i]);
        EV += EV_i[i];
    }
    stress = el_pred(L0, eps_e - EV, ndi);
    Lt = (tangent_mode == tangent_none) ? L0 : st.dSdE;

    // dashpot dissipation on the viscous stresses sigma - L_i EV_i, trapezoidal
    double Dgamma_loc = 0.;
    for (int i=0; i<N_kelvin; i++) {
        Dgamma_loc += 0.5*sum(((stress_start - L_i[i]*EV_i_start[i]) + (stress - L_i[i]*EV_i[i]))
                              %(EV_i[i] - EV_i_start[i]));
    }

    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%DEtot) - Dgamma_loc;
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

