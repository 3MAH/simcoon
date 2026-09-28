///@file Zener_fast.hpp
///@brief User subroutine for Zener viscoelastic model in 3D case, with thermoelastic effect(Poynting-Thomson model)
///@brief This implementation uses a single scalar internal variable for the evaluation of the viscoelastic strain increment
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
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Viscoelasticity/linear_viscoelastic.hpp>

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
    
void umat_zener_fast(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    vec stress_start = stress;
    
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
    L = L0;
    mat L1 = L_iso(E1, nu1, "Enu");
    
    mat H1 = H_iso(etaB1, etaS1);                  //dimension of stiffness tensor
    
    if(start) { //Initialization
        T_init = T;
        EV1 = zeros(6);
        stress = zeros(6);
        stress_start = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
    //Variables at the start of the increment
    vec EV1_start = EV1;
    vec A_v_start = stress_start - L1*EV1;

    // Implicit (backward-Euler) step of the Kelvin branch in closed form: the exact solution of
    // the discrete equations and its consistent tangent (linear_viscoelastic.hpp)
    const LinearViscoStep st = kelvin_series_step(L0, {L1}, {H1}, {EV1_start},
                                                  Etot + DEtot - alpha*(T + DT - T_init), alpha, DTime);
    EV1 = st.EV_i[0];
    const vec DEV1 = EV1 - EV1_start;
    v += norm_strain(DEV1);
    const vec Eel = Etot + DEtot - alpha*(T + DT - T_init) - EV1;
    stress = el_pred(L0, Eel, ndi);
    Lt = (tangent_mode == tangent_none) ? L0 : st.dSdE;

    vec A_v = stress-L1*EV1;
    double Dgamma_loc = 0.5*sum((A_v_start + A_v)%DEV1);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%DEtot) - 0.5*sum((A_v_start + A_v)%DEV1);
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

