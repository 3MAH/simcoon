/* This file is part of simcoon.
 
 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file plastic_chaboche.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Linear Kinematical hardening coupled with a power-law hardenig is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>

using namespace std;
using namespace arma;

namespace simcoon {
    
void umat_plasticity_chaboche(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    
    //From the props to the material properties
    double E = props(0);
    double nu= props(1);
    double alpha_iso = props(2);
    double sigmaY = props(3);
    double Q=props(4);
    double b=props(5);
    double C_1 = props(6);
    double D_1 = props(7);
    double C_2 = props(8);
    double D_2 = props(9);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
       
    ///@brief Temperature initialization
    double T_init = statev(0);
    //From the statev to the internal variables
    double p = statev(1);
    vec EP(6);
    EP(0) = statev(2);
    EP(1) = statev(3);
    EP(2) = statev(4);
    EP(3) = statev(5);
    EP(4) = statev(6);
    EP(5) = statev(7);
    
    ///@brief a is the internal variable associated with kinematical hardening
    vec a_1 = zeros(6);
    a_1(0) = statev(8);
    a_1(1) = statev(9);
    a_1(2) = statev(10);
    a_1(3) = statev(11);
    a_1(4) = statev(12);
    a_1(5) = statev(13);

    vec a_2 = zeros(6);
    a_2(0) = statev(14);
    a_2(1) = statev(15);
    a_2(2) = statev(16);
    a_2(3) = statev(17);
    a_2(4) = statev(18);
    a_2(5) = statev(19);
    
    vec X_1 = zeros(6);
    X_1(0) = statev(20);
    X_1(1) = statev(21);
    X_1(2) = statev(22);
    X_1(3) = statev(23);
    X_1(4) = statev(24);
    X_1(5) = statev(25);

    vec X_2 = zeros(6);
    X_2(0) = statev(26);
    X_2(1) = statev(27);
    X_2(2) = statev(28);
    X_2(3) = statev(29);
    X_2(4) = statev(30);
    X_2(5) = statev(31);
    
    // Voce R(p) in closed form from p; statev(32) is an output only
    VoceHardening hardening;
    { int off = 0; hardening.configure(vec{Q, b}, off); }
    double Hp = hardening.R(p);

    // back-strains rotate strain-like, back-stresses stress-like
    EP = rotate_strain(EP, DR);
    a_1 = rotate_strain(a_1, DR);
    a_2 = rotate_strain(a_2, DR);
    X_1 = rotate_stress(X_1, DR);
    X_2 = rotate_stress(X_2, DR);

    vec X = X_1 + X_2;

    //Elstic stiffness tensor
    L = L_iso(E, nu, "Enu");
        
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        stress = vide;
        EP = vide;
        a_1 = vide;
        a_2 = vide;
        X_1 = vide;
        X_2 = vide;
        X = vide;
        p = 0.;
        Hp = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
    //Additional parameters and variables
    double dHpdp = hardening.dR_dp(p);
    
    //Variables values at the start of the increment
    vec stress_start = stress;
    vec EP_start = EP;
    vec a_1start = a_1;
    vec a_2start = a_2;
    vec X_1start = X_1;
    vec X_2start = X_2;
    
    double A_p_start = -Hp;
    vec A_a1_start = -X_1start;
    vec A_a2_start = -X_2start;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    stress = el_pred(L, Eel);
    
    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
    double dPhidp=0.;
    vec dPhida_1 = zeros(6);
    vec dPhida_2 = zeros(6);
    vec dPhidsigma = zeros(6);
    
    //Compute the explicit flow direction
    vec Lambdap = eta_stress(stress-X);
    vec Lambdaa_1 = Lambdap - D_1*a_1;
    vec Lambdaa_2 = Lambdap - D_2*a_2;
    std::vector<vec> kappa_j(1);
    kappa_j[0] = L*Lambdap;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    // Bhat / K below only feed the cutting-plane tangent, not the closest-point branch
    ReturnMappingResult rm;
    const bool use_cpp = (tangent_mode == simcoon::tangent_closest_point);
    if (use_cpp) {
        PlasticityMechanism mech(YieldType::VON_MISES, IsoHardType::VOCE, KinHardType::CHABOCHE, vec{sigmaY, Q, b, C_1, D_1, C_2, D_2}, 1, 2);
        auto &iv = mech.variables();
        iv.get("p").scalar() = p;
        iv.get("EP").raw_voigt() = EP;
        iv.get("a_0").raw_voigt() = a_1;
        iv.get("a_1").raw_voigt() = a_2;
        mech.set_start();
        if (mech.yield_function(stress) > 0.) {   // an admissible trial is committed as is (Lt = L below)
            rm = closest_point_return_mapping(stress, L, {&mech}, {0}, Ds_j, vec{sigmaY});
            if (rm.converged) {
                stress = rm.sigma;
                p = iv.get("p").scalar();
                EP = iv.get("EP").raw_voigt();
                a_1 = iv.get("a_0").raw_voigt();
                a_2 = iv.get("a_1").raw_voigt();
                mat X_i;
                mech.kinematic_hardening().compute_backstresses(iv, X_i);
                X_1 = X_i.col(0);
                X_2 = X_i.col(1);
                X = X_1 + X_2;
                Hp = hardening.R(p);
            }
            else {   // step cut; the kernel state is still the start state (the mechanism owned the iterate)
                tnew_dt = 0.5;
                stress = stress_start;
            }
        }
    }
    else {
    const vec Ir05_v = Ir05();
    const vec Eel_0 = Etot + DEtot - alpha*(T + DT - T_init);
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        p = s_j(0);
        Hp = hardening.R(p);
        dHpdp = hardening.dR_dp(p);
        const vec sigma_eff = stress-X;
        dPhidsigma = eta_stress(sigma_eff);
        dPhidp = -1.*dHpdp;
        dPhida_1 = -1.*(2./3.)*C_1*(dPhidsigma%Ir05_v);
        dPhida_2 = -1.*(2./3.)*C_2*(dPhidsigma%Ir05_v);
        
        //compute Phi and the derivatives
        Phi(0) = Mises_stress(sigma_eff) - Hp - sigmaY;
        
        Lambdap = dPhidsigma;
        Lambdaa_1 = dPhidsigma - D_1*a_1;
        Lambdaa_2 = dPhidsigma - D_2*a_2;
        kappa_j[0] = L*Lambdap;
        
        K(0,0) = dPhidp + sum(dPhida_1%Lambdaa_1) + sum(dPhida_2%Lambdaa_2);
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        EP = EP + ds_j(0)*Lambdap;
        a_1 = a_1 + ds_j(0)*Lambdaa_1;
        a_2 = a_2 + ds_j(0)*Lambdaa_2;
        X_1 += ds_j(0)*(2./3.)*C_1*(Lambdaa_1%Ir05_v);
        X_2 += ds_j(0)*(2./3.)*C_2*(Lambdaa_2%Ir05_v);
        X = X_1 + X_2;
        
        //the stress is now computed using the relationship stress = L(E-Ep)
        Eel = Eel_0 - EP;
        stress = el_pred(L, Eel);
    }
    // the loop reads p at the top and corrects s_j at the bottom
    p = s_j(0);
    Hp = hardening.R(p);
    }
    
    //Computation of the increments of variables
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];
    vec Da_1 = a_1 - a_1start;
    vec Da_2 = a_2 - a_2start;
        
    // tangent via the shared leading-mechanism helper (doc §7.4)
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const ContinuumTangent ct = rm.converged ? cpp_consistent_tangent(rm, L) : compute_tangent_operator(
        tangent_mode, Bhat, kappa_j[0], dPhidsigma, Ds_j, L,
        [&]() -> mat {  // lazy: evaluated only in algorithmic mode
            // J2 flow on sigma-X; the backstress coupling is only in the closest-point branch
            return mat(deta_stress(stress - X));
        });
    Lt = ct.Lt;
    
    double A_p = -Hp;
    vec A_a1 = -X_1;
    vec A_a2 = -X_2;
    
    double Dgamma_loc = 0.5*sum((stress_start+stress)%DEP) + 0.5*(A_p_start + A_p)*Dp + 0.5*sum((A_a1_start + A_a1)%Da_1)+0.5*sum((A_a2_start + A_a2)%Da_2);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%(DEtot-DEP)) - 0.5*sum((A_a1_start + A_a1)%Da_1) - 0.5*sum((A_a2_start + A_a2)%Da_2);
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;
            
    ///@brief statev evolving variables
    //statev
    statev(0) = T_init;
    statev(1) = p;
    
    statev(2) = EP(0);
    statev(3) = EP(1);
    statev(4) = EP(2);
    statev(5) = EP(3);
    statev(6) = EP(4);
    statev(7) = EP(5);
    
    statev(8) = a_1(0);
    statev(9) = a_1(1);
    statev(10) = a_1(2);
    statev(11) = a_1(3);
    statev(12) = a_1(4);
    statev(13) = a_1(5);
    
    statev(14) = a_2(0);
    statev(15) = a_2(1);
    statev(16) = a_2(2);
    statev(17) = a_2(3);
    statev(18) = a_2(4);
    statev(19) = a_2(5);
    
    statev(20) = X_1(0);
    statev(21) = X_1(1);
    statev(22) = X_1(2);
    statev(23) = X_1(3);
    statev(24) = X_1(4);
    statev(25) = X_1(5);
    
    statev(26) = X_2(0);
    statev(27) = X_2(1);
    statev(28) = X_2(2);
    statev(29) = X_2(3);
    statev(30) = X_2(4);
    statev(31) = X_2(5);

    statev(32) = Hp;
}

    
} //namespace simcoon
