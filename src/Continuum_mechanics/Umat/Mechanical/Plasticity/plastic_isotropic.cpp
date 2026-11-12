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

///@file plastic_isotropic.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Isotropic hardening with a power-law hardenig is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>


using namespace std;
using namespace arma;

namespace simcoon {
    
void umat_plasticity_iso(const string &umat_name, const vec &Etot, const vec &DEtot, vec &stress, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode, const int &ncp_code)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    
    const NCPMethod ncp = ncp_method_of(ncp_code, "EPICP");

    //From the props to the material properties
    double E = props(0);
    double nu= props(1);
    double alpha_iso = props(2);
    double sigmaY = props(3);
    double k=props(4);
    double m=props(5);
    
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
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    
    //Elstic stiffness tensor
    L = L_iso(E, nu, "Enu");
    
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        stress = vide;
        EP = vide;
        p = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;            
    }
    
    double Hp=0.;
    double dHpdp=0.;

    // PowerLawHardening: C1 blend below p = 1e-6 (slope of k p^m infinite at onset for m < 1)
    PowerLawHardening hardening;
    int hardening_offset = 0;
    hardening.configure(vec{k, m}, hardening_offset);
    Hp = hardening.R(p);
    dHpdp = hardening.dR_dp(p);
    
    //Variables values at the start of the increment
    vec stress_start = stress;
    vec EP_start = EP;
    double A_p_start = -Hp;
    
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
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec Lambdap = eta_stress(stress);
    std::vector<vec> kappa_j(1);
    mul6(L, Lambdap, kappa_j[0]);
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    // Bhat / K below only feed the cutting-plane tangent, not the closest-point branch
    ReturnMappingResult rm;
    const bool use_cpp = (tangent_mode == simcoon::tangent_closest_point);
    if (use_cpp) {
        PlasticityMechanism mech(YieldType::VON_MISES, IsoHardType::POWER_LAW, KinHardType::NONE, vec{sigmaY, k, m});
        auto &iv = mech.variables();
        iv.get("p").scalar() = p;
        iv.get("EP").raw_voigt() = EP;
        mech.set_start();
        if (mech.yield_function(stress) > 0.) {   // an admissible trial is committed as is (Lt = L below)
            rm = closest_point_return_mapping(stress, L, {&mech}, {0}, Ds_j, vec{sigmaY});
            if (rm.converged) {
                stress = rm.sigma;
                p = iv.get("p").scalar();
                EP = iv.get("EP").raw_voigt();
                Hp = hardening.R(p);
            }
            else {   // step cut; the kernel state is still the start state (the mechanism owned the iterate)
                tnew_dt = 0.5;
                stress = stress_start;
            }
        }
    }
    else {
    ncp_step local_step(ncp);
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        p = s_j(0);
        Hp = hardening.R(p);
        dHpdp = hardening.dR_dp(p);
        dPhidsigma = eta_stress(stress);
        dPhidp = -1.*dHpdp;
        
        //compute Phi and the derivatives
        Phi(0) = Mises_stress(stress) - Hp - sigmaY;
        
        Lambdap = eta_stress(stress);
        mul6(L, Lambdap, kappa_j[0]);
        
        K(0,0) = dPhidp;
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;

        if (compteur == 0 && Phi(0) < 0. && sqrt(pow(Phi(0),2.)) + Phi(0) == 0.) {
            // admissible trial: FB residual exactly 0, skip the solve
            ds_j.zeros();
            error = 0.;
        }
        else {
            local_step(Phi, Y_crit, B, Ds_j, ds_j, error);
        }
        
        s_j(0) += ds_j(0);
        EP = EP + ds_j(0)*Lambdap;
        
        //the stress is now computed using the relationship stress = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        stress = el_pred(L, Eel);
    }
    }
    
    //Computation of the increments of variables
    vec Dsigma = stress - stress_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];
    
    // tangent via the shared leading-mechanism helper (doc §7.4)
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    // lazy: needed in algorithmic mode with an active multiplier only
    HessianProvider dLambda_dsigma = nullptr;
    if (Ds_j(0) > simcoon::iota) {
        dLambda_dsigma = [&]() -> mat {
            // J2 associated flow: dLambda/dsigma = deta_stress, complete (no internal-variable coupling)
            return mat(deta_stress(stress));
        };
    }
    if (!use_cpp && tangent_mode == simcoon::tangent_algorithmic && Ds_j(0) > simcoon::iota) {
        // closed-form radial-return tangent (formula in the header)
        const double G = E/(2.*(1.+nu));
        const double K_b = E/(3.*(1.-2.*nu));
        const double q_tr = Mises_stress(stress) + 3.*G*Dp;
        const double theta = 1. - 3.*G*Dp/q_tr;
        const double theta_bar = 1./(1. + hardening.dR_dp(p)/(3.*G)) - (1. - theta);
        const vec s_dev = dev(stress);
        const double s_norm = sqrt(s_dev(0)*s_dev(0) + s_dev(1)*s_dev(1) + s_dev(2)*s_dev(2)
                                   + 2.*(s_dev(3)*s_dev(3) + s_dev(4)*s_dev(4) + s_dev(5)*s_dev(5)));
        const vec N = s_dev/s_norm;
        Lt = 3.*K_b*Ivol() + 2.*G*theta*Idev() - 2.*G*theta_bar*(N*N.t());
    }
    else {
        const ContinuumTangent ct = rm.converged ? cpp_consistent_tangent(rm, L) : compute_tangent_operator(
            tangent_mode, Bhat, kappa_j[0], dPhidsigma, Ds_j, L, dLambda_dsigma);
        Lt = ct.Lt;
    }

    std::vector<double> P_theta(1);
    vec L_alpha;
    mul6(L, alpha, L_alpha);
    P_theta[0] = dPhidtheta - sum(dPhidsigma%L_alpha);

    double A_p = -Hp;        
    double Dgamma_loc = 0.5*sum((stress_start+stress)%DEP) + 0.5*(A_p_start + A_p)*Dp;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((stress_start+stress)%DEtot);
    Wm_r += 0.5*sum((stress_start+stress)%(DEtot-DEP));
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
}
    
} //namespace simcoon
