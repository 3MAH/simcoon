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
#include <armadillo>
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/plasticity_mechanism.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void umat_plasticity_iso_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    
	//From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double E = props(2);
	double nu= props(3);
	double alpha_iso = props(4);
	double sigmaY = props(5);
	double k=props(6);
	double m=props(7);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Elastic stiffness tensor
    mat L = L_iso(E, nu, "Enu");
    
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
    
	///@brief Initialization
	if(start)
	{
		T_init = T;
		vec vide = zeros(6);
		sigma = vide;
		EP = vide;
		p = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
        
        Wt = 0.;
        Wt_r = 0.;
        Wt_ir = 0.;
	}
    
    //Additional parameters
    double c_0 = rho*c_p;
    
    double Hp=0.;
    double dHpdp=0.;

    // PowerLawHardening: C1 blend below p = 1e-6 (slope of k p^m infinite at onset for m < 1)
    PowerLawHardening hardening;
    int hardening_offset = 0;
    hardening.configure(vec{k, m}, hardening_offset);
    Hp = hardening.R(p);
    dHpdp = hardening.dR_dp(p);
    
    //Variables values at the start of the increment
    vec sigma_start = sigma;
    vec EP_start = EP;
    double A_p_start = -Hp;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);
    
	///Elastic prediction - Accounting for the thermal prediction
	vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    sigma = el_pred(L, Eel);
    
	//Define the plastic function and the stress
	vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
	double dPhidp=0.;
	vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec Lambdap = eta_stress(sigma);
    std::vector<vec> kappa_j(1);
    mul6(L, Lambdap, kappa_j[0]);
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    // the thermal tail reads dPhidsigma / Lambdap / kappa_j / dHpdp from either branch
    ReturnMappingResult rm;
    const bool use_cpp = (tangent_mode == simcoon::tangent_closest_point);
    if (use_cpp) {
        PlasticityMechanism mech(YieldType::VON_MISES, IsoHardType::POWER_LAW, KinHardType::NONE, vec{sigmaY, k, m});
        auto &iv = mech.variables();
        iv.get("p").scalar() = p;
        iv.get("EP").raw_voigt() = EP;
        mech.set_start();
        if (mech.yield_function(sigma) > 0.) {   // an admissible trial is committed as is (Lt = L below)
            rm = closest_point_return_mapping(sigma, L, {&mech}, {0}, Ds_j, vec{sigmaY});
            if (rm.converged) {
                sigma = rm.sigma;
                p = iv.get("p").scalar();
                EP = iv.get("EP").raw_voigt();
                Hp = hardening.R(p);
                dHpdp = hardening.dR_dp(p);
                const ClosestPointIngredients &ing = (*mech.closest_point_ingredients())[0];
                dPhidsigma = ing.dPhi_dsigma;
                Lambdap = ing.Lambda;
                kappa_j = rm.kappa_j;
            }
            else {   // step cut; the kernel state is still the start state (the mechanism owned the iterate)
                tnew_dt = 0.5;
                sigma = sigma_start;
            }
        }
    }
    else {
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        p = s_j(0);
        Hp = hardening.R(p);
        dHpdp = hardening.dR_dp(p);
        dPhidsigma = eta_stress(sigma);
        dPhidp = -1.*dHpdp;
        
        //compute Phi and the derivatives
        Phi(0) = Mises_stress(sigma) - Hp - sigmaY;

        Lambdap = dPhidsigma;   // associated J2 flow
        mul6(L, Lambdap, kappa_j[0]);
        
        K(0,0) = dPhidp;
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);

        s_j(0) += ds_j(0);
        EP = EP + ds_j(0)*Lambdap;
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        sigma = el_pred(L, Eel);
    }
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];
    
    // tangent via the shared leading-mechanism helper (doc §7.4)
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    const std::vector<vec> dPhidsigma_l = { dPhidsigma };
    // tangent_none -> continuum: P_epsilon / invBhat feed the physical heat source r
    const int tangent_mode_eff = (tangent_mode == tangent_none)
        ? tangent_continuum : tangent_mode;
    const ContinuumTangent ct = rm.converged ? cpp_consistent_tangent(rm, L) : compute_tangent_operator(
        tangent_mode_eff, Bhat, kappa_j, dPhidsigma_l, Ds_j, L,
        [&]() -> std::vector<mat> {  // lazy: evaluated only in algorithmic mode
            // J2: dLambda/dsigma = deta_stress; only dSdE is corrected, thermal blocks stay continuum
            const std::vector<mat> dLambda_dsigma_l = { deta_stress(sigma) };
            return dLambda_dsigma_l;
        });
    dSdE = ct.Lt;
    const std::vector<vec>& P_epsilon = ct.P_epsilon;

    std::vector<double> P_theta(1);
    // P_theta = invBhat (dPhi/dtheta - dPhi/dsigma:L:alpha), 1/K; invBhat masked so elastic gives 0
    vec L_alpha;
    mul6(L, alpha, L_alpha);
    P_theta[0] = ct.invBhat(0, 0) * (dPhidtheta - sum(dPhidsigma%L_alpha));

    dSdT = -1.*L_alpha - (kappa_j[0]*P_theta[0]);
    
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
    
    double A_p = -Hp;
    double dA_pdp = -dHpdp;

    // Dissipation of the increment (the quantity Wm_d accumulates) and its derivatives through
    // dDp = P and dDEP = Lambdap P (flow direction held fixed over the increment)
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp;
    const double f_p = 0.5*(sum((sigma_start+sigma)%Lambdap) + A_p_start + A_p + dA_pdp*Dp);
    const vec dDgamma_dE = 0.5*dSdE.t()*DEP + f_p*P_epsilon[0];
    const double dDgamma_dT = 0.5*sum(dSdT%DEP) + f_p*P_theta[0];

    // heat of the increment and its derivatives (heat_source_increment)
    const double T1 = T + DT;
    heat_source_increment(Dgamma_loc, dDgamma_dE, dDgamma_dT, Deta, dSdE.t()*alpha, c_0/T1 + sum(dSdT%alpha),
                          T1, r, drdE, drdT);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;
    
    Wt += (T+DT)*Deta;
    Wt_r += (T+DT)*Deta_r;
    Wt_ir += (T+DT)*Deta_ir;
    
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
