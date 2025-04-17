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

///@file damage_LLD_m.cpp
///@brief Modified Ladeveze-Le Dantec model that accounts for uniaxial damage and anisotopic plasticity
///@brief Damage evolution is considered as a function of time

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Simulation/Phase/output.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Damage/damage_LLD_0.hpp>

using namespace std;
using namespace arma;

namespace simcoon {
    
///@brief Material properties
///@param props(0) : E Young's modulus
///@param props(1) : nu Poisson's ration
///@param props(2) : alpha CTE
///@param props(3) : alphaD Damage evolution parameter alpha
///@param props(4) : betaD Damage evolution parameter beta
///@param props(5) : lambdaD Damage evolution parameter lambda
///@param props(6) : deltaD Damage evolution parameter delta

void umat_damage_LLD_0(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt) {
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double Tinit = statev(0);
    double d_22 = statev(1);
    double d_12 = statev(2);
    double p_ts = statev(3);    //Accumulated plasticity
        
    vec EP(6);
    EP(0) = statev(4);
    EP(1) = statev(5);
    EP(2) = statev(6);
    EP(3) = statev(7);
    EP(4) = statev(8);
    EP(5) = statev(9);

    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    
    double axis = props(0);
    double EL = props(1);
    double ET = props(2);
    double nuTL = props(3);
    double nuTT = props(4);
    double GLT = props(5);
    double alphaL = props(6);
    double alphaT = props(7);
    UNUSED(axis);  //Hide a warning
    
    //Shear damage
    double Y_12_0 = props(8);
    double Y_12_c = props(9);
    
    //Transverse damage
    double Y_22_0 = props(10);
    double Y_22_c = props(11);
    double Y_22_u = props(12);
    double b = props(13);
    
    //Coupled transverse-shear yield & plasticity
    double A_ts = props(14);
    double sigma_ts_0 = props(15);
    double alpha_ts = props(16);
    double beta_ts = props(17);
    
    //Set the Lagrange multipliers coefficients
    double c_lambda = 1.E-2;
    double p0_lambda = 1.E-1;
    double n_lambda = 1.0;
    double alpha_lambda = 0.;

    //Elastic stiffness tensor
    L = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    
    ///@brief Initialization
    if (start) {

        Tinit = T;
        d_22 = 0.;
        d_12 = 0.;
        
        p_ts = 10.*simcoon::limit;
        EP = zeros(6);
        sigma = zeros(6);
    }
    
    vec sigma_start = sigma;
    vec EP_start = EP;
    
    double E1 = EL;
    double E2_0 = ET;
    double E3_0 = ET;
    double E2 = ET*(1-d_22);
    double E3 = ET*(1-d_22);
    double nu12 = nuTL*(E1/E2_0);
    double nu13 = nuTL*(E1/E3_0);
    double nu32 = nuTT;
    double nu23 = nu32;
    double G12_0 = GLT;
    double G12 = G12_0*(1-d_12);
    double G13 = G12_0*(1-d_12);
    double G23 = E3/(2.*(1.+nu32));
    
    // ######################  Elastic stiffness #################################
    //defines L
    L = L_isotrans(EL, ET, nuTL, nuTT, GLT, 1);
    
    // cout << "EL = " << EL << "/ ET = " << ET << "/ nuTL = " << nuTL << "/ nuTT = " << nuTT << "/ GLT = " << GLT << endl;
    // 
    // cout << "L Iso T = " << L << endl; 
    
    //definition of the CTE tensor
    vec alpha = zeros(6);
    alpha = alphaT*Ith();
    alpha(0) += alphaL-alphaT;                              //WARNING CTE tensor is defined according to L is the direction 1
    
    //Compute the elastic strain and the related stress
    vec Eel_start = Etot - alpha*(T-Tinit) - EP;
    vec DEel_start = DEtot - alpha*DT;
    vec Eel = Eel_start + DEel_start;
    
    vec sigma_eff = zeros(6);
    vec sigma_eff_ts = zeros(6);
    
    ///Compute the Yd_12 and Yd_22 which are the criterium associated with damage
    double Yd_12 = 0.;
    double Yd_13 = 0.;
    double Yd_22 = 0.;
    double Yd_33 = 0.;
    
    //Compute Y_l, Y_t and Y_ts, which are used in the evolution equation of damage
    double Y_ts = 0.;              //transvers/shear coupling
    double Y_t = 0.;               //transverse

    double lambda_22 = lagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double lambda_12 = lagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double dlambda_22 = 0.;
    double dlambda_12 = 0.;
    
    //Compute the Plasticity functions
    //Compute the explicit flow direction
    vec Lambdap_ts = zeros(6);
    
    //Define the plastic function and the stress
    double Hp_ts = 0.;
    
    double Phi_p_ts = 0.;

    double dPhip_tsdp_ts = 0.;
    vec dPhi_p_tsd_sigma = zeros(6);
    
    //Note : The sup function are not required here, since we utilize Kuhn-Tucker conditions instead (dD >= 0)
    vec Phi_d = zeros(2);
    vec Phi_p = zeros(1);
    vec Dd = zeros(2);
    vec dd = zeros(2);
    vec Dp = zeros(1);
    vec dp = zeros(1);
    
    vec Y_dcrit = ones(2);
    vec Y_pcrit = ones(1);
    
    mat denom_d = zeros(2,2);
    mat denom_p = zeros(1,1);
    
    double dYd_22dd = 0.;
    double dYd_33dd = 0.;
    double dYd_12dd = 0.;
    double dYd_13dd = 0.;
    UNUSED(dYd_33dd);  //Hide a warning
    
    double dY_tsdd_22 = 0.;
    double dY_tsdd_12 = 0.;

    int compteur = 0;
    double error = 1.;
    
    mat Theta_ts = zeros(6,6);
    Theta_ts(1,1) = A_ts;
    Theta_ts(2,2) = A_ts;
    Theta_ts(3,3) = 1.;
    Theta_ts(4,4) = 1.;
    
    error = 1.;
    if (Y_t > Y_22_u) {
        d_22 = 0.99;
        d_12 = 0.99;
        error = 0.;
    }
    
    //First we find the plasticity
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        //Plasticity computations
                //Compute the hardening
        if (p_ts > simcoon::limit)
            Hp_ts = beta_ts*pow(p_ts, alpha_ts);
        else
            Hp_ts = simcoon::iota;
        
        //effective stress
        sigma_eff = el_pred(L,Eel,ndi);
        sigma_eff_ts = Theta_ts*sigma_eff;
        
        //Determine the Phi functions and their derivatives
        Phi_p_ts = Mises_stress(sigma_eff_ts) - Hp_ts - sigma_ts_0;

        //Compute the explicit flow direction
        Lambdap_ts = eta_stress(sigma_eff_ts);
        
        if (p_ts > simcoon::limit)
            dPhip_tsdp_ts = -1.*alpha_ts*beta_ts*pow(p_ts, alpha_ts-1.);
        else
            dPhip_tsdp_ts = simcoon::iota;
        
        dPhi_p_tsd_sigma = Theta_ts*eta_stress(sigma_eff_ts); //Here as well
        
        //compute Phi and the derivatives
        Phi_p(0) = Phi_p_ts;
        
        denom_p(0, 0) = -1.*sum(dPhi_p_tsd_sigma % (L*Lambdap_ts)) + dPhip_tsdp_ts;
        
        Y_pcrit(0) = sigma_ts_0;
        
        Fischer_Burmeister_m(Phi_p, Y_pcrit, denom_p, Dp, dp, error);
        
        p_ts += dp(0);
        
        EP = EP + dp(0)*Lambdap_ts;
        Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
    }
    
    if(compteur == maxiter_umat)
        tnew_dt = 0.2;
    
    error = 1.;
    if (Y_t > Y_22_u) {
        d_22 = 0.99;
        d_12 = 0.99;
        error = 0.;
    }
    
    //Derivatives specific to damage:
    double dPhi_d_22d_Yts = 1./Y_22_c;
    double dPhi_d_12d_Yts = 1./Y_12_c;
    
    double dYts_d_Y22 = 0.;
    double dYts_d_Y33 = 0.;
    double dYts_d_Y12 = 0.;
    double dYts_d_Y13 = 0.;
    
    double dY22_d_s22 = 0.;
    double dY33_d_s33 = 0.;
    double dY12_d_s12 = 0.;
    double dY13_d_s13 = 0.;
    
    vec dPhi_d_22d_sigma = zeros(6);
    vec dPhi_d_12d_sigma = zeros(6);
    
    double dPhi_d_12d_12 = 0.;
    double dPhi_d_12d_22 = 0.;
    double dPhi_d_22d_12 = 0.;
    double dPhi_d_22d_22 = 0.;
    
    mat dStildedd22 = zeros(6,6);
    mat dStildedd12 = zeros(6,6);
    
    vec Lambdad_22 = zeros(6);
    vec Lambdad_12 = zeros(6);
    
    
    //So it is forced to enter the damage loop once
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        mat L_tilde = L_ortho(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23, "EnuG");
        
        // cout << "d22 = " << d_22 << "/ d12 = " << d_12 << endl;
        // 
		// cout << "E1 = " << E1 << "/ E2 = " << E2 << "/ E3 = " << E3 << "/ nu12 = " << nu12 << "/ nu13 = " << nu13 << "/ nu23 = " << nu23 << "/ G12 = " << G12 << "/ G13 = " << G13 << "/ G23 = " << G23 << endl;
    // 
		// cout << "L ortho = " << L_tilde << endl;
        
        //damaged modulus
        //Compute the elastic strain and the related stress
        Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
        sigma = el_pred(L_tilde, Eel, ndi);
        
        if (fabs(sigma(1)) < simcoon::iota )
            Yd_22 = 0.;
        else
            Yd_22 = 0.5*(pow(Macaulay_p(sigma(1)),2.)/(E2_0*pow(1.-d_22,2.)));

        if (fabs(sigma(2)) < simcoon::iota )
            Yd_33 = 0.;
        else
            Yd_33 = 0.5*(pow(Macaulay_p(sigma(2)),2.)/(E2_0*pow(1.-d_22,2.)));
        
        if (fabs(sigma(3)) < simcoon::iota )
            Yd_12 = 0.;
        else
            Yd_12 = 0.5*(pow(sigma(3),2.)/(G12_0*pow(1.-d_12,2.)));

        if (fabs(sigma(4)) < simcoon::iota )
            Yd_13 = 0.;
        else
            Yd_13 = 0.5*(pow(sigma(4),2.)/(G12_0*pow(1.-d_12,2.)));
                
        //Compute Y_t and Y_ts, which are used in the evolution equation of damage
        Y_ts = sqrt(Yd_12 + Yd_13 + b*(Yd_22 + Yd_33));     //transvers/shear coupling
        Y_t = sqrt(Yd_22 + Yd_33);               //transverse
        
        lambda_12 = lagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        lambda_22 = lagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        //Preliminaries to compute damage the derivatives of Phi
        if (fabs(sigma(1)) < simcoon::iota )
            dYd_22dd = 0.;
        else
            dYd_22dd = -0.25*(pow(Macaulay_p(sigma(1)),2.)/(E2_0*pow(1-d_22,3.)));

        if (fabs(sigma(2)) < simcoon::iota )
            dYd_33dd = 0.;
        else
            dYd_33dd = -0.25*(pow(Macaulay_p(sigma(2)),2.)/(E2_0*pow(1-d_22,3.)));
        
        if (fabs(sigma(3)) < simcoon::iota )
            dYd_12dd = 0.;
        else
            dYd_12dd = -0.25*(pow(sigma(3),2.)/(G12_0*pow(1-d_12,3.)));

        if (fabs(sigma(4)) < simcoon::iota )
            dYd_13dd = 0.;
        else
            dYd_13dd = -0.25*(pow(sigma(4),2.)/(G12_0*pow(1-d_12,3.)));

        
        
        if (fabs(Yd_12 + Yd_13 + b*(Yd_22 + Yd_33)) < simcoon::iota ) {
            dY_tsdd_22 = 0.;
            dY_tsdd_12 = 0.;
        }
        else {
            dY_tsdd_22 = 0.5*(b*(dYd_22dd+dYd_22dd))*pow(Yd_12 + Yd_13+ b*(Yd_22+Yd_33),-0.5);
            dY_tsdd_12 = 0.5*(dYd_12dd + dYd_13dd)*pow(Yd_12 + Yd_13+ b*(Yd_22+Yd_33),-0.5);
        }
                    
        dlambda_22 = dlagrange_pow_1(d_22, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda_12 = dlagrange_pow_1(d_12, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        
        //Compute the derivatives of sigma
        dPhi_d_22d_Yts = 1./Y_22_c;
        dPhi_d_12d_Yts = 1./Y_12_c;
        
        dYts_d_Y22 = 0.;
        dYts_d_Y33 = 0.;
        dYts_d_Y12 = 0.;
        dYts_d_Y13 = 0.;
        
        if (Y_ts > simcoon::limit) {
            dYts_d_Y22 = b/(2.*Y_ts);
            dYts_d_Y33 = b/(2.*Y_ts);
            dYts_d_Y12 = 1./(2.*Y_ts);
            dYts_d_Y13 = 1./(2.*Y_ts);
        }
        
        dY22_d_s22 = Macaulay_p(sigma(1))/(E2_0*(1.-d_22));
        dY33_d_s33 = Macaulay_p(sigma(2))/(E2_0*(1.-d_22));
        dY12_d_s12 = Macaulay_p(sigma(3))/(G12_0*(1.-d_12));
        dY13_d_s13 = Macaulay_p(sigma(4))/(G12_0*(1.-d_12));
        
        dPhi_d_22d_sigma(0) = 0.;
        dPhi_d_22d_sigma(1) = dPhi_d_22d_Yts*dYts_d_Y22*dY22_d_s22;
        dPhi_d_22d_sigma(2) = dPhi_d_22d_Yts*dYts_d_Y33*dY33_d_s33;
        dPhi_d_22d_sigma(3) = dPhi_d_22d_Yts*dYts_d_Y12*dY12_d_s12;
        dPhi_d_22d_sigma(4) = dPhi_d_22d_Yts*dYts_d_Y13*dY13_d_s13;
        dPhi_d_22d_sigma(5) = 0.;
        
        dPhi_d_12d_sigma(0) = 0.;
        dPhi_d_12d_sigma(1) = dPhi_d_12d_Yts*dYts_d_Y22*dY22_d_s22;
        dPhi_d_12d_sigma(2) = dPhi_d_12d_Yts*dYts_d_Y33*dY33_d_s33;
        dPhi_d_12d_sigma(3) = dPhi_d_12d_Yts*dYts_d_Y12*dY12_d_s12;
        dPhi_d_12d_sigma(4) = dPhi_d_12d_Yts*dYts_d_Y13*dY13_d_s13;
        dPhi_d_12d_sigma(5) = 0.;
        
        //Compute the explicit "damage direction" and flow direction
        dStildedd22 = { {0.,0.,0.,0.,0.,0.},
            {0,E2_0/pow(E2,2.),-nu23*E2_0/pow(E2,2.),0,0,0},
            {0,-nu23*E2_0/pow(E2,2.),E2_0/pow(E2,2.),0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0} };
        
        dStildedd12 = { {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,0,0,0},
            {0,0,0,G12_0/pow(G12,2.),0,0},
            {0,0,0,0,G12_0/pow(G12,2.),0},
            {0,0,0,0,0,0} };

        Lambdad_22 = dStildedd22*sigma;
        Lambdad_12 = dStildedd12*sigma;
        
        //compute Phi and the derivatives
        Phi_d(0) = Macaulay_p(Y_ts - Y_22_0)/Y_22_c - lambda_22 - d_22;
        Phi_d(1) = Macaulay_p(Y_ts - Y_12_0)/Y_12_c - lambda_12 - d_12;
        
        denom_d(0, 0) = -1.*sum(dPhi_d_22d_sigma%Lambdad_22) + Macaulay_p(dY_tsdd_22)/Y_22_c - dlambda_22 - 1.;
        denom_d(0, 1) = -1.*sum(dPhi_d_22d_sigma%Lambdad_12) + Macaulay_p(dY_tsdd_12)/Y_22_c;
        
        denom_d(1, 0) = -1.*sum(dPhi_d_12d_sigma%Lambdad_22) + Macaulay_p(dY_tsdd_22)/Y_12_c;
        denom_d(1, 1) = -1.*sum(dPhi_d_12d_sigma%Lambdad_12) + Macaulay_p(dY_tsdd_12)/Y_12_c - dlambda_12 - 1.;
        
        
        Fischer_Burmeister_m(Phi_d, Y_dcrit, denom_d, Dd, dd, error);
        
        d_22 += dd(0);
        d_12 += dd(1);
        
        //Transverse Damage term
        if (Y_t > Y_22_u) {
            d_22 = 0.99;
            d_12 = 0.99;
            error = 0.;
        }
        
        //Update constitutive parameters
        E2 = ET*(1.-d_22);
        E3 = ET*(1.-d_22);
        G12 = G12_0*(1.-d_12);
        G13 = G12_0*(1.-d_12);
    }
    
    if(compteur == maxiter_umat)
        tnew_dt = 0.2;
    
    //Update constitutive parameters
    E2 = ET*(1.-d_22);
    E3 = ET*(1.-d_22);
    G12 = G12_0*(1.-d_12);
    G13 = G12_0*(1.-d_12);
    
    mat L_tilde = L_ortho(E1,E2,E3,nu12,nu13,nu23,G12,G13,G23, "EnuG");
    
    //damaged modulus
    //Compute the elastic strain and the related stress
    Eel = Etot + DEtot - alpha*(T+DT-Tinit) - EP;
    sigma = el_pred(L_tilde, Eel, ndi);
    
    if((solver_type == 0)||(solver_type == 2)) {
    
		//Tangent modulus
		mat B = L*inv(L_tilde);         //stress "localization factor" in damage
		
		//Compute the derivatives
		dPhi_d_22d_Yts = 1./Y_22_c;
		dPhi_d_12d_Yts = 1./Y_12_c;

		dYts_d_Y22 = 0.;
		dYts_d_Y33 = 0.;
		dYts_d_Y12 = 0.;
		dYts_d_Y13 = 0.;
		
		if (Y_ts > simcoon::limit) {
			dYts_d_Y22 = b/(2.*Y_ts);
			dYts_d_Y33 = b/(2.*Y_ts);
			dYts_d_Y12 = 1./(2.*Y_ts);
			dYts_d_Y13 = 1./(2.*Y_ts);
		}
		
		dY22_d_s22 = Macaulay_p(sigma(1))/(E2_0*(1.-d_22));
		dY33_d_s33 = Macaulay_p(sigma(2))/(E2_0*(1.-d_22));
		dY12_d_s12 = Macaulay_p(sigma(3))/(G12_0*(1.-d_12));
		dY13_d_s13 = Macaulay_p(sigma(4))/(G12_0*(1.-d_12));
		
		dPhi_d_22d_sigma(0) = 0.;
		dPhi_d_22d_sigma(1) = dPhi_d_22d_Yts*dYts_d_Y22*dY22_d_s22;
		dPhi_d_22d_sigma(2) = dPhi_d_22d_Yts*dYts_d_Y33*dY33_d_s33;
		dPhi_d_22d_sigma(3) = dPhi_d_22d_Yts*dYts_d_Y12*dY12_d_s12;
		dPhi_d_22d_sigma(4) = dPhi_d_22d_Yts*dYts_d_Y13*dY13_d_s13;
		dPhi_d_22d_sigma(5) = 0.;
		
		dPhi_d_12d_sigma(0) = 0.;
		dPhi_d_12d_sigma(1) = dPhi_d_12d_Yts*dYts_d_Y22*dY22_d_s22;
		dPhi_d_12d_sigma(2) = dPhi_d_12d_Yts*dYts_d_Y33*dY33_d_s33;
		dPhi_d_12d_sigma(3) = dPhi_d_12d_Yts*dYts_d_Y12*dY12_d_s12;
		dPhi_d_12d_sigma(4) = dPhi_d_12d_Yts*dYts_d_Y13*dY13_d_s13;
		dPhi_d_12d_sigma(5) = 0.;
		
		dPhi_d_22d_22 = Macaulay_p(dY_tsdd_22)/Y_22_c - dlambda_22 - 1.;
		dPhi_d_22d_12 = Macaulay_p(dY_tsdd_12)/Y_22_c;
		
		dPhi_d_12d_22 = Macaulay_p(dY_tsdd_22)/Y_12_c;
		dPhi_d_12d_12 = Macaulay_p(dY_tsdd_12)/Y_12_c - dlambda_12 - 1.;
    
		// dPhi_p_tsd_sigma = (B*Theta_ts*eta_stress(sigma_eff_ts));
		dPhi_p_tsd_sigma = (Theta_ts*eta_stress(sigma_eff_ts));
		
		//Compute the explicit "damage direction" and flow direction
		dStildedd22 = { {0.,0.,0.,0.,0.,0.},
							{0,E2_0/pow(E2,2.),-nu23*E2_0/pow(E2,2.),0,0,0},
							{0,-nu23*E2_0/pow(E2,2.),E2_0/pow(E2,2.),0,0,0},
							{0,0,0,0,0,0},
							{0,0,0,0,0,0},
							{0,0,0,0,0,0} };
		
		dStildedd12 = { {0,0,0,0,0,0},
							{0,0,0,0,0,0},
							{0,0,0,0,0,0},
							{0,0,0,G12_0/pow(G12,2.),0,0},
							{0,0,0,0,G12_0/pow(G12,2.),0},
							{0,0,0,0,0,0} };

		Lambdad_22 = dStildedd22*sigma;
		Lambdad_12 = dStildedd12*sigma;    
		Lambdap_ts = eta_stress(sigma_eff_ts);

		std::vector<vec> kappa_j(3);
		std::vector<vec> kappa_tilde_j(3);
		kappa_j[0] = L_tilde*Lambdad_22;
		kappa_j[1] = L_tilde*Lambdad_12;
		kappa_j[2] = L_tilde*Lambdap_ts;
		
		kappa_tilde_j[0] = L*Lambdad_22;
		kappa_tilde_j[1] = L*Lambdad_12;
		kappa_tilde_j[2] = L*Lambdap_ts;    
		
		mat K = zeros(3,3);
		K(0,0) = dPhi_d_22d_22;
		K(0,1) = dPhi_d_22d_12;
		K(0,2) = 0.;
		
		K(1,0) = dPhi_d_12d_22;
		K(1,1) = dPhi_d_12d_12;
		K(1,2) = 0.;
		
		K(2,0) = 0.;
		K(2,1) = 0.;
		K(2,2) = dPhip_tsdp_ts;

		
		mat Bhat = zeros(3, 3);
		Bhat(0, 0) = sum(dPhi_d_22d_sigma % kappa_j[0]) - K(0,0);
		Bhat(0, 1) = sum(dPhi_d_22d_sigma % kappa_j[1]) - K(0,1);
		Bhat(0, 2) = sum(dPhi_d_22d_sigma % kappa_j[2]) - K(0,2);

		Bhat(1, 0) = sum(dPhi_d_12d_sigma % kappa_j[0]) - K(1,0);
		Bhat(1, 1) = sum(dPhi_d_12d_sigma % kappa_j[1]) - K(1,1);
		Bhat(1, 2) = sum(dPhi_d_12d_sigma % kappa_j[2]) - K(1,2);
		
		Bhat(2, 0) = sum(dPhi_p_tsd_sigma % kappa_tilde_j[0]) - K(2,0);
		Bhat(2, 1) = sum(dPhi_p_tsd_sigma % kappa_tilde_j[1]) - K(2,1);
		Bhat(2, 2) = sum(dPhi_p_tsd_sigma % kappa_tilde_j[2]) - K(2,2);
    
		vec op = zeros(3);
		mat delta = eye(3,3);

		if(Dd(0) > simcoon::iota)
			op(0) = 1.;
		if(Dd(1) > simcoon::iota)
			op(1) = 1.;
		if(Dp(0) > simcoon::iota)
			op(2) = 1.;
		
		
		mat Bbar = zeros(3,3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
			}
		}
		
		mat invBbar = zeros(3, 3);
		mat invBhat = zeros(3, 3);
		invBbar = inv(Bbar);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
			}
		} 
		
		std::vector<vec> P_epsilon(3);
		P_epsilon[0] = invBhat(0, 0)*dPhi_d_22d_sigma + invBhat(1, 0)*dPhi_d_12d_sigma + invBhat(2, 0)*dPhi_p_tsd_sigma;
		P_epsilon[1] = invBhat(0, 1)*dPhi_d_22d_sigma + invBhat(1, 1)*dPhi_d_12d_sigma + invBhat(2, 1)*dPhi_p_tsd_sigma;
		P_epsilon[2] = invBhat(0, 2)*dPhi_d_22d_sigma + invBhat(1, 2)*dPhi_d_12d_sigma + invBhat(2, 2)*dPhi_p_tsd_sigma;
		
		Lt = L_tilde - (kappa_j[0]*P_epsilon[0].t() + kappa_j[1]*P_epsilon[1].t() + kappa_j[2]*P_epsilon[2].t());
    }
    else if(solver_type == 1) {
        sigma_in = -L*(Etot - Eel);
    }
    
    statev(0) = Tinit;
    statev(1) = d_22;
    statev(2) = d_12;
    
    statev(3) = p_ts;
               
    statev(4) = EP(0);
    statev(5) = EP(1);
    statev(6) = EP(2);
    statev(7) = EP(3);
    statev(8) = EP(4);
    statev(9) = EP(5);
            
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    
    //For analysis
    statev(10) = Y_t;
    statev(11) = Y_ts;
    statev(12) = Hp_ts;
	
	
    // vector<mat> moduli_outputed;
    // moduli_outputed = {Lt, L, L_tilde, B};
    // moduli_outputed = {L};
    // output_moduli(moduli_outputed,13, statev);
    
    
    // cout << "************************************" << endl;
}


    
} //namespace simcoon
