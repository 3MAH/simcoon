/* This file is part of Simcoon private.
 
 Only part of Simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file damage_weibull.cpp
///@brief Isotropic damage with a weibull distribution of damage along a static loading

#include <iostream>
#include <fstream>
#include <math.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/lagrange.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/stats.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/damage.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>

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
    
void umat_damage_weibull(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt) {
                
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    UNUSED(DR);
    UNUSED(sigma_in);
    UNUSED(solver_type);
    
    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double alpha_iso = props(2);
    double alpha_weibull = props(3);
    double beta_weibull = props(4);
    
    //Set the Lagrange multipliers coefficients
    double c_lambda = props(5);
    double p0_lambda = props(6);
    double n_lambda = props(7);
    double alpha_lambda = props(8);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    ///@brief Temperature initialization
    double T_init = statev(0);
    ///@brief Damage variable
    double Damage = statev(1);
    
    //Elastic stiffness tensors
    L = L_iso(E, nu, "Enu");
    mat M = M_iso(E, nu, "Enu");    
    
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        Damage = precision_umat;
        sigma = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }

    vec sigma_eff = zeros(6);
    
    double lambda_D = lagrange_pow_1(Damage, c_lambda, p0_lambda, n_lambda, alpha_lambda);
    double dlambda_D = dlagrange_pow_1(Damage, c_lambda, p0_lambda, n_lambda, alpha_lambda);

    double A_Damage = 0.;
    double Rv = 0.;
    if (Mises_stress(sigma) > simcoon::iota) {
        Rv = (2./3.)*(1.+nu)+3.*(1.-2.*nu)*pow(tr(sigma)/(3.*Mises_stress(sigma)),2.);
        A_Damage = pow(Mises_stress(sigma),2.)*Rv/(2.*E*pow(1.-Damage,2.));
    }
    double A_Damage_start = A_Damage;
    
    double Weibull = 0.;
    double dWeibulldsigma_eq = 0.;
    
    //Variables values at the start of the increment
    vec sigma_start = sigma;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = Damage;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init);
    mat L_tilde = (1.-Damage)*L;
    sigma = el_pred(L_tilde, Eel, ndi);
    sigma_eff = el_pred(L, Eel, ndi);
    
    //Define the Damage function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = ones(1);

    vec Lambda_Damage = zeros(6);
    mat dStildedDamage = zeros(6,6);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = zeros(6);
    mat K = zeros(1,1);
    
    mat B_conc = L*inv(L_tilde);         //stress "localization factor" in damage
    
    double dPhidDamage=0.;
    vec dPhidsigma = zeros(6);
    vec dPhidsigma_eff = zeros(6);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    
    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        Damage = s_j(0);
        L_tilde = (1.-s_j(0))*L;
        
        Weibull = cumul_distrib_weibull(Mises_stress(sigma_eff), alpha_weibull, beta_weibull);
        dWeibulldsigma_eq = proba_distrib_weibull(Mises_stress(sigma_eff), alpha_weibull, beta_weibull);
        
        lambda_D = lagrange_pow_1(Damage, c_lambda, p0_lambda, n_lambda, alpha_lambda);
        dlambda_D = dlagrange_pow_1(Damage, c_lambda, p0_lambda, n_lambda, alpha_lambda);

        dPhidsigma_eff = dWeibulldsigma_eq*eta_stress(sigma_eff);
        dPhidDamage = -1. - dlambda_D;
        //dPhidDamage = -1.;
        
        //compute Phi and the derivatives
        Phi(0) = Weibull - Damage - lambda_D;
        //Phi(0) = Weibull - Damage;
        
        //Compute the explicit "damage direction" and flow direction
        dStildedDamage = 1./(pow(1.-Damage,2.))*M;
        Lambda_Damage = dStildedDamage*sigma;
        kappa_j[0] = L_tilde*Lambda_Damage;
        
        K(0,0) = dPhidDamage;
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        //B(0, 0) = K(0,0);
        Y_crit(0) = 1.;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init);
        sigma = el_pred(L_tilde, Eel, ndi);
        sigma_eff = el_pred(L, Eel, ndi);
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    double DDamage = Ds_j[0];
    
    dPhidsigma = dWeibulldsigma_eq*(B_conc.t()*eta_stress(sigma_eff));
    dStildedDamage = 1./(pow(1.-Damage,2.))*M;
    Lambda_Damage = dStildedDamage*sigma;
    kappa_j[0] = L_tilde*Lambda_Damage;
    
    //Computation of the tangent modulus
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);
    
    vec op = zeros(1);
    mat delta = eye(1,1);
    
    for (int i=0; i<1; i++) {
        if(Ds_j[i] > simcoon::iota)
            op(i) = 1.;
    }
    
    mat Bbar = zeros(1,1);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1-op(i)*op(j));
        }
    }
    
    mat invBbar = zeros(1, 1);
    mat invBhat = zeros(1, 1);
    invBbar = inv(Bbar);
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 1; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    std::vector<vec> P_epsilon(1);
    P_epsilon[0] = invBhat(0, 0)*(L_tilde*dPhidsigma);
    
    Lt = L_tilde - (kappa_j[0]*P_epsilon[0].t());
    //Lt = L_tilde;
    
    if (Mises_stress(sigma) > simcoon::iota) {
        Rv = (2./3.)*(1.+nu)+3.*(1.-2.*nu)*pow(tr(sigma)/(3.*Mises_stress(sigma)),2.);
        A_Damage = pow(Mises_stress(sigma),2.)*Rv/(2.*E*pow(1.-Damage,2.));
    }
    double Dgamma_loc = 0.5*(A_Damage_start + A_Damage)*DDamage;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot)) - 0.5*(A_Damage_start + A_Damage)*DDamage;
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    ///@brief statev evolving variables
    //statev
    statev(0) = T_init;
    statev(1) = Damage;    
}
    
} //namespace simcoon
