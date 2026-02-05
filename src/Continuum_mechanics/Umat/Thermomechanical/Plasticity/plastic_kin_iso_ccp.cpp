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

///@file plastic_kin_iso_ccp.cpp
///@brief User subroutine for elastic-plastic materials in 1D-2D-3D case
///@brief This subroutines uses a convex cutting plane algorithm
///@brief Linear Kinematical hardening coupled with a power-law hardenig is considered
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Plasticity/plastic_kin_iso_ccp.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
///@brief The thermomechanical elastic-plastic UMAT with kinematic + isotropic hardening requires 9 constants:

///@brief props[0] : density
///@brief props[1] : specific heat capacity
///@brief props[2] : Young modulus
///@brief props[3] : Poisson ratio
///@brief props[4] : CTE
///@brief props[5] : J2 equivalent yield stress limit : sigmaY
///@brief props[6] : hardening parameter k
///@brief props[7] : exponent m
///@brief props[8] : linear kinematical hardening h

///@brief The elastic-plastic UMAT with kinematic + isotropic hardening requires 14 statev:
///@brief statev[0] : T_init : Initial temperature
///@brief statev[1] : Accumulative plastic parameter: p
///@brief statev[2] : Plastic strain 11: EP(0,0)
///@brief statev[3] : Plastic strain 22: EP(1,1)
///@brief statev[4] : Plastic strain 33: EP(2,2)
///@brief statev[5] : Plastic strain 12: EP(0,1) (*2)
///@brief statev[6] : Plastic strain 13: EP(0,2) (*2)
///@brief statev[7] : Plastic strain 23: EP(1,2) (*2)
///@brief statev[8] : Backstress 11: X(0,0)
///@brief statev[9] : Backstress 11: X(1,1)
///@brief statev[10] : Backstress 11: X(2,2)
///@brief statev[11] : Backstress 11: X(0,1)
///@brief statev[12] : Backstress 11: X(0,2)
///@brief statev[13] : Backstress 11: X(1,2)


void umat_plasticity_kin_iso_CCP_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double E = props(2);
    double nu= props(3);
    double alpha_iso = props(4);
    double sigmaY = props(5);
    double k=props(6);
    double m=props(7);
    double kX = props(8);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Elastic stiffness tensor and thermal tensor
    mat L = L_iso(E, nu, "Enu");
    mat M = M_iso(E, nu, "Enu");
    
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
    vec a = zeros(6);
    a(0) = statev(8);
    a(1) = statev(9);
    a(2) = statev(10);
    a(3) = statev(11);
    a(4) = statev(12);
    a(5) = statev(13);
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    a = rotate_stress(a, DR);
    
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        sigma = vide;
        EP = vide;
        a = vide;
        p = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
        
        Wt = 0.;
        Wt_r = 0.;
        Wt_ir = 0.;        
    }
    
    //Additional parameters and variables
    double c_0 = rho*c_p;
    vec X = kX*(a%Ir05());
    
    double Hp=0.;
    double dHpdp=0.;
    
    if (p > simcoon::iota)	{
        dHpdp = m*k*pow(p, m-1);
        Hp = k*pow(p, m);
    }
    else {
        dHpdp = 0.;
        Hp = 0.;
    }
    
    //Variables values at the start of the increment
    vec sigma_start = sigma;
    vec EP_start = EP;
    vec a_start = a;
    vec X_start = X;
    
    double A_p_start = -Hp;
    vec A_a_start = -X_start;
    
    //Variables required for the loop
    vec s_j = zeros(1);
    s_j(0) = p;
    vec Ds_j = zeros(1);
    vec ds_j = zeros(1);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    sigma = el_pred(L, Eel, ndi);
    
    //Define the plastic function and the stress
    vec Phi = zeros(1);
    mat B = zeros(1,1);
    vec Y_crit = zeros(1);
    
    double dPhidp=0.;
    vec dPhida = zeros(6);
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec Lambdap = eta_stress(sigma-X);
    vec Lambdaa = eta_stress(sigma-X);
    std::vector<vec> kappa_j(1);
    kappa_j[0] = L*Lambdap;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    //Loop
    for (compteur = 0; ((compteur < simcoon::maxiter_umat) && (error > simcoon::precision_umat)); compteur++) {
        
        p = s_j(0);
        if (p > simcoon::iota)	{
            dHpdp = m*k*pow(p, m-1);
            Hp = k*pow(p, m);
        }
        else {
            dHpdp = 0.;
            Hp = 0.;
        }
        dPhidsigma = eta_stress(sigma-X);
        dPhidp = -1.*dHpdp;
        dPhida = -1.*kX*(eta_stress(sigma - X)%Ir05());
        
        //compute Phi and the derivatives
        Phi(0) = Mises_stress(sigma-X) - Hp - sigmaY;
        
        Lambdap = eta_stress(sigma-X);
        Lambdaa = eta_stress(sigma-X);
        kappa_j[0] = L*Lambdap;
        
        K(0,0) = dPhidp + sum(dPhida%Lambdaa);
        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);
        EP = EP + ds_j(0)*Lambdap;
        a = a + ds_j(0)*Lambdaa;
        X = kX*(a%Ir05());
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        sigma = el_pred(L, Eel, ndi);
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];
    vec Da = a - a_start;
    
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
    P_epsilon[0] = invBhat(0, 0)*(L*dPhidsigma);
    std::vector<double> P_theta(1);
    P_theta[0] = dPhidtheta - sum(dPhidsigma%(L*alpha));
    
    dSdE = L - (kappa_j[0]*P_epsilon[0].t());
    dSdT = -1.*L*alpha - (kappa_j[0]*P_theta[0]);
    
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
    
    vec Gamma_epsilon = zeros(6);
    double Gamma_theta = 0.;
    
    vec N_epsilon = zeros(6);
    double N_theta = 0.;
    
    double A_p = -Hp;
    double dA_pdp = -dHpdp;
    vec A_a = -X;
    mat dA_ada = -1.*kX*Ireal();
    
    //    double A_theta = 0;
    
    if(DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        Gamma_epsilon = (dSdE*DEP)*(1./DTime) + (dA_pdp)*P_epsilon[0]*(Dp/DTime) + sum((dA_ada*Lambdaa)%P_epsilon[0])*(Da/DTime) + A_p/DTime*P_epsilon[0] + sum(A_a%Lambdaa)*P_epsilon[0]*(1./DTime) + sum(sigma%Lambdap)*P_epsilon[0]/DTime;
        Gamma_theta = dA_pdp*P_theta[0]*(Dp/DTime) + sum((dA_ada*Lambdaa)%Da)*P_theta[0]/DTime + A_p/DTime*P_theta[0] + sum(A_a%Lambdaa)*P_theta[0]*(1./DTime) + sum(dSdT%DEP)*(1./DTime) + sum(sigma%Lambdap)*P_theta[0]/DTime;
        
        N_epsilon = -1./DTime*(T + DT)*(dSdE*alpha);
        N_theta = -1./DTime*(T + DT)*sum(dSdT%alpha) -1.*Deta/DTime - rho*c_p*(1./DTime);
        
        drdE = N_epsilon + Gamma_epsilon;
        drdT = N_theta + Gamma_theta;
        
        r = sum(N_epsilon%DEtot) + N_theta*DT + sum(Gamma_epsilon%DEtot) + Gamma_theta*DT;
    }
    
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp + 0.5*sum((A_a_start + A_a)%Da);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP)) - 0.5*sum((A_a_start + A_a)%Da);
    Wm_ir += -0.5*(A_p_start + A_p)*Dp;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir = (T+0.5*DT)*Deta_ir;
    
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
    
    statev(8) = a(0);
    statev(9) = a(1);
    statev(10) = a(2);
    statev(11) = a(3);
    statev(12) = a(4);
    statev(13) = a(5);
}
    
} //namespace simcoon