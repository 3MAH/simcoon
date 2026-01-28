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
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Plasticity/Generic_chaboche_ccp.hpp>

using namespace std;
using namespace arma;

namespace simcoon {
    
///@brief The mechanical elastic-plastic UMAT with kinematic + isotropic hardening requires 7 constants:

///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE
///@brief props[3] : J2 equivalent yield stress limit : sigmaY
///@brief props[4] : hardening parameter k
///@brief props[5] : exponent m
///@brief props[6] : linear kinematical hardening h

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


void umat_generic_chaboche_CCP(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    
    double E = props(0);
    double nu = props(1);
    double G = props(2);    
    double alpha_iso = props(3);
    double sigmaY = props(4);

    int N_iso_hard = int(props(5));
    int N_kin_hard = int(props(6));
    int criteria = int(props(7));    

    vec Q = zeros(N_iso_hard);
    vec b = zeros(N_iso_hard);
    for (int i=0; i<N_iso_hard; i++) {
        Q(i)=props(8+i*2);
        b(i)=props(8+i*2+1);
    }    

    vec C = zeros(N_kin_hard);
    vec D = zeros(N_kin_hard);    
    for (int i=0; i<N_kin_hard; i++) {
        C(i)=props(8+N_iso_hard*2+i*2);
        D(i)=props(8+N_iso_hard*2+i*2+1);
    }    

    vec criteria_params;

    switch (criteria)
    {
        case 0: {//Mises surface
            break;       
        } 
        case 1: {//Hill
            double F_hill = props(8+N_iso_hard*2+N_kin_hard*2);
            double G_hill = props(9+N_iso_hard*2+N_kin_hard*2);
            double H_hill = props(10+N_iso_hard*2+N_kin_hard*2);
            double L_hill = props(11+N_iso_hard*2+N_kin_hard*2);
            double M_hill = props(12+N_iso_hard*2+N_kin_hard*2);
            double N_hill = props(13+N_iso_hard*2+N_kin_hard*2);
            criteria_params = {F_hill,G_hill,H_hill,L_hill,M_hill,N_hill};
            break;
        }
        case 2: {//DFA
            double F_dfa = props(8+N_iso_hard*2+N_kin_hard*2);
            double G_dfa = props(9+N_iso_hard*2+N_kin_hard*2);
            double H_dfa = props(10+N_iso_hard*2+N_kin_hard*2);
            double L_dfa = props(11+N_iso_hard*2+N_kin_hard*2);
            double M_dfa = props(12+N_iso_hard*2+N_kin_hard*2);
            double N_dfa = props(13+N_iso_hard*2+N_kin_hard*2);
            double K_dfa = props(14+N_iso_hard*2+N_kin_hard*2);
            criteria_params = {F_dfa,G_dfa,H_dfa,L_dfa,M_dfa,N_dfa,K_dfa}; 
            break;   
        }
        case 3: {//ANI
            double P11 = props(8+N_iso_hard*2+N_kin_hard*2);
            double P22 = props(9+N_iso_hard*2+N_kin_hard*2);
            double P33 = props(10+N_iso_hard*2+N_kin_hard*2);
            double P12 = props(11+N_iso_hard*2+N_kin_hard*2);
            double P13 = props(12+N_iso_hard*2+N_kin_hard*2);
            double P23 = props(13+N_iso_hard*2+N_kin_hard*2);
            double P44 = props(14+N_iso_hard*2+N_kin_hard*2);
            double P55 = props(15+N_iso_hard*2+N_kin_hard*2);
            double P66 = props(16+N_iso_hard*2+N_kin_hard*2);        
            criteria_params = {P11,P22,P33,P12,P13,P23,P44,P55,P66};
            break;   
        }             
        default: {
            break;
        }
    }

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
    
    double Hp = statev(8);

    ///@brief a is the internal variable associated with kinematical hardening
    std::vector<vec> a_N(N_kin_hard);
    std::vector<vec> X_N(N_kin_hard);

    for (int i=0; i<N_kin_hard; i++) {
        a_N[i] = zeros(6);
        X_N[i] = zeros(6);

        a_N[i](0) = statev(8+i*12+1);
        a_N[i](1) = statev(8+i*12+2);
        a_N[i](2) = statev(8+i*12+3);
        a_N[i](3) = statev(8+i*12+4);
        a_N[i](4) = statev(8+i*12+5);
        a_N[i](5) = statev(8+i*12+6);

        X_N[i](0) = statev(8+i*12+7);
        X_N[i](1) = statev(8+i*12+8);
        X_N[i](2) = statev(8+i*12+9);
        X_N[i](3) = statev(8+i*12+10);                
        X_N[i](4) = statev(8+i*12+11);
        X_N[i](5) = statev(8+i*12+12);
    }
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);

    vec X = zeros(6);
    for (int i=0; i<N_kin_hard; i++) {
        a_N[i] = rotate_strain(a_N[i], DR);
        X += X_N[i];
    }
    
    //Elstic stiffness tensor
    L = L_cubic(E, nu, G, "EnuG");
        
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        vec vide = zeros(6);
        sigma = vide;
        EP = vide;

        for (int i=0; i<N_kin_hard; i++) {
            a_N[i] = vide;
        }

        p = 0.;
        Hp = 0.;
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
    //Additional parameters and variables
    double dHpdp=0.;
    
    if (p > sim_iota)	{
        for (int i=0; i<N_iso_hard; i++) {        
            dHpdp += b[i]*(Q[i]-Hp);
        }
    }
    else {
        dHpdp = 0.;
    }
    
    //Variables values at the start of the increment
    vec sigma_start = sigma;
    vec EP_start = EP;

    std::vector<vec> a_N_start(N_kin_hard);
    std::vector<vec> X_N_start(N_kin_hard);
    for (int i=0; i<N_kin_hard; i++) {
        a_N_start[i] = a_N[i];
        X_N_start[i] = X_N[i];        
    }
    
    double A_p_start = -Hp;

    std::vector<vec> A_a_N_start(N_kin_hard);    
    for (int i=0; i<N_kin_hard; i++) {
        A_a_N_start[i] = -X_N_start[i];
    }
    
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
    std::vector<vec> dPhida_N(N_kin_hard);  
    for (int i=0; i<N_kin_hard; i++) {
        dPhida_N[i] = zeros(6);
    }    
    vec dPhidsigma = zeros(6);
    double dPhidtheta = 0.;
    
    //Compute the explicit flow direction
    vec Lambdap = zeros(6);
    std::vector<vec> Lambdaa_N(N_kin_hard);      
    switch (criteria)
    {
        case 0: {//Mises surface
            Lambdap = eta_stress(sigma-X);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = eta_stress(sigma-X) - D[i]*a_N[i];
            }            
            break;   
        }     
        case 1: {//Hill
            Lambdap = dHill_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dHill_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }        
            break;
        }
        case 2: {//DFA
            Lambdap = dDFA_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dDFA_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }        
            break;
        }
        case 3: {//ANI
            Lambdap = dAni_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dAni_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }        
            break;
        }
        default: {
            break;
        }
    }

    std::vector<vec> kappa_j(1);
    kappa_j[0] = L*Lambdap;
    mat K = zeros(1,1);
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;
    
    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        p = s_j(0);
        if (p > sim_iota)	{
            dHpdp = 0.;
            for (int i=0; i<N_iso_hard; i++) {        
                dHpdp += b[i]*(Q[i]-Hp);
            }            
        }
        else {
            dHpdp = 0.;
        }
        dPhidp = -1.*dHpdp;

    switch (criteria)
    {
        case 0: {//Mises surface
            dPhidsigma = eta_stress(sigma-X);
            for (int i=0; i<N_kin_hard; i++) {
                dPhida_N[i] = -1.*(2./3.)*C[i]*(eta_stress(sigma-X)%Ir05());
            }    
            //compute Phi and the derivatives
            Phi(0) = Mises_stress(sigma-X) - Hp - sigmaY;        
            Lambdap = eta_stress(sigma-X);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = eta_stress(sigma-X) - D[i]*a_N[i];
            }
            break;   
        }     
        case 1: {//Hill
            dPhidsigma = dHill_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                dPhida_N[i] = -1.*(2./3.)*C[i]*(dHill_stress(sigma-X, criteria_params)%Ir05());
            }    
            //compute Phi and the derivatives
            Phi(0) = Hill_stress(sigma-X, criteria_params) - Hp - sigmaY;        
            Lambdap = dHill_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dHill_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }
            break;
        }
        case 2: {//DFA
            dPhidsigma = dDFA_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                dPhida_N[i] = -1.*(2./3.)*C[i]*(dDFA_stress(sigma-X, criteria_params)%Ir05());
            }    
            //compute Phi and the derivatives
            Phi(0) = DFA_stress(sigma-X, criteria_params) - Hp - sigmaY;        
            Lambdap = dDFA_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dDFA_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }
            break;
        }
        case 3: {//ANI
            dPhidsigma = dAni_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                dPhida_N[i] = -1.*(2./3.)*C[i]*(dAni_stress(sigma-X, criteria_params)%Ir05());
            }    
            //compute Phi and the derivatives
            Phi(0) = Ani_stress(sigma-X, criteria_params) - Hp - sigmaY;        
            Lambdap = dAni_stress(sigma-X, criteria_params);
            for (int i=0; i<N_kin_hard; i++) {
                Lambdaa_N[i] = dAni_stress(sigma-X, criteria_params) - D[i]*a_N[i];
            }
            break;
        }
        default: {
            break;
        }
    }
//        dPhida = 0.*(eta_stress(sigma - X)%Ir05());
        
        kappa_j[0] = L*Lambdap;
        K(0,0) = dPhidp ;
        for (int i=0; i<N_kin_hard; i++) {
             K(0,0) += sum(dPhida_N[i]%Lambdaa_N[i]);
        }                

        B(0, 0) = -1.*sum(dPhidsigma%kappa_j[0]) + K(0,0);
        Y_crit(0) = sigmaY;
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        
        s_j(0) += ds_j(0);

        for (int i=0; i<N_iso_hard; i++) {        
            Hp += b[i]*(Q[i]-Hp)*ds_j(0);
        }            

        EP = EP + ds_j(0)*Lambdap;

        X = zeros(6);
        for (int i=0; i<N_kin_hard; i++) {
            a_N[i] = a_N[i] + ds_j(0)*Lambdaa_N[i];
            X_N[i] += ds_j(0)*(2./3.)*C[i]*(Lambdaa_N[i]%Ir05());             
            X+= X_N[i];
        }                

        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        sigma = el_pred(L, Eel, ndi);
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    double Dp = Ds_j[0];

    std::vector<vec> Da_N(N_kin_hard);
    for (int i=0; i<N_kin_hard; i++) {
        Da_N[i] = a_N[i] - a_N_start[i];
    }                
        
    //Computation of the tangent modulus
    mat Bhat = zeros(1, 1);
    Bhat(0, 0) = sum(dPhidsigma%kappa_j[0]) - K(0,0);

    vec op = zeros(1);
    mat delta = eye(1,1);

    for (int i=0; i<1; i++) {
        if(Ds_j[i] > sim_iota)
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

    Lt = L - (kappa_j[0]*P_epsilon[0].t());
    
    double A_p = -Hp;

    std::vector<vec> A_a_N(N_kin_hard);
    for (int i=0; i<N_kin_hard; i++) {
        A_a_N[i] = -X_N[i];
    }                
    
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP) + 0.5*(A_p_start + A_p)*Dp;
    for (int i=0; i<N_kin_hard; i++) {
        Dgamma_loc+=0.5*sum((A_a_N_start[i] + A_a_N[i])%Da_N[i]);
    }                
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    for (int i=0; i<N_kin_hard; i++) {
        Wm_r+= -0.5*sum((A_a_N_start[i] + A_a_N[i])%Da_N[i]);
    }    
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
    
    statev(8) = Hp;

    for (int i=0; i<N_kin_hard; i++) {
        statev(8+i*12+1) = a_N[i](0);
        statev(8+i*12+2) = a_N[i](1);
        statev(8+i*12+3) = a_N[i](2);
        statev(8+i*12+4) = a_N[i](3);
        statev(8+i*12+5) = a_N[i](4);
        statev(8+i*12+6) = a_N[i](5);

        statev(8+i*12+7) = X_N[i](0);
        statev(8+i*12+8) = X_N[i](1);
        statev(8+i*12+9) = X_N[i](2);
        statev(8+i*12+10) = X_N[i](3);                
        statev(8+i*12+11) = X_N[i](4);
        statev(8+i*12+12) = X_N[i](5);
    }
}

    
} //namespace simcoon
