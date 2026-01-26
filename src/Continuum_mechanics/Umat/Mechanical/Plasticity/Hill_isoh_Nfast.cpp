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

///@file Iso_Nfast.cpp
///@brief User subroutine for Multiple plasticity model in 3D case
///@author Chemisky, Chatzigeorgiou
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>

using namespace std;
using namespace arma;

///@brief The Multiple plastic series model requires a number of constants
//      -------------------
///@brief      props(0) = E0            Young modulus
///@brief      props(1) = nu0           Poisson ratio
///@brief      props(2) = alpha_iso     CTE
///@brief      props(3) = N_plas        Number of plastic branches
///@brief      props(4+i*3) = sigmaY(i)   yield stress limit sigmaY of branch i
///@brief      props(4+i*3+1) = k(i)   hardening parameter k of branch i
///@brief      props(4+i*3+2) = m(i)   exponent m of branch i

///@brief Number of statev required for thermoelastic constitutive law

namespace simcoon {
    
void umat_plasticity_hill_isoh_CCP_N(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{
    
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double alpha_iso = props(2);
    int N_plas = int(props(3));
    
    vec sigmaY = zeros(N_plas);
    vec k = zeros(N_plas);
    vec m = zeros(N_plas);
    vec F_hill = zeros(N_plas);
    vec G_hill = zeros(N_plas);
    vec H_hill = zeros(N_plas);
    vec L_hill = zeros(N_plas);
    vec M_hill = zeros(N_plas);
    vec N_hill = zeros(N_plas);
    
    for (int i=0; i<N_plas; i++) {
        sigmaY(i) = props(4+i*9);
        k(i) = props(4+i*9+1);
        m(i) = props(4+i*9+2);
        F_hill = props(4+i*9+3);
        G_hill = props(4+i*9+4);
        H_hill = props(4+i*9+5);
        L_hill = props(4+i*9+6);
        M_hill = props(4+i*9+7);
        N_hill = props(4+i*9+8);
    }

    std::vector<vec> Hill_params(N_plas);
    for (int i=0; i<N_plas; i++) {
        Hill_params[i] = {F_hill[i],G_hill[i],H_hill[i],L_hill[i],M_hill[i],N_hill[i]};
    }
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Define the elastic stiffness
    mat L = L_iso(E, nu, "Enu");
    mat M = M_iso(E, nu, "Enu");
    ///@brief Temperature initialization
    double T_init = statev(0);
    
    //From the statev to the internal variables
    vec EP = zeros(6);
    EP(0) = statev(1);
    EP(1) = statev(2);
    EP(2) = statev(3);
    EP(3) = statev(4);
    EP(4) = statev(5);
    EP(5) = statev(6);
    
    std::vector<vec> EP_branch(N_plas);
    vec p = zeros(N_plas);
    
    for (int i=0; i<N_plas; i++) {
        p(i) = statev(i*7+7);
        EP_branch[i] = zeros(6);
        EP_branch[i](0) = statev(i*7+7+1);
        EP_branch[i](1) = statev(i*7+7+2);
        EP_branch[i](2) = statev(i*7+7+3);
        EP_branch[i](3) = statev(i*7+7+4);
        EP_branch[i](4) = statev(i*7+7+5);
        EP_branch[i](5) = statev(i*7+7+6);
    }
    
    //Rotation of internal variables (tensors)
    EP = rotate_strain(EP, DR);
    for (int i=0; i<N_plas; i++) {
        EP_branch[i] = rotate_strain(EP_branch[i], DR);
    }
    
    vec sigma_start = sigma;
    vec EP_start = sigma;
    
    if(start) { //Initialization
        T_init = T;
        EP = zeros(6);
        for (int i=0; i<N_plas; i++) {
            EP_branch[i] = zeros(6);
        }
        sigma = zeros(6);
        sigma_start = zeros(6);
    }
    
    vec Hp = zeros(N_plas);
    vec dHpdp = zeros(N_plas);

    for (int i=0; i<N_plas; i++) {
        if (p[i] > simcoon::iota)	{
            dHpdp[i] = m[i]*k[i]*pow(p[i], m[i]-1.0);
            Hp[i] = k[i]*pow(p[i], m[i]);
        }
        else {
            dHpdp[i] = 0.;
            Hp[i] = 0.;
        }
    }
    std::vector<double> A_p_start(N_plas);
    for (int i = 0; i < N_plas; i++) {
        A_p_start[i] = -Hp[i];
    }
    
    //Variables required for the loop
    vec s_j = p;
    vec Ds_j = zeros(N_plas);
    vec ds_j = zeros(N_plas);
    
    ///Elastic prediction - Accounting for the thermal prediction
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init) - EP;
    sigma = el_pred(L, Eel, ndi);

    //Define the plastic function and the stress
    vec Phi = zeros(N_plas);
    mat B = zeros(N_plas,N_plas);
    vec Y_crit = zeros(N_plas);
    
    vec dPhidp = zeros(N_plas);
    //Compute the explicit flow direction
    std::vector<vec> Lambdap(N_plas);
    std::vector<vec> kappa_j(N_plas);
    std::vector<vec> dPhidsigma(N_plas);
    mat K = zeros(N_plas,N_plas);
    
    for (int i=0; i<N_plas; i++) {
        Lambdap[i] = dHill_stress(sigma,Hill_params[i]);
        kappa_j[i] = L*Lambdap[i];
        dPhidsigma[i] = zeros(6);
    }
    
    //Loop parameters
    int compteur = 0;
    double error = 1.;

    //Loop
    for (compteur = 0; ((compteur < maxiter_umat) && (error > precision_umat)); compteur++) {
        
        p = s_j;
        for (int i=0; i<N_plas; i++) {
            if (p[i] > simcoon::iota)	{
                dHpdp[i] = m[i]*k[i]*pow(p[i], m[i]-1.0);
                Hp[i] = k[i]*pow(p[i], m[i]);
            }
            else {
                dHpdp[i] = 0.;
                Hp[i] = 0.;
            }
            
            Phi(i) = Hill_stress(sigma,Hill_params[i]) - Hp[i] - sigmaY[i];
            Lambdap[i] = dHill_stress(sigma,Hill_params[i]);
            
            dPhidsigma[i] = Lambdap[i];
            dPhidp[i] = -1.*dHpdp[i];
            
            kappa_j[i] = L*Lambdap[i];
            K(i,i) = dPhidp[i];
            for (int j=0; j<N_plas; j++) {
                B(i,j) = -1.*sum(dPhidsigma[i]%kappa_j[j]) + K(i,j);
            }
            Y_crit(i) = sigmaY[i];
        }
        
        Fischer_Burmeister_m(Phi, Y_crit, B, Ds_j, ds_j, error);
        s_j += ds_j;
        
        for (int i=0; i<N_plas; i++) {
            EP_branch[i] += ds_j(i)*Lambdap[i];
            EP += ds_j(i)*Lambdap[i];
        }
        
        //the stress is now computed using the relationship sigma = L(E-Ep)
        Eel = Etot + DEtot - alpha*(T + DT - T_init) - EP;
        sigma = el_pred(L, Eel, ndi);
    }
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    vec DEP = EP - EP_start;
    vec Dp = Ds_j;
    
    //Computation of the tangent modulus
    mat Bhat = zeros(N_plas, N_plas);
    for (int i=0; i<N_plas; i++) {
        for (int j=0; j<N_plas; j++) {
            Bhat(i, j) = sum(dPhidsigma[i]%kappa_j[j]) - K(i,j);
        }
    }
    
    vec op = zeros(N_plas);
    mat delta = eye(N_plas,N_plas);
    
    for (int i=0; i<N_plas; i++) {
        if(Ds_j[i] > simcoon::iota)
            op(i) = 1.;
    }
    
    mat Bbar = zeros(N_plas,N_plas);
    for (int i = 0; i < N_plas; i++) {
        for (int j = 0; j < N_plas; j++) {
            Bbar(i, j) = op(i)*op(j)*Bhat(i, j) + delta(i,j)*(1.-op(i)*op(j));
        }
    }
    
    mat invBbar = zeros(N_plas, N_plas);
    mat invBhat = zeros(N_plas, N_plas);
    invBbar = inv(Bbar);
    for (int i = 0; i < N_plas; i++) {
        for (int j = 0; j < N_plas; j++) {
            invBhat(i, j) = op(i)*op(j)*invBbar(i, j);
        }
    }
    
    std::vector<vec> P_epsilon(N_plas);

    Lt = L;
    for (int i = 0; i < N_plas; i++) {
        P_epsilon[i] = zeros(6);
        for (int j = 0; j < N_plas; j++) {
            P_epsilon[i] += invBhat(j, i)*(L*dPhidsigma[j]);
        }
        Lt = Lt - (kappa_j[i]*P_epsilon[i].t());
    }
    
    std::vector<double> A_p(N_plas);
    for (int i = 0; i < N_plas; i++) {
        A_p[i] = -Hp[i];
    }
    
    double Dgamma_loc = 0.5*sum((sigma_start+sigma)%DEP);
    for (int i = 0; i < N_plas; i++) {
        Dgamma_loc += 0.5*(A_p_start[i] + A_p[i])*Dp[i];
    }
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot-DEP));
    for (int i = 0; i < N_plas; i++) {
        Wm_ir += -0.5*(A_p_start[i] + A_p[i])*Dp[i];
    }
    Wm_d += Dgamma_loc;
    
    ///@brief statev evolving variables
    //statev
    statev(0) = T_init;
    
    statev(1) = EP(0);
    statev(2) = EP(1);
    statev(3) = EP(2);
    statev(4) = EP(3);
    statev(5) = EP(4);
    statev(6) = EP(5);
    
    for (int i=0; i<N_plas; i++) {
        statev(i*7+7) = p(i);
        statev(i*7+7+1) = EP_branch[i](0);
        statev(i*7+7+2) = EP_branch[i](1);
        statev(i*7+7+3) = EP_branch[i](2);
        statev(i*7+7+4) = EP_branch[i](3);
        statev(i*7+7+5) = EP_branch[i](4);
        statev(i*7+7+6) = EP_branch[i](5);
    }
}

} //namespace simcoon

