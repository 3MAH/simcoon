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

///@file Elastic_orthotropic.cpp
///@brief User subroutine for ortothropic elastic materials in 3D case

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_orthotropic.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The thermomechanical elastic orthoptropic UMAT requires 13 constants:
///@brief props[0] : density
///@brief props[1] : specific heat capacity
///@brief props[2-4] : 3 Young modulus
///@brief props[5-7] : 3 Poisson ratio
///@brief props[8-10] : 3 Shear ratio ??
///@brief props[11] : CTE x
///@brief props[12] : CTE y
///@brief props[13] : CTE z

///@brief No statev is required for thermoelastic constitutive law

void umat_elasticity_ortho_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{  	

    UNUSED(Etot);
    UNUSED(DR);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(T);
    UNUSED(Time);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
	//From the props to the material properties
    double rho = props(0);
    double c_p = props(1);
    double Ex = props(2);
	double Ey = props(3);
	double Ez = props(4);
	double nuxy = props(5);
	double nuxz = props(6);
    double nuyz = props(7);
	double Gxy = props(8);
	double Gxz = props(9);
    double Gyz = props(10);
	double alphax = props(11);
	double alphay = props(12);
    double alphaz = props(13);
    
    double T_init = statev(0);
    
    //definition of the CTE tensor
    vec alpha = zeros(6);
    alpha(0) = alphax;
    alpha(1) = alphay;
    alpha(2) = alphaz;    
    
    //Elastic stiffness tensor and thermal tensor
    dSdE = L_ortho(Ex,Ey,Ez,nuxy,nuxz,nuyz,Gxy,Gxz,Gyz, "EnuG");
    dSdT = -1.*dSdE*alpha;
    
    if(start) { //Initialization
        T_init = T;
        sigma = zeros(6);
        
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
    
    vec sigma_start = sigma;
    
    //Compute the elastic strain and the related stress
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init);
    sigma = el_pred(dSdE, Eel, ndi);
    
    //Computation of the increments of variables
    vec Dsigma = sigma - sigma_start;
    
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
    
    if(DTime < 1.E-12) {
        r = 0.;
        drdE = zeros(6);
        drdT = 0.;
    }
    else {
        Gamma_epsilon = zeros(6);
        Gamma_theta = 0.;
        
        N_epsilon = -1./DTime*(T + DT)*(dSdE*alpha);
        N_theta = -1./DTime*(T + DT)*sum(dSdT%alpha) -1.*Deta/DTime - rho*c_p*(1./DTime);
        
        drdE = N_epsilon;// + Gamma_epsilon;
        drdT = N_theta;// + Gamma_theta;
        
        r = sum(N_epsilon%DEtot) + N_theta*DT + sum(Gamma_epsilon%DEtot) + Gamma_theta*DT;
    }
    
    double Dgamma_loc = 0.;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot));
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir = (T+0.5*DT)*Deta_ir;
    
    statev(0) = T_init;
}

} //namespace simcoon
