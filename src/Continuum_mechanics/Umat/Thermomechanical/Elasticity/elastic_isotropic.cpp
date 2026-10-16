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

///@file elastic_isotropic.cpp
///@brief User thermomechanical subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>
#include <simcoon/Continuum_mechanics/Umat/tangent_assembly.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_elasticity_iso_T(const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, mat &dSdT, mat &drdE, mat &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{  	

    UNUSED(Etot);
    UNUSED(DR);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(T);
    UNUSED(Time);
    UNUSED(tnew_dt);    
    
	//From the props to the material properties
    double rho = props(0);
    double c_p = props(1); // Make sure c_p has been identified at T = T_init
    double E = props(2);
	double nu = props(3);
	double alpha_iso = props(4);
    
    double T_init = statev(0);
    
    //definition of the CTE tensor
    vec alpha = alpha_iso*Ith();
    
    //Elastic stiffness tensor and thermal tensor
	dSdE = L_iso(E, nu, "Enu");
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
    sigma = el_pred(dSdE, Eel, ndi, nshr);
    
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
    
    // heat of the increment R = -(T+DT) Deta (no dissipation) and its exact derivatives
    const double T1 = T + DT;
    heat_source_increment(0., zeros(6), 0., Deta, dSdE.t()*alpha, c_0/T1 + sum(dSdT%alpha), T1, r, drdE, drdT);
    
    double Dgamma_loc = 0.;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%(DEtot));
    Wm_ir += 0.;
    Wm_d += Dgamma_loc;
    
    Wt += (T+0.5*DT)*Deta;
    Wt_r += (T+0.5*DT)*Deta_r;
    Wt_ir += (T+0.5*DT)*Deta_ir;
    
    statev(0) = T_init;
}
    
} //namespace simcoon
