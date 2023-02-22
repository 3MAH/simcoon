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
#include <simcoon/Continuum_mechanics/Umat/Mechanical/Elasticity/elastic_orthotropic.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 9 constants:
///@brief props[0-2] : 3 Young modulus
///@brief props[3-5] : 3 Poisson ratio
///@brief props[6-8] : 3 Shear ratio ??
///@brief props[9] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_elasticity_ortho(const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
{  	

    UNUSED(Etot);
    UNUSED(DR);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);
    
    //From the props to the material properties
    double Ex = props(0);
    double Ey = props(1);
    double Ez = props(2);
    double nuxy = props(3);
    double nuxz = props(4);
    double nuyz = props(5);
    double Gxy = props(6);
    double Gxz = props(7);
    double Gyz = props(8);
    double alphax = props(9);
    double alphay = props(10);
    double alphaz = props(11);

    //Elastic stiffness tensor
    L = L_ortho(Ex,Ey,Ez,nuxy,nuxz,nuyz,Gxy,Gxz,Gyz, "EnuG");    
    
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        sigma = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
	vec sigma_start = sigma;
	
	//definition of the CTE tensor
	vec alpha = zeros(6);
	alpha(0) = alphax;
	alpha(1) = alphay;
	alpha(2) = alphaz;
    
	//Compute the elastic strain and the related stress	
    vec Eel = Etot + DEtot - alpha*(T+DT-T_init);
    sigma = el_pred(L, Eel, ndi);
    
    if((solver_type == 0)||(solver_type==2)) {
        Lt = L;
	}
    else if(solver_type == 1) {
        sigma_in = zeros(6);
    }
        
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
