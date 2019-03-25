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
///@brief User subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/Neo_hookean.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_elasticity_iso(const mat &F, vec &sigma, mat &Lt, mat &L, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
{  	

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);
    
    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double alpha = props(2);

    mat C = R_Cauchy_Green(F);
    mat B = L_Cauchy_Green(F);
    
    vec lambda = eig_sym(C);
    
    double I1 = pow(lambda(2),2.) + pow(lambda(1),2.) + pow(lambda(0),2.); //ascending order
    double J = lambda(2)*lambda(1)*lambda(0);
    double I1_bar = pow(J,-2./3.)*I1;
    
    double W = C10*(I1_bar-3.) + D1*pow(J-2.,2.);

    double p = -2.*D1*J*(J-1);

    mat I = eye(3,3);
    mat PKII = (-2./3.)*I1_bar*inc(C) (1./J)*(-p*I)+2.*C1*dev()
    mat sigma_Cauchy = (1./J)*(-p*I)+2.*C1*dev();
    vec sigma = m2t_stress(sigma_Cauchy);

    //Elastic stiffness tensor
    L = L_iso(E, nu, "Enu");
    
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
