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
#include <math.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/Neo_hookean_comp.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_generic_hyper(const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
{  	

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);
    
    int N_invariants = int(props(0));
    vec invariants = zeros(N_invariants);
    //From the props to the material properties
    for (unsigned int i = 0; i < N_invariants; i++)
    {
        invariants(i) = props(i+1);
    }
    
    double phi_1 = 0;
    double phi_2 = 0;
    double phi_3 = 0;        
    
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
    
    //definition of the Right Cauchy-Green tensor
    mat C = R_Cauchy_Green(F1);
    
    //definition of the Right Cauchy-Green tensor
    mat Eel = Green_Lagrange(F1);    
        
    //Compute the PKII stress and then the Cauchy stress
    mat S = el_pred(L, Eel, ndi);    
    sigma = t2v_stress(PKII2Cauchy(S, F1));

    Lt = DSDE_2_DtauDe(L, get_BBBB(F1), F1, v2t_stress(sigma));    
        
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_r += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
