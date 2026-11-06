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

///@file saint_venant.cpp
///@brief User subroutine for Isotropic elastic materials in 3D case
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
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
#include <simcoon/Continuum_mechanics/Umat/Finite/saint_venant.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void umat_saint_venant(const string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);
    
    //From the props to the material properties
    double E = props(0);
    double nu = props(1);
    double alpha = props(2);

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
    
    //definition of the Right Cauchy-Green tensor
    vec Eel = t2v_strain(Green_Lagrange(F1));
        
    //Compute the PKII stress and then the Kirchhoff stress (the route stress)
    mat S = v2t_stress(el_pred(L, Eel));
    sigma = t2v_stress(PKII2Kirchoff(S, F1));

    // not box_DtauDe_from_dSdE: it expects Cauchy, and would square J on tau
    Lt = DSDE_2_DtauDe_corate(L, corate_type, F1, v2t_stress(sigma));

    // Kirchhoff work per reference volume, tau : d(lnV): no J
    Wm   += 0.5*sum((sigma_start + sigma)%Detot);
    Wm_r += 0.5*sum((sigma_start + sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
