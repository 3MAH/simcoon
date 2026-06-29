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

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_saint_venant(const string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{

    UNUSED(umat_name);
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
        
    //Compute the PKII stress and then the Cauchy stress
    mat S = v2t_stress(el_pred(L, Eel, ndi));
    sigma = t2v_stress(PKII2Cauchy(S, F1));

    // Box tangent in the canonical convention Lt = d(tau_hat)/d(De). For SVK the elastic
    // stiffness L IS the material tangent dS/dE, so map it straight to the box tangent.
    Lt = box_DtauDe_from_dSdE(L, F1, sigma);
        
    //Computation of the mechanical and thermal work quantities.
    // Wm is the Kirchhoff work per REFERENCE volume: tau:d(lnV), conjugate to the log-strain
    // increment Detot. The box stress 'sigma' is the true Cauchy, so use tau = J*sigma (NOT
    // sigma, which would give per-current-volume work) -- consistent with the small-strain
    // (Kirchhoff-route) boxes whose route stress IS tau.
    double J0 = det(F0), J1 = det(F1);
    Wm   += 0.5*sum((J0*sigma_start + J1*sigma)%Detot);
    Wm_r += 0.5*sum((J0*sigma_start + J1*sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
