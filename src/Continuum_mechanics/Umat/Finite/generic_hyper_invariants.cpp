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

///@file generic_hyper_invariants.cpp
///@brief User subroutine for hyperelastic materials using isochoric invariants
///@version 1.0

#include <iostream>
#include <fstream>
#include <map>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_generic_hyper_invariants(const std::string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
{  	

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);    
	vec sigma_start = sigma;

    //definition of the Right Cauchy-Green tensor
    mat b = L_Cauchy_Green(F1);

    double J;
    try {
        J = det(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside umat_generic_hyper_invariants.");
    }     
    vec I_bar = isochoric_invariants(b, J);

    static const std::map<string, HyperPotential> list_potentials = {
        {"NEOHC", HyperPotential::NEOHC}, {"MOORI", HyperPotential::MOORI},
        {"YEOHH", HyperPotential::YEOHH}, {"ISHAH", HyperPotential::ISHAH},
        {"GETHH", HyperPotential::GETHH}, {"SWANH", HyperPotential::SWANH}};

    auto it_potential = list_potentials.find(umat_name);
    if (it_potential == list_potentials.end()) {
        throw std::invalid_argument("The choice of hyperelastic potential could not be found in the simcoon library: " + umat_name);
    }
    const hyper_invariants_dW dW = hyper_potential_derivatives(it_potential->second, props, I_bar, J);
    
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

    hyper_invariants_response(dW, b, J, F1, sigma, Lt);

    if(start) {
        L = Lt;
    }

    //Computation of the mechanical and thermal work quantities.
    // Kirchhoff work per reference volume: tau:d(lnV) with tau = J*sigma (see saint_venant).
    double J0 = det(F0);
    Wm   += 0.5*sum((J0*sigma_start + J*sigma)%Detot);
    Wm_r += 0.5*sum((J0*sigma_start + J*sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
