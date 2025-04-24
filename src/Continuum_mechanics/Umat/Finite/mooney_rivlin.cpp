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
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/mooney_rivlin.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_mooney_rivlin(const vec &Etot, const vec &DEtot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
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

    double C_10 = E/(4.*(1+nu));
    double D_1 = 6.*(1-2.*nu)/E;
    
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
    
    //Invariants of C
    double I1 = trace(C); // pow(lambda_alpha(2),2.) + pow(lambda_alpha(1),2.) + pow(lambda_alpha(0),2.); //ascending order
    double J;
    try {
        J = det(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside umat_mooney_rivlin.");
    }    
    double I1_bar = pow(J,-2./3.)*I1;
    
    double W = C_10*(I1_bar-3.) + (1./D_1)*pow(J-1.,2.);

    mat invC;
    try {
        invC = inv(C);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside umat_mooney_rivlin.");
    }      
    mat I = eye(3,3);

    //Compute the PKII stress and then the Cauchy stress
    mat S = (-2./3.)*C_10*I1_bar*invC + 2.*C_10*(1./J)*pow(J,-2./3.)*I + (2./D_1)*(J-1)*J*invC;
    mat sigma_Cauchy = PKII2Cauchy(S, F1, J);
    sigma = t2v_stress(sigma_Cauchy);
	
//    L = (-2./3.)*C_10*pow(J,-2./3.)*sym_dyadic(invC,I)+(2./9.)*C_10*I1_bar*sym_dyadic(invC,invC)-(2./3.)*C_10*I1_bar*dinvSdSsym(C)
//    -(2./3.)*C_10*pow(J,-2./3.)*sym_dyadic(I,invC)
//    +(1./D_1)*(J-1.)*J*sym_dyadic(invC,invC)+(2./D_1)*(J-1)*J*dinvSdSsym(C);
    L = (-2./3.)*C_10*pow(J,-2./3.)*dyadic(invC,I)+(2./9.)*C_10*I1_bar*auto_dyadic(invC)-(2./3.)*C_10*I1_bar*dinvSdSsym(C)
    -(2./3.)*C_10*pow(J,-2./3.)*dyadic(I,invC)
    +(1./D_1)*(J-1.)*J*auto_dyadic(invC)+(2./D_1)*(J-1)*J*dinvSdSsym(C);

    
    if((solver_type == 0)||(solver_type==2)) {
        Lt = L;
	}
    
    //Computation of the mechanical and thermal work quantities
    /*
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_ir += 0.;
    Wm_d += 0.;
    */
    
    statev(0) = T_init;
}

} //namespace simcoon
