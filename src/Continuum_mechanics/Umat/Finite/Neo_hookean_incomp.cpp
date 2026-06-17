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
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/derivatives.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/neo_hookean_incomp.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_neo_hookean_incomp(const string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode)
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
        throw simcoon::exception_det("Error in det function inside umat_neo_hookean_incomp.");
    }   
    double I1_bar = pow(J,-2./3.)*I1;
    
    double W = C_10*(I1_bar-3.) + (1./D_1)*pow(J-1.,2.);

    mat invC;
    try {
        invC = inv(C);
    } catch (const std::runtime_error &e) {
        cerr << "Error in inv: " << e.what() << endl;
        throw simcoon::exception_inv("Error in inv function inside umat_neo_hookean_incomp.");
    }
    mat I = eye(3,3);

    //Compute the PKII stress and then the Cauchy stress.
    // S = 2 dW/dC for W = C_10*(I1_bar-3) + (1/D_1)*(J-1)^2, with I1_bar = J^(-2/3) tr(C):
    //   S_iso = 2 C_10 J^(-2/3) I - (2/3) C_10 I1_bar invC ;  S_vol = (2/D_1)(J-1) J invC.
    // (The previous deviatoric I-term carried a stray 1/J -> J^(-5/3) instead of J^(-2/3),
    //  a ~ (1-1/J) stress error; fixed so stress and the rebuilt tangent both match dS/dE.)
    mat S = (-2./3.)*C_10*I1_bar*invC + 2.*C_10*pow(J,-2./3.)*I + (2./D_1)*(J-1)*J*invC;
    mat sigma_Cauchy = PKII2Cauchy(S, F1, J);
    sigma = t2v_stress(sigma_Cauchy);
	
    // Tangent. The previous hand-built dyadic material tangent did NOT match dS/dE
    // (FD ~100% off) -- a pre-existing bug. Rebuild it from the same Neo-Hookean
    // potential W = C_10*(I1_bar-3) + (1/D_1)*(J-1)^2 using the verified invariant
    // hyperelastic machinery (Cauchy/Oldroyd-Lie spatial elasticity), then standardize
    // to the canonical box convention Lt = d(tau_hat)/d(De) (Kirchhoff log-rate, XBM),
    // identical to generic_hyper_invariants / saint_venant.
    mat b = L_Cauchy_Green(F1);
    double dWdI_1_bar = C_10;       // dW/dI1_bar; dW/dI2_bar = 0 (no I2 term), all 2nd deviatoric derivs = 0
    double dUdJ   = (2./D_1)*(J-1.);
    double dU2dJ2 = 2./D_1;
    mat Lt_spatial = L_iso_hyper_invariants(dWdI_1_bar, 0., 0., 0., 0., b, J) + L_vol_hyper(dUdJ, dU2dJ2, b, J);
    mat dSdE = Dtau_LieDD_2_DSDE(J*Lt_spatial, F1);
    Lt = DSDE_2_DtauDe(dSdE, get_BBBB(F1), F1, J*v2t_stress(sigma));

    if(start) {
        L = Lt;
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
