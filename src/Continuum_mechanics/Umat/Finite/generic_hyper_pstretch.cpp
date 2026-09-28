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

///@file generic_hyper_pstretch.cpp
///@brief User subroutine for hyperelastic materials using isochoric principal stretches
///@version 1.0

#include <iostream>
#include <fstream>
#include <armadillo>
#include <math.h>
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
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void umat_generic_hyper_pstretch(const std::string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode)
{  	

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
        throw simcoon::exception_det("Error in det function inside umat_generic_hyper_pstretch.");
    }     

    vec dWdlambda_bar = zeros(3);
    mat dW2dlambda_bar2 = zeros(3,3);
    vec lambda_bar = zeros(3);
    mat n_pvectors = zeros(3,3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", J);

    // single principal-stretch potential so far; the potential itself is documented in the .hpp
    if (umat_name != "OGDEN") {
        throw std::invalid_argument("The choice of hyperelastic potential could not be found in the simcoon library: " + umat_name);
    }

    if (nprops < 2) {
        throw std::invalid_argument("OGDEN expects props = {N, kappa, mu_1, alpha_1, ...}, got nprops = " + std::to_string(nprops));
    }
    int N_Ogden = int(props(0));
    double kappa = props(1);
    if (N_Ogden < 1) {
        throw std::invalid_argument("OGDEN expects N >= 1 as props(0) (a zero/negative N would silently yield a volumetric-only law), got N = " + std::to_string(N_Ogden));
    }
    if (nprops < 2 + 2*N_Ogden) {
        throw std::invalid_argument("OGDEN expects nprops >= 2 + 2*N, got nprops = " + std::to_string(nprops) + " for N = " + std::to_string(N_Ogden));
    }

    for (int i=0; i<N_Ogden; i++) {
        const double mu_i = props(2+i*2);
        const double alpha_i = props(2+i*2+1);
        if (std::abs(alpha_i) < simcoon::iota) {
            throw std::invalid_argument("OGDEN exponent alpha_" + std::to_string(i+1) + " must be nonzero (c_i = 2 mu_i / alpha_i)");
        }
        const double c_i = 2.*mu_i/alpha_i;
        const vec p = pow(lambda_bar, alpha_i-1.);
        dWdlambda_bar += c_i*p;
        dW2dlambda_bar2.diag() += c_i*(alpha_i-1.)*(p/lambda_bar);
    }

    // U(J): kappa (J ln J - J + 1) unless props(2 + 2N) = 1 selects kappa/2 (J - 1)^2
    double dUdJ = 0., dU2dJ2 = 0.;
    volumetric_derivatives(volumetric_potential_of(props, 2 + 2*N_Ogden), kappa, J, dUdJ, dU2dJ2);
    
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

    // Kirchhoff-native: the route stress IS tau, Cauchy is tau/J formed at the output boundary.
    // The 3-argument tau_iso_hyper_pstretch needs no J at all.
    mat m_tau_iso = tau_iso_hyper_pstretch(dWdlambda_bar, lambda_bar, N_projectors);
    mat m_tau_vol = tau_vol_hyper(dUdJ, b, J);
    mat m_tau = m_tau_iso + m_tau_vol;
    sigma = t2v_stress(m_tau);

    // These return the spatial elasticity c = (1/J) d(L_v tau)/dD, not d(L_v sigma)/dD (they
    // differ by sigma (x) I). The J that turns c into the Kirchhoff-Lie tangent is applied once,
    // below. L_iso carries an explicit 1/J and L_vol is J-free by cancellation, so they must be
    // scaled as a sum -- see the note on L_vol_hyper.
    mat Lt_iso = L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, lambda_bar, n_pvectors, J);
    mat Lt_vol = L_vol_hyper(dUdJ, dU2dJ2, b, J);
    mat Lt_spatial = Lt_iso + Lt_vol;

    // Standardize to the canonical box convention Lt = d(tau_hat)/d(De) -- see generic_hyper_invariants.
    Lt = Dtau_LieDD_2_DtauDe_corate(J*Lt_spatial, corate_type, F1, m_tau);

    if(start) {
        L = Lt;
    }

/*    cout << "L = " << L << endl;
    cout << "Lt = " << Lt << endl;
    cout << "Lt_iso = " << Lt_iso << endl;    
    cout << "Lt_vol = " << Lt_vol << endl;        
    cout << "eig(Lt)" << eig_sym(Lt);
*/
    
    //Computation of the mechanical and thermal work quantities.
    // Kirchhoff work per reference volume: tau:d(lnV). Both ends are ALREADY tau -- the kernel
    // is Kirchhoff-native and the stored state is tau_n -- so no J enters here any more.
    Wm   += 0.5*sum((sigma_start + sigma)%Detot);
    Wm_r += 0.5*sum((sigma_start + sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
