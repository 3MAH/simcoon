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
#include <map>
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
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_pstretch.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_generic_hyper_pstretch(const std::string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
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

    double dUdJ = 0.;
    double dU2dJ2 = 0.;

    double J = det(F1);  

    vec dWdlambda_bar = zeros(3);
    mat dW2dlambda_bar2 = eye(3,3);    
    vec lambda_bar = zeros(3);
    vec n_pvectors = zeros(3);
    std::vector<mat> N_projectors(3);
    isochoric_pstretch(lambda_bar, n_pvectors, N_projectors, b, "b", J);

    std::map<string, int> list_potentials;
    list_potentials = {{"OGDEN",0}};

    switch (list_potentials[umat_name]) {
        case 0: {

            int N_Ogden = int(props(0));
            double kappa = props(1);            
            
            vec mu = zeros(N_Ogden);
            vec alpha = zeros(N_Ogden);
            
            for (int i=0; i<N_Ogden; i++) {
                mu(i) = props(2+i*2);
                alpha(i) = props(2+i*2+1);
            }

            // \f$ W = \sum_i=0^N{ \frac{2 \mu_i}{\alpha_i^2} \left( \lambda_1^{\alpha_i} + \lambda_2^{\alpha_i} + \lambda_3^{\alpha_i} \right) } \f$ 

            dWdlambda_bar = 2./lambda_bar;
            for (int i=0; i<N_Ogden; i++) {
                dWdlambda_bar = mu(i)/alpha(i)*(dWdlambda_bar%pow(lambda_bar,alpha(i)));
            }
            
            for (int i=0; i<N_Ogden; i++) {
                for(int a=0; a<3; a++) {
                    dW2dlambda_bar2(a,a) = 2./(pow(lambda_bar(a),2.)*mu(i)*((alpha(i)-1.)/alpha(i))*pow(lambda_bar(a),alpha(i)));
                }
            }

            dUdJ = kappa*log(J);
            dU2dJ2 = kappa/J;
            break;
        }
        default: {
            cout << "Error: The choice of hyperelastic potential could not be found in the simcoon library :" << umat_name
             << "\n";
                exit(0);
        }
    }
    
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

    mat m_sigma_iso = sigma_iso_hyper_pstretch(dWdlambda_bar, lambda_bar, N_projectors, J);
    mat m_sigma_vol = sigma_vol_hyper(dUdJ, b, J);    
    mat m_sigma = m_sigma_iso + m_sigma_vol;
    sigma = t2v_stress(m_sigma);  

    mat Lt_iso = L_iso_hyper_pstretch(dWdlambda_bar, dW2dlambda_bar2, lambda_bar, n_pvectors, J);  
    mat Lt_vol = L_vol_hyper(dUdJ, dU2dJ2, b, J);
    Lt = Lt_iso + Lt_vol;

    if(start) {
        L = Lt;
    }

/*    cout << "L = " << L << endl;
    cout << "Lt = " << Lt << endl;
    cout << "Lt_iso = " << Lt_iso << endl;    
    cout << "Lt_vol = " << Lt_vol << endl;        
    cout << "eig(Lt)" << eig_sym(Lt);
*/
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_r += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
