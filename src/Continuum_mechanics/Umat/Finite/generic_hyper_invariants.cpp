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
///@brief User subroutine for hyperelastic materials using isochoric invariants
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
#include <simcoon/Continuum_mechanics/Umat/Finite/generic_hyper_invariants.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The elastic UMAT requires 2 constants:
///@brief props[0] : Young modulus
///@brief props[1] : Poisson ratio
///@brief props[2] : CTE

///@brief No statev is required for thermoelastic constitutive law

void umat_generic_hyper_invariants(const std::string &umat_name, const vec &etot, const vec &Detot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
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

    double dWdI_1_bar = 0.;
    double dWdI_2_bar = 0.;    

    double dW2dI_11_bar = 0.;    
    double dW2dI_12_bar = 0.;    
    double dW2dI_22_bar = 0.;    

    double dUdJ = 0.;
    double dU2dJ2 = 0.;

    double J = det(F1);  
    vec I_bar = isochoric_invariants(b, J);

    std::map<string, int> list_potentials;
    list_potentials = {{"NEOHC",0},{"MOORI",1},{"YEOHH",2}};

    switch (list_potentials[umat_name]) {
        case 0: {
            // \f$ W = \frac{\mu}{2}*left(\bar{I}_1 -3\right) + \kappa \left( J textrm{ln} J - J +1 \right) \f$ 
            double mu = props(0);
            double kappa = props(1);            
            dWdI_1_bar = 0.5*mu;
            dUdJ = kappa*log(J);
            dU2dJ2 = kappa/J;
            break;
        }
        case 1: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{01} left(\bar{I}_2 -3\right) + \kappa \left( J textrm{ln} J - J +1 \right) \f$ 
            double C_10 = props(0);
            double C_01 = props(1);            
            double kappa = props(2);                        
            dWdI_1_bar = C_10;
            dWdI_2_bar = C_01;            
            dUdJ = kappa*log(J);
            dU2dJ2 = kappa/J;
            break;
        }
        case 2: {
            // \f$ W = C_{10} left(\bar{I}_1 -3\right) + C_{20} left(\bar{I}_1 -3\right)^2 + C_{30} left(\bar{I}_1 -3\right)^3 + \kappa \left( J textrm{ln} J - J +1 \right) \f$             
            double C_10 = props(0);
            double C_20 = props(1);            
            double C_30 = props(2);            
            double kappa = props(3);     
            dWdI_1_bar = C_10 + 2.*C_20*(I_bar(0)-3.) + 3.*C_30*pow((I_bar(0)-3.),2.);
            dW2dI_11_bar = 2.*C_20 + 6.*C_30*(I_bar(0)-3.);    
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

    mat m_sigma_iso = sigma_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, b, J);
    mat m_sigma_vol = sigma_vol_hyper(dUdJ, b, J);    
    mat m_sigma = m_sigma_iso + m_sigma_vol;
    sigma = t2v_stress(m_sigma);  

    mat Lt_iso = L_iso_hyper_invariants(dWdI_1_bar, dWdI_2_bar, dW2dI_11_bar, dW2dI_12_bar, dW2dI_22_bar, b, J);    
    mat Lt_vol = L_vol_hyper(dUdJ, dU2dJ2, b, J);
    Lt = Lt_iso + Lt_vol;

    if(start) {
        L = Lt;
    }

    cout << "L = " << L << endl;
    cout << "Lt = " << Lt << endl;
    cout << "Lt_iso = " << Lt_iso << endl;    
    cout << "Lt_vol = " << Lt_vol << endl;        
    cout << "eig(Lt)" << eig_sym(Lt);
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_r += 0.5*sum((sigma_start+sigma)%Detot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
