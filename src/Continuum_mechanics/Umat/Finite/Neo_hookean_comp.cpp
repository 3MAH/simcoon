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

void umat_neo_hookean_comp(const vec &Etot, const vec &DEtot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, vec &sigma_in, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, const int &solver_type, double &tnew_dt)
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

    double mu = E/(2.*(1+nu));
    double lambda = E*nu/((1.+nu)*(1.-2.*nu));
    
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
    double J = det(F1); //lambda(2)*lambda(1)*lambda(0)
    
    double W = mu/2.0*(I1-3.) - mu*log(J) + lambda/2.0*pow(log(J),2.);

    mat invC = inv(C);
    mat I = eye(3,3);

    //Compute the PKII stress and then the Cauchy stress
    mat S = mu*(I-invC) + lambda*log(J)*invC;

    cout << "J = " << J;
    cout << "invC = " << invC;
    cout << "S = " << S;
    
    sigma = t2v_stress(PKII2Cauchy(S, F1));
    vec sigma2 = t2v_stress(mu/J*(L_Cauchy_Green(F1) - I) + lambda*log(J)/J*I);
  
    cout << "sigma = " << sigma.t();
    cout << "sigma2 = " << sigma2.t();
    
//    L = lambda*sym_dyadic(invC,invC)+2.0*(mu-lambda*log(J))*dinvSdS(C);
    L = lambda*auto_dyadic(invC)+2.0*(mu-lambda*log(J))*dinvSdS(C);
    Lt = DSDE_2_DtauDe(L, get_BBBB(F1), F1, v2t_stress(sigma));
//    Lt = (1./J)*(DSDE_2_Dtau_LieDD(L, F1));
    
    
    double lambdap = lambda/J;
    double mup =(mu-lambda*log(J))/J;
    
//    mat Lt2 = lambdap*sym_dyadic(I,I)+2.0*mup*dinvSdS(I);
    mat Lt2 = lambdap*auto_dyadic(I)+2.0*mup*dinvSdS(I);
    
    cout << "L = " << L << endl;
    cout << "auto_dyadic(invC,invC)" << auto_dyadic(invC) << endl;
    cout << "dinvSdS(C)" << dinvSdS(C) << endl;
    cout << "Lt = " << Lt << endl;
    cout << "Lt2 = " << Lt2 << endl;
    cout << "Lt_iso = " << L_iso(mup, lambdap, "mulambda") << endl;
    
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
