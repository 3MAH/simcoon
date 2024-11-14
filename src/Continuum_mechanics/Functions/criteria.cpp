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

///@file criteria.cpp
///@brief Provide function for yield surfaces
///@version 1.0

#include <iostream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//This function returns the Prager equivalent stress.
double Prager_stress(const vec &v, const double &b, const double &n)
{
     assert(v.size() == 6);
     assert(b >= 0.);
     assert(n > 0.);
     double temp;
     double m;
    
     if (Mises_stress(v) > 0.) {
         if (n < 10.)
         {
             m = 1. / n;
             temp = Mises_stress(v) * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), m);
         }
         else {
             m = 0.;
             temp = Mises_stress(v);
         }
     }
     else {
         m = 0.;
         temp = 0.;
     }
     return temp;
}

//This function returns the derivative of the Prager equivalent stress.
vec dPrager_stress(const vec &v, const double &b, const double &n)
{
     assert(v.size() == 6);
     assert(b >= 0.);
     assert(n > 0.);
     vec vdev = dev(v);
     mat devstress_t = v2t_stress(vdev);
     vec square_stressdev = t2v_stress(devstress_t * devstress_t);
     double m;
     
     vec temp;
    
     if (Mises_stress(v) > 0.)
     {
         if (n < 10.) {
             m = 1. / n;
             temp = sqrt(3.) * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), (m - 1.)) * (0.5 / sqrt(J2_stress(v)) * vdev + b * m / (6. * pow(J2_stress(v), 2.)) * (6. * J2_stress(v) * square_stressdev - 4. * pow(J2_stress(v), 2.) * Ith() + (3. / m - 9.) * J3_stress(v) * vdev));
             
             for (int i = 3; i < 6; i++)
             {
             temp(i) = 2. * temp(i);
             }
         }
         else {
             m = 0.;
             temp = eta_stress(v);
         }         
     }
     else {
         m = 0.;
         temp = zeros(6);
         
     }
     return temp;

}

//This function returns the Prager equivalent stress.
double Tresca_stress(const vec &v)
{
    mat sigma = v2t_stress(v);
    vec lambda = eig_sym(sigma);

    //eig_sym is in ascending order
    return lambda(2) - lambda(0);
}

//This function returns the Tresca equivalent stress.
vec dTresca_stress(const vec &v)
{
    return eta_stress(v);
}

mat P_Ani(const vec &params) {
    assert(params.n_elem == 9); //P_11,P_22,P_33,P_12,P_13,P_23,P_44,P_55,P_66
    mat P = zeros(6,6);
    P(0,0) = params(0); //P_11
    P(1,1) = params(1); //P_22
    P(2,2) = params(2); //P_33
    P(0,1) = params(3); //P_12
    P(1,0) = params(3); //P_12
    P(0,2) = params(4); //P_13
    P(2,0) = params(4); //P_13
    P(1,2) = params(5); //P_23
    P(2,1) = params(5); //P_23
    P(3,3) = 2.*params(6); //P_44 = P_1212
    P(4,4) = 2.*params(7); //P_55 = P_1313
    P(5,5) = 2.*params(8); //P_66 = P_2323
    
    return P;
}
    
mat P_Hill(const vec &params) {
    assert(params.n_elem == 6); //F,G,H,L,M,N
    mat P = zeros(6,6);
    //param(0) = F
    //param(1) = G
    //param(2) = H
    //param(3) = L
    //param(4) = M
    //param(5) = N
    P(0,0) = params(1) + params(2); //P_11 = G+H
    P(1,1) = params(0) + params(2); //P_22 = F+H
    P(2,2) = params(0) + params(1); //P_33 = F+G
    P(0,1) = -1.*params(2); //P_12 = -H
    P(1,0) = -1.*params(2); //P_12 = -H
    P(0,2) = -1.*params(1); //P_13 = -G
    P(2,0) = -1.*params(1); //P_13 = -G
    P(1,2) = -1.*params(0); //P_23 = -F
    P(2,1) = -1.*params(0); //P_23 = -F
    P(3,3) = 2.*params(5); //P_44 = N
    P(4,4) = 2.*params(4); //P_55 = M
    P(5,5) = 2.*params(3); //P_66 = L
    
    return P;
}

mat P_DFA(const vec &params) //Deshpande–Fleck–Ashby {
    assert(params.n_elem == 7); //F,G,H,L,M,N,K
    mat P = zeros(6,6);
    //param(0) = F
    //param(1) = G
    //param(2) = H
    //param(3) = L
    //param(4) = M
    //param(5) = N
    P(0,0) = params(1) + params(2) + params(6); //P_11 = G+H+K
    P(1,1) = params(0) + params(2) + params(6); //P_22 = F+H+K
    P(2,2) = params(0) + params(1) + params(6); //P_33 = F+G+K
    P(0,1) = -1.*params(2) + params(6); //P_12 = -H+K
    P(1,0) = -1.*params(2) + params(6); //P_12 = -H+K
    P(0,2) = -1.*params(1) + params(6); //P_13 = -G+K
    P(2,0) = -1.*params(1) + params(6); //P_13 = -G+K
    P(1,2) = -1.*params(0) + params(6); //P_23 = -F+K
    P(2,1) = -1.*params(0) + params(6); //P_23 = -F+K
    P(3,3) = 2.*params(5); //P_44 = N
    P(4,4) = 2.*params(4); //P_55 = M
    P(5,5) = 2.*params(3); //P_66 = L
    
    return P;
}

double Eq_stress_P(const vec &v, const mat &H) {
    
    if (norm(v,2) > sim_iota) {
        return pow(sum(v%(H*v)),0.5);
    }
    else {
        return 0.;
    }
}
                   
vec dEq_stress_P(const vec &v, const mat &H) {
   
   if (norm(v,2) > sim_iota) {
       return (H*v)/Eq_stress_P(v,H);
   }
   else {
       return zeros(6);
   }
}
    
double Hill_stress(const vec &v, const vec &params) {
    mat P = P_Hill(params);
    return Eq_stress_P(v,P);
}

vec dHill_stress(const vec &v, const vec &params) {
   mat P = P_Hill(params);
   return dEq_stress_P(v,P);
}

double DFA_stress(const vec &v, const vec &params) {
    mat P = P_DFA(params);
    return Eq_stress_P(v,P);
}

vec dDFA_stress(const vec &v, const vec &params) {
   mat P = P_DFA(params);
   return dEq_stress_P(v,P);
}
                   
double Ani_stress(const vec &v, const vec &params) {
    mat P = P_Ani(params);
    return Eq_stress_P(v,P);
}

vec dAni_stress(const vec &v, const vec &params) {
    mat P = P_Ani(params);
    return dEq_stress_P(v,P);
}
                   
double Eq_stress(const vec &v, const string &eq_type, const vec &param)
{
    if(eq_type == "Mises") {
        return Mises_stress(v);
    }
    else if(eq_type == "Tresca") {
        return Tresca_stress(v);
    }
    else if(eq_type == "Prager") {
        return Prager_stress(v, param(0), param(1));
    }
    else if(eq_type == "Hill") {
        return Hill_stress(v, param);
    }
    else if(eq_type == "Ani") {
        return Ani_stress(v, param);
    }
    else {
        cout << "Error in Eq_stress : No valid arguement is given\n";
        exit(0);
    }
}

vec dEq_stress(const vec &v, const string &eq_type, const vec &param)
{
    if(eq_type == "Mises") {
        return eta_stress(v);
    }
    else if(eq_type == "Tresca") {
        return dTresca_stress(v);
    }
    else if(eq_type == "Prager") {
        return dPrager_stress(v, param(0), param(1));
    }
    else if(eq_type == "Hill") {
        return dHill_stress(v, param);
    }
    else if(eq_type == "Ani") {
        return dAni_stress(v, param);
    }
    else {
        cout << "Error in dEq_stress : No valid arguement is given\n";
        exit(0);
    }
}
    
} //namespace simcoon
