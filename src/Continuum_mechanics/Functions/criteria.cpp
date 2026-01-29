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
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>

using namespace std;
using namespace arma;

namespace simcoon
{

    // This function returns the Drucker equivalent stress.
    double Drucker_stress(const vec &v, const double &b, const double &n)
    {
        assert(v.size() == 6);
        assert(b >= 0.);
        assert(n > 0.);
        double temp;
        double m;

        if (Mises_stress(v) > 0.)
        {
            if (n < 10.)
            {
                m = 1. / n;
                temp = Mises_stress(v) * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), m);
            }
            else
            {
                m = 0.;
                temp = Mises_stress(v);
            }
        }
        else
        {
            m = 0.;
            temp = 0.;
        }
        return temp;
    }

    // This function returns the derivative of the Drucker equivalent stress.
    vec dDrucker_stress(const vec &v, const double &b, const double &n)
    {
        assert(v.size() == 6);
        assert(b >= 0.);
        assert(n > 0.);  

        double Mises = Mises_stress(v);
        vec dMises = eta_stress(v);

        double J2 = J2_stress(v);
        double J3 = J3_stress(v);
        double J2_3_2 = pow(J2, 1.5);
        double J2_5_2 = pow(J2, 2.5);
        vec dJ3 = dJ3_stress(v);
        vec dJ2 = dJ2_stress(v);

        double m;
        vec temp;

        if (Mises > 0.)
        {
            if (n < 10.)
            {
                m = 1. / n;

                temp = pow((1. + b * J3 / J2_3_2), m) * dMises + Mises * b * m * (pow((1. + b * J3 / J2_3_2), (m-1))* (dJ3/J2_3_2-1.5*J3*dJ2/J2_5_2));
            }
            else
            {
                m = 0.;
                temp = eta_stress(v);
            }
        }
        else
        {
            m = 0.;
            temp = zeros(6);
        }
        return temp;
    }

    // This function returns the Drucker equivalent stress.
    double Tresca_stress(const vec &v)
    {
        mat sigma = v2t_stress(v);
        vec lambda = eig_sym(sigma);

        // eig_sym is in ascending order
        return lambda(2) - lambda(0);
    }

    // This function returns the Tresca equivalent stress.
    vec dTresca_stress(const vec &v)
    {
        return eta_stress(v);
    }

    // This function returns the derivative of J2
    vec dJ2_stress(const vec &v)
    {
        vec vdev = (dev(v)%Ir2());  

        return vdev;
    }

    // This function returns the derivative of J3
    vec dJ3_stress(const vec &v)
    {
        assert(v.size() == 6);

        vec vdev = dev(v);  
        mat S = v2t_stress(vdev);

        mat SS = S * S;

        double trS2 = accu(S % S);

        mat dJ3_mat = SS - (1.0/3.0) * trS2 * eye<mat>(3,3);

        vec dJ3 = t2v_strain(dJ3_mat);

        return dJ3;
    }
    // This function returns the combination of the Drucker equivalent stress by replacing VM by DFA
    double Drucker_ani_stress(const vec &v, const vec &params, const double &b, const double &n)
    {
        assert(v.size() == 6);
        assert(b >= 0.);
        assert(n > 0.);
        double temp;
        double m;

        double dfa_stress = DFA_stress(v, params);

        if (dfa_stress > 0.)
        {
            if (n < 10.)
            {
                m = 1. / n;
                temp = dfa_stress * pow((1. + b * J3_stress(v) / pow(J2_stress(v), 1.5)), m);
            }
            else
            {
                m = 0.;
                temp = dfa_stress;
            }
        }
        else
        {
            m = 0.;
            temp = 0.;
        }
        return temp;
    }

    // This function returns the derivative of the Drucker equivalent stress.
    vec dDrucker_ani_stress(const vec &v,  const vec &params, const double &b, const double &n)
    {
        assert(v.size() == 6);
        assert(b >= 0.);
        assert(n > 0.);  

        double Dfa_stress = DFA_stress(v, params);
        vec dDfa_stress = dDFA_stress(v,params);

        double J2 = J2_stress(v);
        double J3 = J3_stress(v);
        double J2_3_2 = pow(J2, 1.5);
        double J2_5_2 = pow(J2, 2.5);
        vec dJ3 = dJ3_stress(v);
        vec dJ2 = dJ2_stress(v);

        double m;
        vec temp;

        if (Dfa_stress > 0.)
        {
            if (n < 10.)
            {
                m = 1. / n;

                temp = pow((1. + b * J3 / J2_3_2), m) * dDfa_stress + Dfa_stress * b * m * (pow((1. + b * J3 / J2_3_2), (m-1))* (dJ3/J2_3_2-1.5*J3*dJ2/J2_5_2));
            }
            else
            {
                m = 0.;
                temp = eta_stress(v);
            }
        }
        else
        {
            m = 0.;
            temp = zeros(6);
        }
        return temp;
    }

    mat P_Ani(const vec &params)
    {
        assert(params.n_elem == 9); // P_11,P_22,P_33,P_12,P_13,P_23,P_44,P_55,P_66
        mat P = zeros(6, 6);
        P(0, 0) = params(0);      // P_11
        P(1, 1) = params(1);      // P_22
        P(2, 2) = params(2);      // P_33
        P(0, 1) = params(3);      // P_12
        P(1, 0) = params(3);      // P_12
        P(0, 2) = params(4);      // P_13
        P(2, 0) = params(4);      // P_13
        P(1, 2) = params(5);      // P_23
        P(2, 1) = params(5);      // P_23
        P(3, 3) = 2. * params(6); // P_44 = P_1212
        P(4, 4) = 2. * params(7); // P_55 = P_1313
        P(5, 5) = 2. * params(8); // P_66 = P_2323

        return P;
    }

    mat P_Hill(const vec &params)
    {
        assert(params.n_elem == 6); // F,G,H,L,M,N
        mat P = zeros(6, 6);
        // param(0) = F
        // param(1) = G
        // param(2) = H
        // param(3) = L
        // param(4) = M
        // param(5) = N
        P(0, 0) = params(1) + params(2); // P_11 = G+H
        P(1, 1) = params(0) + params(2); // P_22 = F+H
        P(2, 2) = params(0) + params(1); // P_33 = F+G
        P(0, 1) = -1. * params(2);       // P_12 = -H
        P(1, 0) = -1. * params(2);       // P_12 = -H
        P(0, 2) = -1. * params(1);       // P_13 = -G
        P(2, 0) = -1. * params(1);       // P_13 = -G
        P(1, 2) = -1. * params(0);       // P_23 = -F
        P(2, 1) = -1. * params(0);       // P_23 = -F
        P(3, 3) = 2. * params(5);        // P_44 = N
        P(4, 4) = 2. * params(4);        // P_55 = M
        P(5, 5) = 2. * params(3);        // P_66 = L

        return P;
    }

    mat P_DFA(const vec &params)
    {                               // Deshpande–Fleck–Ashby
        assert(params.n_elem == 7); // F,G,H,L,M,N,K
        mat P = zeros(6, 6);
        // param(0) = F
        // param(1) = G
        // param(2) = H
        // param(3) = L
        // param(4) = M
        // param(5) = N
        // param(6) = K
        P(0, 0) = params(1) + params(2) + params(6) / 9.; // P_11 = G+H+K/9
        P(1, 1) = params(0) + params(2) + params(6) / 9.; // P_22 = F+H+K/9
        P(2, 2) = params(0) + params(1) + params(6) / 9.; // P_33 = F+G+K/9
        P(0, 1) = -1. * params(2) + params(6) / 9.;       // P_12 = -H+K/9
        P(1, 0) = -1. * params(2) + params(6) / 9.;       // P_12 = -H+K/9
        P(0, 2) = -1. * params(1) + params(6) / 9.;       // P_13 = -G+K/9
        P(2, 0) = -1. * params(1) + params(6) / 9.;       // P_13 = -G+K/9
        P(1, 2) = -1. * params(0) + params(6) / 9.;       // P_23 = -F+K/9
        P(2, 1) = -1. * params(0) + params(6) / 9.;       // P_23 = -F+K/9
        P(3, 3) = 2. * params(5);                         // P_44 = N
        P(4, 4) = 2. * params(4);                         // P_55 = M
        P(5, 5) = 2. * params(3);                         // P_66 = L

        return P;
    }

    double Eq_stress_P(const vec &v, const mat &H)
    {

        if (norm(v, 2) > iota)
        {
            return pow(sum(v % (H * v)), 0.5);
        }
        else
        {
            return 0.;
        }
    }

    vec dEq_stress_P(const vec &v, const mat &H)
    {

        if (norm(v, 2) > iota)
        {
            return (H * v) / Eq_stress_P(v, H);
        }
        else
        {
            return zeros(6);
        }
    }

    double Hill_stress(const vec &v, const vec &params)
    {
        mat P = P_Hill(params);
        return Eq_stress_P(v, P);
    }

    vec dHill_stress(const vec &v, const vec &params)
    {
        mat P = P_Hill(params);
        return dEq_stress_P(v, P);
    }

    double DFA_stress(const vec &v, const vec &params)
    {
        mat P = P_DFA(params);
        return Eq_stress_P(v, P);
    }

    vec dDFA_stress(const vec &v, const vec &params)
    {
        mat P = P_DFA(params);
        return dEq_stress_P(v, P);
    }

    double Ani_stress(const vec &v, const vec &params)
    {
        mat P = P_Ani(params);
        return Eq_stress_P(v, P);
    }

    vec dAni_stress(const vec &v, const vec &params)
    {
        mat P = P_Ani(params);
        return dEq_stress_P(v, P);
    }

    double Eq_stress(const vec &v, const string &eq_type, const vec &param)
    {
        if (eq_type == "Mises")
        {
            return Mises_stress(v);
        }
        else if (eq_type == "Tresca")
        {
            return Tresca_stress(v);
        }
        else if (eq_type == "Drucker")
        {
            return Drucker_stress(v, param(0), param(1));
        }
        else if (eq_type == "Hill")
        {
            return Hill_stress(v, param);
        }
        else if (eq_type == "Ani")
        {
            return Ani_stress(v, param);
        }
        else
        {
            cout << "Error in Eq_stress : No valid arguement is given\n";
            exit(0);
        }
    }

    vec dEq_stress(const vec &v, const string &eq_type, const vec &param)
    {
        if (eq_type == "Mises")
        {
            return eta_stress(v);
        }
        else if (eq_type == "Tresca")
        {
            return dTresca_stress(v);
        }
        else if (eq_type == "Drucker")
        {
            return dDrucker_stress(v, param(0), param(1));
        }
        else if (eq_type == "Hill")
        {
            return dHill_stress(v, param);
        }
        else if (eq_type == "Ani")
        {
            return dAni_stress(v, param);
        }
        else
        {
            cout << "Error in dEq_stress : No valid arguement is given\n";
            exit(0);
        }
    }

} // namespace simcoon
