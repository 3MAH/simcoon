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

///@file solve.hpp
///@brief random number generators
///@version 1.0

#include <iostream>
#include <assert.h>
#include <cmath>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

namespace simcoon{

vec quadratic(const double &a, const double &b, const double &c)
{
    assert(fabs(a)+fabs(b) > 0.);
    double delta;
    vec t = zeros(2);
    if (a == 0.)
    {
        t.resize(1);
        t(0) = -c/b;
    }
    else
    {
        delta=pow(b,2.)-4.*a*c;
        assert(delta > 0.);
        
        t(0) = (-b+sqrt(delta))/(2.*a);
        t(1) = (-b-sqrt(delta))/(2.*a);
    }
    
    return t;
}


cx_vec cx_quadratic(const double &a, const double &b, const double &c)
{
    cx_double cplx_a = cx_double (a);
    cx_double cplx_b = cx_double (b);
    cx_double cplx_c = cx_double (c);
    
    cx_double delta;
    cx_vec t = zeros<cx_vec>(2);
    delta=pow(cplx_b,2.)-4.*cplx_a*cplx_c;
    if (a == 0.)
    {
        t.resize(1);
        t(0) = -cplx_c/cplx_b;
    }
    else if ( delta == complex<double>(0) )
    {
        t.resize(1);
        t(0) = -cplx_b/(2.*cplx_a);
	}
    else
    {
        t(0) = (-cplx_b+sqrt(delta))/(2.*cplx_a);
        t(1) = (-cplx_b-sqrt(delta))/(2.*cplx_a);
    }
    return t;
}


cx_vec cx_quadratic(const cx_double &cplx_a, const cx_double &cplx_b, const cx_double &cplx_c)
{
    cx_double delta;
    cx_vec t = zeros<cx_vec>(2);
    delta=pow(cplx_b,2.)-4.*cplx_a*cplx_c;
    if (cplx_a == 0.)
    {
        t.resize(1);
        t(0) = -cplx_c/cplx_b;
    }
    else if ( delta == complex<double>(0) )
    {
        t.resize(1);
        t(0) = -cplx_b/(2.*cplx_a);
	}
    else
    {
        t(0) = (-cplx_b+sqrt(delta))/(2.*cplx_a);
        t(1) = (-cplx_b-sqrt(delta))/(2.*cplx_a);
    }
    return t;
}

} //namespace simcoon    

