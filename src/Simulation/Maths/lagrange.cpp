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

///@file lagrange.cpp
///@brief Function that are utilized in constrained problems
///@version 1.0

#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

namespace simcoon{

double lagrange_exp(const double &h, const double &c, const double &p0) {
	if(h <= -c)
	{
		return 0.;
	}
	else
	{
		return (p0/(exp(1.) - 1.)) * ((h/c + 1.)*(exp(h/c + 1.) - 1.));
	}
}

double dlagrange_exp(const double &h, const double &c, const double &p0) {
	if(h <= -c)
	{
		return 0.;
	}
	else
	{
		return (p0/(exp(1.) - 1.)) * ((1./c)*(h/c +2.)*exp(h/c + 1.) - (1./c));
	}
    
}

double lagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha) {
    
	double h = c;
    
	if(x >= h)
	{
		return p0*(1.-x)*pow(x,-n);
	}
	else
	{
		return p0*((-1.*pow(h,-n) - n*(1.-h)*pow(h,-n-1.))*(x-h) + (1.-h)*pow(h,-n)) + alpha*pow(x-h,2.);
	}
}

double dlagrange_pow_0(const double &x, const double &c, const double &p0, const double &n, const double &alpha)
{
    
	double h = c;
    
	if(x >= h)
	{
		return p0*(-1.*pow(x,-n) - n*(1.-x)*pow(x,-n-1.));
	}
	else
	{
		return p0*(-1.*pow(h,-n) - n*(1.-h)*pow(h,-n-1.)) + 2.*alpha*(x-h);
	}
}

double lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha)
{
	double h = 1.-c;
	
	if(x <= h)
	{
		return p0*x*pow(1.-x,-n);
	}
	else
	{
		return p0*((pow(1.-h,-n) + n*h*pow(1.-h,-n-1.))*(x-h) + h*pow(1.-h,-n)) + alpha*pow(x-h,2.);
	}
}

double dlagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha)
{
	double h = 1.-c;
	
	if(x <= h)
	{
		return p0*(pow(1.-x,-n) + n*x*pow(1.-x,-n-1.));
	}
	else
	{
		return p0*(pow(1.-h,-n) + n*h*pow(1.-h,-n-1.)) + 2.*alpha*(x-h);
	}
}

double d2lagrange_pow_1(const double &x, const double &c, const double &p0, const double &n, const double &alpha)
{
	double h = 1.-c;
	
	if(x <= h)
	{
		return p0*(2.*n*pow(1. - x, -n - 1.) + n*(n + 1.)*x*pow(1. - x, -n - 2.));
	}
	else
	{
		return 2. * alpha;
	}			
} 
}//namespace simcoon
