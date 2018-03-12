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

///@file stats.cpp
///@brief Usefull statistical functions
///@version 1.0

#include <math.h>
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_Mechanics/Functions/contimech.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
///Approximation of a normal distribution
double normal_distrib(const double &x, const double &mean, const double &dev){
	
	double x_norm = (x-mean)/dev;
	
	if (x_norm < 0.){
		return (1. - normal_distrib(-1.*x_norm, 0., 1.));
	}
	else {
		double k = 1.0/(1.0 + 0.2316419*x_norm);
		double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
		return (1.0 - (1.0/(sqrt(2*pi)))*exp(-0.5*x_norm*x_norm) * k_sum);
	}	
}

double proba_distrib_weibull(const double &x, const double &alpha, const double &beta) {
    //This is the implementation of a Wibull cumulative distribution function, where
    // alpha > 0 is the shape parameter and λ > 0 is the scale parameter of the distribution.
    // Its complementary cumulative distribution function is a stretched exponential function. The Weibull distribution is related to a number of other probability distributions; in particular, it interpolates between the exponential distribution (alpha = 1) and the Rayleigh distribution (alpha = 2 and \beta = \sqrt{2}\sigma)
    
    assert(alpha > 0.);
    assert(beta > 0.);
    
    return Macaulay_p((alpha/beta)*pow((x/beta),alpha-1.)*exp(-1.* pow(x/beta, alpha)));
}
    
double cumul_distrib_weibull(const double &x, const double &alpha, const double &beta) {
    //This is the implementation of a Wibull cumulative distribution function, where
    // alpha > 0 is the shape parameter and λ > 0 is the scale parameter of the distribution.
    // Its complementary cumulative distribution function is a stretched exponential function. The Weibull distribution is related to a number of other probability distributions; in particular, it interpolates between the exponential distribution (alpha = 1) and the Rayleigh distribution (alpha = 2 and \beta = \sqrt{2}\sigma)
    
    assert(alpha > 0.);
    assert(beta > 0.);
    
    return Macaulay_p(1. - exp(-1.* pow(x/beta, alpha)));
}

//tri_sum of a and b
int tri_sum(const int &a, const int &b) {
    assert(a>0);
    assert(b>0);
    
    return b*(2*a+(b-1))/2;
}

///1. Classic ODF: a1 * cos(Theta)^(2*p1) + a2 * cos(Theta)^(2*p2 + 1) * sin(Theta)^(2*p2)
    
double ODF_sd(const double &Theta, const double &mean, const vec &params) {
	  
    double alpha1 = params(0);
    double alpha2 = params(1);
    double pow1 = params(2);
    double pow2 = params(3);
    
	if (fabs(Theta - mean) < 1.E-6)
		return alpha1;  
	else if (fabs((Theta - mean)-0.5*pi) < 1.E-6)
		return 0.;
	else
		return fabs((alpha1*pow(cos(Theta - mean),2.*pow1) + alpha2*pow(cos(Theta - mean),2.*pow2)*pow(sin(Theta - mean),2.*pow2))*cos(Theta - mean));
}

///2. Classic ODF - hardening-like
double ODF_hard(const double &Theta, const double &mean, const double &std_dev, const double &ampl) {

    assert(ampl>0);
    assert(std_dev>0);
	  
	return ampl*exp(-0.5*pow(fabs(Theta - mean)/std_dev,2.));
}

///3. Gaussian
double Gaussian(const double &X, const double &mean, const double &std_dev, const double &ampl){

    assert(std_dev>0);
	
	return ampl/(std_dev * sqrt(2.*pi)) * exp(-1./2. * pow((X - mean)/std_dev, 2.));
}

///4. Lorentzian
double Lorentzian(const double &X, const double &mean, const double &width, const double &ampl){

    assert(width>0);
    
	return ampl * width/(2.*pi*(pow(X - mean, 2.)+ pow(width / 2., 2.)));
}

///5. Pseudo-Voigt
double PseudoVoigt(const double &X, const double &mean, const double &std_dev, const double &width_lor, const double &ampl, const vec &params){
    
    double eta = params(0);
    
    assert(width_lor>0);
    assert(std_dev>0);
	  
	return eta * Lorentzian(X, mean, width_lor, ampl) + (1.-eta) * Gaussian(X, mean, std_dev, ampl);
}

///6. Pearson VII
double Pearson7(const double& X, const double& mean, const double &inv_width, vec &params){
	
    double max = params(0);
    double shape = params(1);
    assert(shape>0);
    if (fabs(max) < limit) {
        max = 1.;
    }
    
	return max * pow(1. + pow(inv_width*(X - mean), 2.)/shape, -1.*shape);
}

} //namespace simcoon
