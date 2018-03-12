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

///@file stats.hpp
///@brief Usefull statistical functions
///@version 1.0

#pragma once

#include <armadillo>

namespace simcoon{

///Approximation of a normal distribution
double normal_distrib(const double &, const double &, const double &);

//Compute the probabilistic value for a given x, providing the shape parameter alpha and the scale parameter beta
double proba_distrib_weibull(const double &, const double &, const double &);

//Compute the probabilistic value for a given x, providing the shape parameter alpha and the scale parameter beta
double cumul_distrib_weibull(const double &, const double &, const double &);
    
//tri_sum of a and b
int tri_sum(const int &, const int &);

///1. Classic ODF: a1 * cos(Theta)^(2*p1) + a2 * cos(Theta)^(2*p2 + 1) * sin(Theta)^(2*p2) 
double ODF_sd(const double &, const double &, const arma::vec &);

///2. Classic ODF - hardening-like
double ODF_hard(const double &, const double &, const double &, const double &);

///3. Gaussian
double Gaussian(const double &, const double &, const double &, const double & = 1.);

///4. Lorentzian
double Lorentzian(const double &, const double &, const double &, const double & = 1.);

///5. Pseudo-Voigt
double PseudoVoigt(const double &, const double &, const double &, const double &, const double &, const arma::vec & = arma::ones(1));

///6. Pearson VII
double Pearson7(const double &, const double &, const double &, arma::vec &);

} //namespace simcoon
