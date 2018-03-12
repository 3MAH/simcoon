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

///@file criteria.hpp
///@brief Provide function for yield surfaces
///@version 1.0

#pragma once
#include <string.h>
#include <armadillo>

namespace simcoon{
    
//This function returns the Prager equivalent stress.
double Prager_stress(const arma::vec &, const double &, const double &);

//This function returns the derivative of the Prager equivalent stress.
arma::vec dPrager_stress(const arma::vec &, const double &, const double &);

//This function returns the Tresca equivalent stress.
double Tresca_stress(const arma::vec &);

//This function returns the the derivative of the Tresca equivalent stress.
arma::vec dTresca_stress(const arma::vec &);
    
//This function returns the Full anisotropic tensor from components
arma::mat P_ani(const arma::vec &);

//This function returns the Hill anisotropic tensor, providing F,G,H,L,M,N
arma::mat P_hill(const arma::vec &);

//This function returns the anisotropic equivalent stress, given a matrix H
double Ani_stress(const arma::vec &, const arma::mat &);

//This function returns the derivative of the anisotropic equivalent stress, given a matrix H
arma::vec dAni_stress(const arma::vec &, const arma::mat &H);
    
//This function returns the anisotropic equivalent stress, providing the anisotropic tensor H components
double Ani_stress(const arma::vec &, const arma::vec &);

//This function returns the derivative of the anisotropic equivalent stress, providing the anisotropic tensor H components
arma::vec dAni_stress(const arma::vec &, const arma::vec &);
    
//This function returns the anisotropic equivalent stress, providing F,G,H,L,M,N
double Hill_stress(const arma::vec &v, const arma::vec &);

//This function returns the derivative of the anisotropic equivalent stress, providing F,G,H,L,M,N
arma::vec dHill_stress(const arma::vec &, const arma::vec &);
    
//This function computes the selected equivalent stress function
double Eq_stress(const arma::vec &, const std::string &, const arma::vec &);

//This function computes the derivative of the selected equivalent stress function
arma::vec dEq_stress(const arma::vec &, const std::string &, const arma::vec &);
    
} //namespace simcoon
