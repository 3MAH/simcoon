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

///@file stress.hpp
///@brief A set of functions that transforms stress measures into anothers (considering Finite strains)
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

//This function returns the first Piola-Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Cauchy2PKI(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the second Piola-Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Cauchy2PKII(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Cauchy2Kirchoff(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Kirchoff stress tensor from the Cauchy stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::vec Cauchy2Kirchoff(const arma::vec &, const arma::mat &, const double & = 0.);

//This function returns the Cauchy stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Kirchoff2Cauchy(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Cauchy stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::vec Kirchoff2Cauchy(const arma::vec &, const arma::mat &, const double & = 0.);

//This function returns the second Piola-Kirchoff stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Kirchoff2PKI(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the second Piola-Kirchoff stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat Kirchoff2PKII(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the second Piola-Kirchoff stress tensor from the Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::vec Kirchoff2PKII(const arma::vec &, const arma::mat &, const double & = 0.);

//This function returns the Kirchoff stress tensor from the second Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat PKI2Kirchoff(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Kirchoff stress tensor from the first Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat PKII2Kirchoff(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Cauchy stress tensor from the first Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat PKI2Cauchy(const arma::mat &, const arma::mat &, const double & = 0.);

//This function returns the Cauchy stress tensor from the second Piola-Kirchoff stress tensor, the transformation gradient F and its determinant (optional, if not indicated, it will be computed)
arma::mat PKII2Cauchy(const arma::mat &, const arma::mat &, const double & = 0.);
    
} //namespace simcoon
