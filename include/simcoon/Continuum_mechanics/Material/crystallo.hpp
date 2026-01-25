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

///@file crystallo.hpp
///@brief Some definitions coming from crystallography
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file crystallo.hpp
 * @brief Crystallography functions.
 */

/** @addtogroup material
 *  @{
 */


//This function returns the Shmidt tensor (3x3 matrix)
arma::mat Schmid(const arma::vec &n, const arma::vec &m);

//This function returns the Shmidt tensor (6 vector), with the convention of strain
arma::vec Schmid_v(const arma::vec &n, const arma::vec &m);

//This function returns a matrix utilized for the Hill interfacial operator
arma::mat F_nm(const arma::vec &N);

//This function returns the Hill interfacial operator for an isotropic material
arma::mat Q_nm(const arma::vec &N, const double &mu, const double &lambda);


/** @} */ // end of material group

} //namespace simcoon
