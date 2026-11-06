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

///@file output.hpp
///@brief Function to output phases data

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file output.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */

    
/**
 * @brief Copy a matrix into the state variables, row by row.
 * @param L matrix to store (n_rows x n_cols)
 * @param i0 index of statev where the first entry is written
 * @param statev state variables, entries i0 .. i0 + n_rows*n_cols - 1 overwritten
 */
void output_modulus(const arma::mat &L, const int &i0, arma::vec &statev);

/**
 * @brief Copy a list of matrices into the state variables, one after the other, each row by row.
 * @param list matrices to store
 * @param i0 index of statev where the first entry of the first matrix is written
 * @param statev state variables, overwritten from i0 on by the total number of entries of the list
 */
void output_moduli(const std::vector<arma::mat> &list, const int &i0, arma::vec &statev);


/** @} */ // end of phase group

} //namespace simcoon
