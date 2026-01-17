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

///@file doe.hpp
///@brief Design of Experiments library
///@version 1.0

#include <iostream>
#include <armadillo>
#include "parameters.hpp"
#include "generation.hpp"

namespace simcoon{

/**
 * @file doe.hpp
 * @brief Parameter identification functions.
 */

/** @addtogroup identification
 *  @{
 */


//This function computes the test matrix with the parameters of a uniform multidimensional distribution
arma::mat doe_uniform(const int &, const int &, const std::vector<parameters> &);

//This function computes the test matrix with the parameters of a uniform multidimensional distribution, where the extreme samples are in the bounds
arma::mat doe_uniform_limit(const int &, const int &, const std::vector<parameters> &);
    
//This function computes the test matrix with the parameters of a random sampling
arma::mat doe_random(const int &, const int &, const std::vector<parameters> &);

//This function is utilized to initialize the first generation
void gen_initialize(generation &, int &, int&, int &, const int &, const int &, const std::vector<parameters> &, const double &);
    

/** @} */ // end of identification group

} //namespace simcoon
