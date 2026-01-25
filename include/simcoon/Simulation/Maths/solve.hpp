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
#pragma once

#include <armadillo>

namespace simcoon{

/**
 * @file solve.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */


arma::vec quadratic(const double &, const double &, const double &);

arma::cx_vec cx_quadratic(const double &, const double &, const double &);

arma::cx_vec cx_quadratic(const arma::cx_double &, const arma::cx_double &, const arma::cx_double &);


/** @} */ // end of maths group

} //namespace simcoon    
