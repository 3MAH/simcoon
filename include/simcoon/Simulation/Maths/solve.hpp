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
///@brief Roots of quadratic equations (real and complex)
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


/**
 * @brief Real roots of \f$ a t^2 + b t + c = 0 \f$.
 * @param a quadratic coefficient (if 0, the linear root \f$ -c/b \f$ is returned)
 * @param b linear coefficient
 * @param c constant coefficient
 * @return the two roots \f$ (-b \pm \sqrt{\Delta})/(2a) \f$, \f$ \Delta = b^2 - 4ac \f$ (asserted \f$ > 0 \f$), or one root when \f$ a = 0 \f$
 */
arma::vec quadratic(const double &a, const double &b, const double &c);

/**
 * @brief Complex roots of \f$ a t^2 + b t + c = 0 \f$ with real coefficients.
 * @param a quadratic coefficient (if 0, the linear root \f$ -c/b \f$ is returned)
 * @param b linear coefficient
 * @param c constant coefficient
 * @return the two roots \f$ (-b \pm \sqrt{\Delta})/(2a) \f$, or a single root when \f$ a = 0 \f$ or \f$ \Delta = 0 \f$
 */
arma::cx_vec cx_quadratic(const double &a, const double &b, const double &c);

/**
 * @brief Complex roots of \f$ a t^2 + b t + c = 0 \f$ with complex coefficients.
 * @param a quadratic coefficient (if 0, the linear root \f$ -c/b \f$ is returned)
 * @param b linear coefficient
 * @param c constant coefficient
 * @return the two roots \f$ (-b \pm \sqrt{\Delta})/(2a) \f$, or a single root when \f$ a = 0 \f$ or \f$ \Delta = 0 \f$
 */
arma::cx_vec cx_quadratic(const arma::cx_double &a, const arma::cx_double &b, const arma::cx_double &c);


/** @} */ // end of maths group

} //namespace simcoon    
