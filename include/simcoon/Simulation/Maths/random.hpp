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

///@file random.hpp
///@brief random number generators
///@version 1.0

#pragma once

namespace simcoon{

/**
 * @file random.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */

/**
 * @brief Generates a random integer between 0 and a-1.
 * @param a The upper bound (exclusive) (int)
 * @return A random integer in the range [0, a-1] (int)
 * @details Example:
 * @code
 *      int r = alea(10);  // Returns a random integer between 0 and 9
 * @endcode
 */
int alea(const int &a);

/**
 * @brief Generates a random integer between a and b.
 * @param a The lower bound (inclusive) (int)
 * @param b The upper bound (inclusive) (int)
 * @return A random integer in the range [a, b] (int)
 * @details Example:
 * @code
 *      int r = aleab(5, 15);  // Returns a random integer between 5 and 15
 * @endcode
 */
int aleab(const int &a, const int &b);

/**
 * @brief Generates a random double between a and b.
 * @param a The lower bound (double)
 * @param b The upper bound (double)
 * @return A random double in the range [a, b] (double)
 * @details Example:
 * @code
 *      double r = alead(0.0, 1.0);  // Returns a random double between 0.0 and 1.0
 * @endcode
 */
double alead(const double &a, const double &b);


/** @} */ // end of maths group

} //namespace simcoon