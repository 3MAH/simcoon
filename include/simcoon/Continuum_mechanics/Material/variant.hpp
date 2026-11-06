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

///@file variant.hpp
///@brief Class that defines characteristics of a variant of martensite
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file variant.hpp
 * @brief Martensite variant characteristics.
 */

/** @addtogroup material
 *  @{
 */


//======================================
class variant
//======================================
{
	private:

	protected:

	public :
		
		arma::vec n;
		arma::vec m;
		arma::mat R;
		double g;
		arma::vec ETn;	//Strain associated to each variant
		
		/** @brief Default constructor: zero normal, direction, Schmid tensor and transformation strain. */
		variant();
		/**
		 * @brief Constructor with parameters; builds the Schmid tensor and the transformation strain.
		 * @param n habit plane normal (size 3)
		 * @param m transformation direction (size 3)
		 * @param g transformation magnitude, \f$ \varepsilon^{T}_n = g\,\mathrm{Schmid\_v}(n, m) \f$
		 */
		variant(arma::vec n, arma::vec m, double g);
		/** @brief Copy constructor. */
		variant(const variant &);
		/** @brief Destructor. */
		~variant();
		
		/**
		 * @brief Rebuild the Schmid tensor R and the transformation strain ETn from the stored n, m.
		 * @param g transformation magnitude (must be > 0)
		 */
		void build(const double &g);
		/** @brief Copy assignment. */
		virtual variant& operator = (const variant&);
		
        friend std::ostream& operator << (std::ostream&, const variant&);
};


/** @} */ // end of material group

} //namespace simcoon
