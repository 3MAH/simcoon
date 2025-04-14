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

/**
* @file natural_basis.hpp
* @author Yves Chemisky 
* @section The natural_basis library contains a set of function that help to define a natural basis in curvilinear coordinates
*/

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class natural_basis
//======================================
{
	private:

	protected:

	public :
	    /**
	     * @brief covariant vectors
	     * 
	     */
        std::vector<arma::vec> g_i; 

	    /**
	     * @brief contravariant vectors
	     * 
	     */
        std::vector<arma::vec> g0i;
    
	    /**
	     * @brief covariant components of the metric tensor
	     * 
	     */
        arma::mat g_ij;

	    /**
	     * @brief contravariant components of the metric tensor
	     * 
	     */
        arma::mat g0ij; // contravariant components of the metric tensor

		/**
		 * @brief Default constructor.
		 * 
		 * This constructor creates an empty natural_basis object.
		 * 
		 * Example usage:
		 * @code
		 * simcoon::natural_basis nb;
		 * @endcode
		 */
		natural_basis();

		/**
		 * @brief Constructor with parameters.
		 * 
		 * This constructor creates a natural_basis object from a std::vector of arma::vec containing the covariant vectors \f$\mathbf{g}_i\f$.
		 * 
		 * @param[in] g_i A std::vector of arma::vec containing the covariant vectors.
		 * 
		 * Example usage:
		 * @code
		 * std::vector<arma::vec> covariant_vectors = {arma::vec({1, 0, 0}), arma::vec({0, 1, 0}), arma::vec({0, 0, 1})};
		 * simcoon::natural_basis nb(covariant_vectors);
		 * @endcode
		 */
		natural_basis(const std::vector<arma::vec> &g_i);

		/**
		 * @brief Copy constructor.
		 * 
		 * This constructor creates a copy of a natural_basis object.
		 * 
		 * @param[in] other The natural_basis object to be copied.
		 * 
		 * Example usage:
		 * @code
		 * simcoon::natural_basis nb1;
		 * simcoon::natural_basis nb2(nb1);
		 * @endcode
		 */
		natural_basis(const natural_basis &other);

		/**
		 * @brief Destructor.
		 * 
		 * This destructor frees the memory used by the natural_basis object.
		 * 
		 * Example usage:
		 * @code
		 * {
		 *     simcoon::natural_basis nb;
		 * } // nb goes out of scope and is destroyed here
		 * @endcode
		 */
		virtual ~natural_basis();
    
		/**
		 * @brief update with a new set of covariant vectors.
		 * 
		 * This function updates the natural_basis object with a set of covariant vectors \f$\mathbf{g}_i\f$.
		 * 
		 * @param[in] g_i A std::vector of arma::vec containing the new covariant vectors.
		 * 
		 * Example usage:
		 * @code
		 * std::vector<arma::vec> new_covariant_vectors = {arma::vec({1, 0, 0}), arma::vec({0, 1, 0}), arma::vec({0, 0, 1})};
		 * simcoon::natural_basis nb;
		 * nb.update(new_covariant_vectors);
		 * @endcode
		 */
		virtual void update(const std::vector<arma::vec> &g_i);

		/**
		 * @brief Update using the transformation gradient.
		 * 
		 * This function updates the natural_basis object using the transformation gradient \f$\mathbf{F}\f$.
		 * 
		 * @param[in] F An arma::mat containing the transformation gradient.
		 * 
		 * Example usage:
		 * @code
		 * arma::mat F = arma::eye(3, 3); // Identity matrix as an example
		 * simcoon::natural_basis nb;
		 * nb.from_F(F);
		 * @endcode
		 */
		virtual void from_F(const arma::mat &F);
    
		/**
		 * @brief Assignment operator.
		 * 
		 * This operator assigns the values of another natural_basis object to the current object.
		 * 
		 * @param[in] other The natural_basis object from which the values are to be copied.
		 * 
		 * @return A reference to the current natural_basis object.
		 * 
		 * Example usage:
		 * @code
		 * simcoon::natural_basis nb1;
		 * simcoon::natural_basis nb2;
		 * nb2 = nb1;
		 * @endcode
		 */
		virtual natural_basis& operator = (const natural_basis& other);
        
		/**
		 * @brief Output stream operator.
		 * 
		 * This operator overloads the output stream operator to print the natural_basis object to an output stream.
		 * 
		 * @param[in] os The output stream.
		 * @param[in] nb The natural_basis object to be printed.
		 * 
		 * @return A reference to the output stream.
		 * 
		 * Example usage:
		 * @code
		 * simcoon::natural_basis nb;
		 * std::cout << nb;
		 * @endcode
		 */
		friend std::ostream& operator << (std::ostream& os, const natural_basis& nb);
};

} //namespace simcoon
