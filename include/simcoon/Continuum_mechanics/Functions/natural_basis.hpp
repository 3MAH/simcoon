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

///@file natural_basis.hpp
///@brief Natural basis objects
///@version 1.0

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
    
        std::vector<arma::vec> g_i; // Covariant Vectors
        std::vector<arma::vec> g0i; // Contravariant Vectors
    
        arma::mat g_ij; // Covariant components of the metric tensor
        arma::mat g0ij; // Contravariant components of the metric tensor
        
		natural_basis(); 	//default constructor
		natural_basis(const std::vector<arma::vec> &); //Constructor with parameters
		natural_basis(const natural_basis &);	//Copy constructor
		virtual ~natural_basis();
    
        virtual void update(const std::vector<arma::vec> &); //update with a new set of covariant vectors
        virtual void from_F(const arma::mat &F); //update using the transformation gradient
    
		virtual natural_basis& operator = (const natural_basis&);
        
        friend std::ostream& operator << (std::ostream&, const natural_basis&);
};

} //namespace simcoon
