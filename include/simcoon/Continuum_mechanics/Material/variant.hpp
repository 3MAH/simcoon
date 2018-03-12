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
		
		variant(); 	//default constructor
		variant(arma::vec, arma::vec, double); //Constructor with parameters
		variant(const variant &);	//Copy constructor
		~variant();
		
		void build(const double &);
		virtual variant& operator = (const variant&);
		
        friend std::ostream& operator << (std::ostream&, const variant&);
};

} //namespace simcoon
