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

///@file isv.hpp
///@brief Internal state variables in a defined coordinate system:
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class isv
//======================================
{
	private:

	protected:

	public :
		
        double p ;           //V1
        arma::vec E_ir;      //V2
        arma::vec X;      //V3, the most common
        arma::field<arma::vec>  //V4+

        arma::vec Lambda    //Lambda
    
		isv(); 	//default constructor
        isv(const double &, const arma::vec &, const arma::vec &, const arma::field<arma::vec> &; //Constructor with parameters
		isv(const isv &);	//Copy constructor
		virtual ~isv();
		
		virtual isv& operator = (const isv&);

        friend std::ostream& operator << (std::ostream&, const isv&);
};

} //namespace simcoon
