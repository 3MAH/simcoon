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

///@file individual.hpp
///@brief individual for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class individual
//======================================
{
	private:

	protected:

	public :
		int np;
        double cout;
		int id;
		int rank;		
        arma::vec p;
        double lambda;   //The step or Lambda applied
		
		individual(); 	//default constructor
		individual(const int&,const int&,const double&);	//constructor - allocates memory for statev
		individual(const individual&);	//Copy constructor
		~individual();

		void construct();

		virtual individual& operator = (const individual&);
		
        friend std::ostream& operator << (std::ostream&, const individual&);
};
    
} //namespace simcoon
