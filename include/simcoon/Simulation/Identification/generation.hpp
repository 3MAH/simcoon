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

///@file generation.hpp
///@brief generation for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>
#include "individual.hpp"

namespace simcoon{

//======================================
class generation
//======================================
{
	private:

	protected:

	public :
    
        std::vector<individual> pop;
		
		generation(); 	//default constructor
		generation(const int&, const int&, int&, const double & = 0.);	//constructor - allocates memory for statev
		generation(const generation&);	//Copy constructor
		~generation();
		
		int size() const {return pop.size();}       // returns the number of individuals

		void construct(const int&, const int&, int&, const double & = 0.);
		void classify();
		void newid(int &);
		void destruct();

		virtual generation& operator = (const generation&);
		
		friend std::ostream& operator << (std::ostream&, const generation&);
};

} //namespace simcoon
