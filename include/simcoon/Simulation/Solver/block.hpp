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

///@file block.hpp
///@brief object that defines a block
///@version 1.0

#pragma once
#include <iostream>
#include <memory>
#include "step.hpp"

namespace simcoon{

//======================================
class block
//======================================
{
	private:

	protected:

	public :
		int number;
        int nstep;
        int ncycle;
        int type;
    
        std::vector<std::shared_ptr<step> > steps;
    
        block(); 	//default constructor
		block(int, int, int, int);	//constructor - allocates memory for the step vector
        block(int, int, int, int, const std::vector<std::shared_ptr<step> > &); //Constructor with parameters
		block(const block&);	//Copy constructor
		~block();
		
		void generate();
//		void initialize();
    
		virtual block& operator = (const block&);
		
        friend std::ostream& operator << (std::ostream&, const block&);
};

} //namespace simcoon
