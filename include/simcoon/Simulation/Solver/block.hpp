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
		unsigned int number;
        unsigned int nstep;
        unsigned int ncycle;
        unsigned int type;
        unsigned int control_type;
    
        std::vector<std::shared_ptr<step> > steps;
    
        block(); 	//default constructor
		block(const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&);	//constructor - allocates memory for the step vector
        block(const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&, const unsigned int&, const std::vector<std::shared_ptr<step> > &); //Constructor with parameters
		block(const block&);	//Copy constructor
		virtual ~block();
		
		void generate();
//		void initialize();
    
		virtual block& operator = (const block&);
		
        friend std::ostream& operator << (std::ostream&, const block&);
};

} //namespace simcoon
