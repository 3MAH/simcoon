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

///@file geometry.hpp
///@brief Class that define the geometry of a phase
///@version 1.0

#pragma once

#include <iostream>
#include <string>

namespace simcoon{

//======================================
class geometry
//======================================
{
	private:

	protected:

	public :

        double concentration;
    
		geometry(); 	//default constructor
		geometry(const double &);	//constructor with parameters - allocates memory for statev

		geometry(const geometry&);	//Copy constructor
        virtual ~geometry();
    
		virtual geometry& operator = (const geometry&);
		
        friend std::ostream& operator << (std::ostream&, const geometry&);
};

} //namespace simcoon
