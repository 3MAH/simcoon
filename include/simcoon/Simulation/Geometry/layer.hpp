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

///@file layer.hpp
///@brief Characteristics of an layer geometry, which hereditates from:
///-geometry
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace simcoon{

//======================================
class layer : public geometry
//======================================
{
	private:

	protected:

	public :
    
        int layerup;
        int layerdown;
    
        //a stack of layer have the same direction, it is attributed to each one of them since ther might be several stacks (one per inclusion, for instance..)
        double psi_geom;
        double theta_geom;
        double phi_geom;    
    
		layer(); 	//default constructor
        layer(const double &, const int &, const int &, const double &,const double &, const double &); //Constructor with parameters

		layer(const layer&);	//Copy constructor
        virtual ~layer();
    
		virtual layer& operator = (const layer&);
		
        friend std::ostream& operator << (std::ostream&, const layer&);
};

} //namespace simcoon
