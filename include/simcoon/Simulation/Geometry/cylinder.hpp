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

///@file ellipsoid.hpp
///@brief Characteristics of an ellipsoidal phase, which hereditates from:
///-phase characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace simcoon{

/**
 * @file cylinder.hpp
 * @brief Inclusion geometry functions.
 */

/** @addtogroup geometry
 *  @{
 */


//======================================
class cylinder : public geometry
//======================================
{
	private:

	protected:

	public :

        int coatingof;
        int coatedby;
    
        double L;  //geometric parameter of the cylinder (Length)
        double R;  //geometric parameter of the cylinder (Radius)
        
        double psi_geom;    //geometric orientation of the cylinder psi
        double theta_geom;  //geometric orientation of the cylinder theta
        double phi_geom;    //geometric orientation of the cylinder phi
    
		cylinder(); 	//default constructor
    		cylinder(const double &Lval, const int &coatingof, const int &coatedby, const double &Rval, const double &psi_geom, const double &theta_geom, const double &phi_geom, const double &dummy);

		cylinder(const cylinder&);	//Copy constructor
        virtual ~cylinder();
    
		virtual cylinder& operator = (const cylinder&);
		
        friend std::ostream& operator << (std::ostream& os, const cylinder &cyl);
};


/** @} */ // end of geometry group

} //namespace simcoon
