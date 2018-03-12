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

//======================================
class ellipsoid : public geometry
//======================================
{
	private:

	protected:

	public :

        int coatingof;
        int coatedby;
    
        double a1;  //geometric parameter of the ellipsoid (relative semi-axis 1)
        double a2;  //geometric parameter of the ellipsoid (relative semi-axis 2)
        double a3;  //geometric parameter of the ellipsoid (relative semi-axis 3)
        
        double psi_geom;    //geometric orientation of the ellipsoid psi
        double theta_geom;  //geometric orientation of the ellipsoid theta
        double phi_geom;    //geometric orientation of the ellipsoid phi
    
		ellipsoid(); 	//default constructor
        ellipsoid(const double &, const int &, const int &, const double &,const double &, const double &, const double &, const double &, const double &);

		ellipsoid(const ellipsoid&);	//Copy constructor
        virtual ~ellipsoid();
    
		virtual ellipsoid& operator = (const ellipsoid&);
		
        friend std::ostream& operator << (std::ostream&, const ellipsoid&);
};

} //namespace simcoon
