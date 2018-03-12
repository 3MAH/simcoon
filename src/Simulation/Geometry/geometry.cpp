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

///@file geometry.cpp
///@brief Characteristics of an geometry
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Simulation/Geometry/geometry.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for geometry===================================

//=====Public methods for geometry============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
geometry::geometry()
//-------------------------------------------------------------
{
    concentration = 0.;
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
geometry::geometry(const double &mconcentration)
//-------------------------------------------------------------
{	
    concentration = mconcentration;
}

/*!
  \brief Copy constructor
  \param s geometry characteristics object to duplicate
*/

//------------------------------------------------------
geometry::geometry(const geometry& sv)
//------------------------------------------------------
{
    concentration = sv.concentration;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
geometry::~geometry() {}
//-------------------------------------

/*!
  \brief Standard operator = for geometry
*/
    
    
//----------------------------------------------------------------------
geometry& geometry::operator = (const geometry& sv)
//----------------------------------------------------------------------
{
    concentration = sv.concentration;
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const geometry& sv)
//--------------------------------------------------------------------------
{
    s << "Display geometry properties\n";
    s << "Volume fraction = " << sv.concentration << "\n";
    s << "\n\n";

	return s;
}

} //namespace simcoon
