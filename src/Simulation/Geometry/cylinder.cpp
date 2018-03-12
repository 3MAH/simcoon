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

///@file ellipsoid_characteristics.cpp
///@brief Characteristics of an cylinder thermomechanical phase, which hereditates from:
///- phase_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for cylinder===================================

//=====Public methods for cylinder============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
cylinder::cylinder() : geometry()
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
    
    L = 0.;
    R = 0.;
    
    psi_geom = 0.;
    theta_geom = 0.;
    phi_geom = 0.;
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
    cylinder::cylinder(const double &mconcentration, const int &mcoatingof, const int &mcoatedby, const double &mL, const double &mR, const double &mpsi_geom, const double &mtheta_geom, const double &mphi_geom) : geometry(mconcentration)
//-------------------------------------------------------------
{	
    coatingof = mcoatingof;
    coatedby = mcoatedby;
    
    L = mL;
    R = mR;
    
    psi_geom = mpsi_geom;
    theta_geom = mtheta_geom;
    phi_geom = mphi_geom;
}

/*!
  \brief Copy constructor
  \param s cylinder_characteristics object to duplicate
*/

//------------------------------------------------------
cylinder::cylinder(const cylinder& sv) : geometry(sv)
//------------------------------------------------------
{
    coatingof = sv.coatingof;
    coatedby = sv.coatedby;
    
    L = sv.L;
    R = sv.R;
    
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
cylinder::~cylinder() {}
//-------------------------------------

/*!
  \brief Standard operator = for cylinder_characteristics
*/
    
    
//----------------------------------------------------------------------
cylinder& cylinder::operator = (const cylinder& sv)
//----------------------------------------------------------------------
{
    concentration = sv.concentration;
    
    coatingof = sv.coatingof;
    coatedby = sv.coatedby;

    L = sv.L;
    R = sv.R;
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const cylinder& sv)
//--------------------------------------------------------------------------
{    
    s << "Display cylinder properties\n";
    s << "Volume fraction = " << sv.concentration << "\n";
    s << "local geometrical orientation: psi = " << sv.psi_geom << "\t theta = " << sv.theta_geom << "\t phi = " << sv.phi_geom << "\n";
    s << "dimensions: Length = " << sv.L << "\t radius = " << sv.R << "\n";
    s << "is a coating of the phase: " << sv.coatingof << "\n";
    s << "is coated by the phase: " << sv.coatedby << "\n";
    s << "\n\n";

	return s;
}

} //namespace simcoon
