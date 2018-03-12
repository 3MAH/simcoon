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
///@brief Characteristics of an ellipsoidal thermomechanical phase, which hereditates from:
///- phase_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for ellipsoid===================================

//=====Public methods for ellipsoid============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
ellipsoid::ellipsoid() : geometry()
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
    
    a1 = 0.;
    a2 = 0.;
    a3 = 0.;
    
    psi_geom = 0.;
    theta_geom = 0.;
    phi_geom = 0.;
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
    ellipsoid::ellipsoid(const double &mconcentration, const int &mcoatingof, const int &mcoatedby, const double &ma1, const double &ma2, const double &ma3, const double &mpsi_geom, const double &mtheta_geom, const double &mphi_geom) : geometry(mconcentration)
//-------------------------------------------------------------
{	
    coatingof = mcoatingof;
    coatedby = mcoatedby;
    
    a1 = ma1;
    a2 = ma2;
    a3 = ma3;
    
    psi_geom = mpsi_geom;
    theta_geom = mtheta_geom;
    phi_geom = mphi_geom;
}

/*!
  \brief Copy constructor
  \param s ellipsoid_characteristics object to duplicate
*/

//------------------------------------------------------
ellipsoid::ellipsoid(const ellipsoid& sv) : geometry(sv)
//------------------------------------------------------
{
    coatingof = sv.coatingof;
    coatedby = sv.coatedby;
    
    a1 = sv.a1;
    a2 = sv.a2;
    a3 = sv.a3;
    
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
ellipsoid::~ellipsoid() {}
//-------------------------------------

/*!
  \brief Standard operator = for ellipsoid_characteristics
*/
    
    
//----------------------------------------------------------------------
ellipsoid& ellipsoid::operator = (const ellipsoid& sv)
//----------------------------------------------------------------------
{
    concentration = sv.concentration;
    
    coatingof = sv.coatingof;
    coatedby = sv.coatedby;

    a1 = sv.a1;
    a2 = sv.a2;
    a3 = sv.a3;
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const ellipsoid& sv)
//--------------------------------------------------------------------------
{    
    s << "Display ellipsoid properties\n";
    s << "Volume fraction = " << sv.concentration << "\n";
    s << "local geometrical orientation: psi = " << sv.psi_geom << "\t theta = " << sv.theta_geom << "\t phi = " << sv.phi_geom << "\n";
    s << "semi-axes: " << sv.a1 << "\t" << sv.a2 << "\t" << sv.a3 << "\n";
    s << "is a coating of the phase: " << sv.coatingof << "\n";
    s << "is coated by the phase: " << sv.coatedby << "\n";
    s << "\n\n";

	return s;
}

} //namespace simcoon
