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

///@file layer.cpp
///@brief Characteristics of an layeral thermomechanical phase, which hereditates from:
///- phase_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for layer===================================

//=====Public methods for layer============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
layer::layer() : geometry()
//-------------------------------------------------------------
{
    layerup = 0;
    layerdown = 0;
    
    psi_geom = 0.;
    theta_geom = 0.;
    phi_geom = 0.;    
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
layer::layer(const double &mconcentration, const int &mlayerup, const int &mlayerdown, const double &mpsi_geom, const double &mtheta_geom, const double &mphi_geom) : geometry(mconcentration)
//-------------------------------------------------------------
{	
    layerup = mlayerup;
    layerdown = mlayerdown;
    
    psi_geom = mpsi_geom;
    theta_geom = mtheta_geom;
    phi_geom = mphi_geom;    
}

/*!
  \brief Copy constructor
  \param s layer characteristics object to duplicate
*/

//------------------------------------------------------
layer::layer(const layer& sv) : geometry(sv)
//------------------------------------------------------
{
    layerup = sv.layerup;
    layerdown = sv.layerdown;
    
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;    
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
layer::~layer() {}
//-------------------------------------

/*!
  \brief Standard operator = for layer
*/
    
    
//----------------------------------------------------------------------
layer& layer::operator = (const layer& sv)
//----------------------------------------------------------------------
{
    concentration = sv.concentration;
    
    layerup = sv.layerup;
    layerdown = sv.layerdown;
    
    psi_geom = sv.psi_geom;
    theta_geom = sv.theta_geom;
    phi_geom = sv.phi_geom;    
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const layer& sv)
//--------------------------------------------------------------------------
{    
    s << "Display layer properties\n";
    s << "Volume fraction = " << sv.concentration << "\n";
    s << "local geometrical orientation: psi = " << sv.psi_geom << "\t theta = " << sv.theta_geom << "\t phi = " << sv.phi_geom << "\n";
    s << "upper layer: " << sv.layerup << "\n";
    s << "lower layer: " << sv.layerdown << "\n";
    s << "\n\n";

	return s;
}

} //namespace simcoon
