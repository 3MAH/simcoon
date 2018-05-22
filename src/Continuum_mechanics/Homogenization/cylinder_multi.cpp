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

///@file cylinder_multi.cpp
///@brief Micromechanical characteristics of a phase
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/cylinder_multi.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
    
//=====Private methods for cylinder_multi===================================

//=====Public methods for cylinder_multi====================================

/*!
  \brief default constructor
*/
    
//-------------------------------------------------------------
cylinder_multi::cylinder_multi() : phase_multi(), T_loc(6,6), T(6,6), A_loc(6,6), B_loc(6,6)
//-------------------------------------------------------------
{
    //This calls only the constructor of the two matrix A & B
    
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
cylinder_multi::cylinder_multi(const mat &mA, const mat &mA_start, const mat &mB, const mat &mB_start, const vec &mA_in, const mat &mA_loc, const mat &mB_loc, const mat &mT_loc, const mat &mT) : phase_multi(mA, mA_start, mB, mB_start, mA_in), T_loc(6,6), T(6,6), A_loc(6,6), B_loc(6,6)
//-------------------------------------------------------------
{
    T_loc = mT_loc;
    T = mT;
    
    A_loc = mA_loc;
    B_loc = mB_loc;
}

/*!
  \brief Copy constructor
  \param s cylinder_multi object to duplicate
*/
    
//------------------------------------------------------
cylinder_multi::cylinder_multi(const cylinder_multi& pc) : phase_multi(pc), T_loc(6,6), T(6,6), A_loc(6,6), B_loc(6,6)
//------------------------------------------------------
{
    T_loc = pc.T_loc;
    T = pc.T;
    
    A_loc = pc.A_loc;
    B_loc = pc.B_loc;
}

/*!
  \brief Destructor

  Deletes phase_multi (the arma::mat).
*/

//-------------------------------------
cylinder_multi::~cylinder_multi() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_multi
*/

        
//----------------------------------------------------------------------
cylinder_multi& cylinder_multi::operator = (const cylinder_multi& pc)
//----------------------------------------------------------------------
{
    A = pc.A;
    B = pc.B;

    A_start = pc.A_start;
    B_start = pc.B_start;
    
    A_loc = pc.A_loc;
    B_loc = pc.B_loc;
    
    T_loc = pc.T_loc;
    T = pc.T;
    
	return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const cylinder_multi& pc)
//--------------------------------------------------------------------------
{
	s << "Display phase multi:\n";
	s << "Display strain concentration tensor (global coordinates):\n";
    s << pc.A;
    s << "Display strain concentration tensor (local coordinates):\n";
    s << pc.A_loc;
    s << "Display stress concentration tensor (global coordinates):\n";
    s << pc.B;
    s << "Display stress concentration tensor (local coordinates):\n";
    s << pc.B_loc;

    s << "Display Interaction concentration tensor (local coordinates):\n";
    s << pc.T_loc;
    s << "Display Interaction concentration tensor (global coordinates):\n";
    s << pc.T;
    
    
    s << "\n\n";

	return s;
}

} //namespace simcoon
