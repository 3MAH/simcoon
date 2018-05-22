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

///@file phase_multi.cpp
///@brief Micromechanical characteristics of a phase
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for phase_multi===================================

//=====Public methods for phase_multi====================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
phase_multi::phase_multi() : A(6,6), A_start(6,6), B(6,6), B_start(6,6), A_in(6)
//-------------------------------------------------------------
{
    
}

/*!
  \brief Constructor with parameters
*/

//-------------------------------------------------------------
phase_multi::phase_multi(const mat &mA, const mat &mA_start, const mat &mB, const mat &mB_start, const vec &mA_in) : A(6,6), A_start(6,6), B(6,6), B_start(6,6), A_in(6)
//-------------------------------------------------------------
{
    A = mA;
    B = mB;
    
    A_start = mA_start;
    B_start = mB_start;
    
    A_in = mA_in;
}

/*!
  \brief Copy constructor
  \param s phase_multi object to duplicate
*/
    
//------------------------------------------------------
phase_multi::phase_multi(const phase_multi& pc) : A(6,6), A_start(6,6), B(6,6), B_start(6,6), A_in(6)
//------------------------------------------------------
{
    A = pc.A;
    B = pc.B;
    
    A_start = pc.A_start;
    B_start = pc.B_start;
    
    A_in = pc.A_in;
}

/*!
  \brief Destructor

  Deletes phase_multi (the arma::mat).
*/

//-------------------------------------
phase_multi::~phase_multi() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_multi
*/

//-------------------------------------------------------------
void phase_multi::to_start()
//-------------------------------------------------------------
{
    A = A_start;
    B = B_start;
}

//-------------------------------------------------------------
void phase_multi::set_start()
//-------------------------------------------------------------
{
    A_start = A;
    B_start = B;
}
    
//----------------------------------------------------------------------
phase_multi& phase_multi::operator = (const phase_multi& pc)
//----------------------------------------------------------------------
{
    A = pc.A;
    B = pc.B;
    
    A_start = pc.A_start;
    B_start = pc.B_start;
    
    A_in = pc.A_in;
    
	return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const phase_multi& pc)
//--------------------------------------------------------------------------
{
	s << "Display phase multi:\n";
	s << "Display strain concentration tensor:\n";
    s << pc.A;
    s << "Display stress concentration tensor:\n";
    s << pc.B;
    s << "Display inelastic strain concentration vector:\n";
    s << pc.A_in;
    
    
    s << "\n\n";

	return s;
}

} //namespace simcoon
