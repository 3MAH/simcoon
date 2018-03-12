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

///@file variant.cpp
///@brief Class that defines characteristics of a variant of martensite
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <armadillo>
#include <simcoon/Continuum_Mechanics/Material/variant.hpp>
#include <simcoon/Continuum_Mechanics/Material/crystallo.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for variant===================================

//=====Public methods for variant====================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
variant::variant() : n(3), m(3), R(3,3), ETn(6)
//-------------------------------------------------------------
{

	n = zeros(3);
	m = zeros(3);
	R = zeros(3,3);
	g = 0.;
	ETn = zeros(6);

}

/*!
  \brief Constructor with parameters
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
variant::variant(vec mn, vec mm, double mg) : n(3), m(3), R(3,3), ETn(6)
//-------------------------------------------------------------
{	
	assert (mn.size() == 3);
	assert (mm.size() == 3);
	
	n = mn;
	m = mm;
	R = Schmid(n, m);
	g = mg;
	ETn = g*Schmid_v(n, m);
}

/*!
  \brief Copy constructor
  \param s variant object to duplicate
*/

//------------------------------------------------------
variant::variant(const variant& sv) : n(3), m(3), R(3,3), ETn(6)
//------------------------------------------------------
{
	n = sv.n;
	m = sv.m;
	R = sv.R;
	g = sv.g;
	ETn = sv.ETn;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
variant::~variant()
//-------------------------------------
{

}

/*!
  \brief Standard operator = for variant
*/

//-------------------------------------------------------------
void variant::build(const double& g)
//-------------------------------------------------------------
{
	assert (n.size() == 3);
	assert (m.size() == 3);	
	assert (g > 0.);
	
	R = Schmid(n, m);
	ETn = g*Schmid_v(n, m);
}

//----------------------------------------------------------------------
variant& variant::operator = (const variant& sv)
//----------------------------------------------------------------------
{
	n = sv.n;
	m = sv.m;
	R = sv.R;
	g = sv.g;
	ETn = sv.ETn;

	return *this;

}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const variant& sv)
//--------------------------------------------------------------------------
{
	s << "n: \n" << sv.n << "\n";
	s << "m: \n" << sv.m << "\n";
	s << "R: \n" << sv.R << "\n";
	s << "g: \n" << sv.g << "\n";
	s << "ETn: \n" << sv.ETn << "\n";
	s << "\n";

	return s;
}

} //namespace simcoon
