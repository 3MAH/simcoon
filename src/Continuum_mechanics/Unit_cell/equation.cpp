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

///@file equation.cpp
///@brief Characteristics of an equation between nodes
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for equation===================================

//=====Public methods for equation============================================

/*!
 \brief default constructor
 */
    
//-------------------------------------------------------------
equation::equation()
//-------------------------------------------------------------
{
}
    
/*!
 \brief Constructor
 \param nodes : list of nodes that are linked in the equation
 \param coefs : coefs that links the nodes in the equation
 \n\n
 */

//-------------------------------------------------------------
equation::equation(const std::vector<component> &mcomponents)
//-------------------------------------------------------------
{
    components = mcomponents;
}

/*!
 \brief Copy constructor
 \param eq equation object to duplicate
 */
    
//------------------------------------------------------
equation::equation(const equation& eq)
//------------------------------------------------------
{
    components = eq.components;
}

/*!
 \brief Destructor
 
 Deletes equation
 */

//-------------------------------------
equation::~equation() {}
//-------------------------------------
    
//----------------------------------------------------------------------
equation& equation::operator = (const equation& eq)
//----------------------------------------------------------------------
{
    components = eq.components;
    
    return *this;
}
    
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const equation& eq)
//--------------------------------------------------------------------------
{
    s << "Display info on the equation:\n";
    for (auto n:eq.components) {
        s << n;
    }
    s << endl;
    return s;
}
    
} //namespace simcoon
