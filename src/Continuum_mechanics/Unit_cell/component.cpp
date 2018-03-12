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

///@file component.cpp
///@brief Implementation of the component of an equation
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Unit_cell/component.hpp>
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
component::component()
//-------------------------------------------------------------
{
    dof = 0;
    coef = 0.;
}
    
/*!
 \brief Constructor
 \param nodes : list of nodes that are linked in the equation
 \param coefs : coefs that links the nodes in the equation
 \n\n
 */

//-------------------------------------------------------------
component::component(const Node &mnode, const unsigned int &mdof, const double &mcoef)
//-------------------------------------------------------------
{
    node = mnode;
    dof = mdof;
    coef = mcoef;
}

/*!
 \brief Copy constructor
 \param eq equation object to duplicate
 */
    
//------------------------------------------------------
component::component(const component& eqc)
//------------------------------------------------------
{
    node = eqc.node;
    dof = eqc.dof;
    coef = eqc.coef;
}

/*!
 \brief Destructor
 
 Deletes equation
 */

//-------------------------------------
component::~component() {}
//-------------------------------------
    
//----------------------------------------------------------------------
component& component::operator = (const component& eqc)
//----------------------------------------------------------------------
{
    node = eqc.node;
    dof = eqc.dof;
    coef = eqc.coef;
    
    return *this;
}
    
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const component& eqc)
//--------------------------------------------------------------------------
{
    s << "Display info on the component:\n";
    s << "Node: " << eqc.node << "\n";
    s << "Degree of freedom: " << eqc.dof << "\n";
    s << "coef: " << eqc.coef << endl;
    return s;
}
    
} //namespace simcoon
