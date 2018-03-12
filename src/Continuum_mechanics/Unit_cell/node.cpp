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

///@file node.cpp
///@brief Implementation of functions relative to a node (in the Finite Element sense)
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for material_characteristics===================================

//=====Public methods for material_characteristics============================================

/*!
 \brief default constructor
 */
    
//-------------------------------------------------------------
Node::Node()
//-------------------------------------------------------------
{
    number = 0;
    coords = Point(0.,0.,0.);
}

/*!
 \brief constructor with parameters
 */
    
//-------------------------------------------------------------
Node::Node(const unsigned int &mnumber, const Point &mcoords)
//-------------------------------------------------------------
{
    number = mnumber;
    coords = mcoords;
}

/*!
 \brief copy constructor
 */

//-------------------------------------------------------------
Node::Node(const Node & n)
//-------------------------------------------------------------
{
    number = n.number;
    coords= n.coords;
}


/*!
 \brief Destructor
 
 Deletes Node, the vectors and matrix.
 */

//-------------------------------------
Node::~Node() {}
//-------------------------------------
    
//-------------------------------------------------------------
void Node::initialize(const unsigned int &mnumber, const Point &mcoords)
//-------------------------------------------------------------
{
    number = mnumber;
    coords = mcoords;
}

//----------------------------------------------------------------------
Node& Node::operator = (const Node& n)
//----------------------------------------------------------------------
{
    number = n.number;
    coords= n.coords;
    
    return *this;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream &s, const Node &n)
//--------------------------------------------------------------------------
{
    s << n.number << ",\t" << n.coords.x() << ",\t" << n.coords.y() << ",\t" << n.coords.z() << "\n";
	return s;
}
    
} //namespace simcoon
