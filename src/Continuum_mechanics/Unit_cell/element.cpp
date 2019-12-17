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

///@file element.cpp
///@brief Implementation of functions relative to an element (in the Finite Element sense)
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/Continuum_mechanics/Unit_cell/element.hpp>

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
Element::Element()
//-------------------------------------------------------------
{
    number = 0;
    type = "";
}

/*!
 \brief constructor with parameters
 */
    
//-------------------------------------------------------------
Element::Element(const std::string &mtype, const unsigned int &mnumber, const std::vector<Node> &mnodes)
//-------------------------------------------------------------
{
    type = mtype;
    number = mnumber;
    nodes = mnodes;
}

/*!
 \brief copy constructor
 */

//-------------------------------------------------------------
Element::Element(const Element & e)
//-------------------------------------------------------------
{
    type= e.type;
    number = e.number;
    nodes= e.nodes;
}


/*!
 \brief Destructor
 
 Deletes Node, the vectors and matrix.
 */

//-------------------------------------
Element::~Element() {}
//-------------------------------------
    
//-------------------------------------------------------------
void Element::initialize(const std::string &mtype, const unsigned int &mnumber, const std::vector<Node> &mnodes)
//-------------------------------------------------------------
{
    type = mtype;
    number = mnumber;
    nodes = mnodes;
}

//----------------------------------------------------------------------
Element& Element::operator = (const Element& e)
//----------------------------------------------------------------------
{
    type= e.type;
    number = e.number;
    nodes= e.nodes;
    
    return *this;
}

//-------------------------------------------------------------
void write_aba_format(ostream &s, const Element &e)
//-------------------------------------------------------------
{
    s << e.number;
    for(auto n : e.nodes) {
        s << ", " << n.number;
    }
    s << endl;
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream &s, const Element &e)
//--------------------------------------------------------------------------
{
    s << "element type: " << e.type << "\t" << e.number << "\t nodes list: ";
    for(auto n : e.nodes) {
        cout << n << ", ";
    }
    s << endl;
    return s;
}
    
} //namespace simcoon
