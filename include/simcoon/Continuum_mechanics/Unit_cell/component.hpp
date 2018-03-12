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

///@file component.hpp
///@brief Characteristics of an equation component (node, DOF, coef)

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace simcoon{
    
//======================================
class component
//======================================
{
    private:
            
    protected:
            
    public :

        Node node;
        unsigned int dof;
        double coef;
    
        component(); 	//default constructor
        component(const Node &, const unsigned int &, const double &);	//constructor with parameters
        
        component(const component &);	//Copy constructor
        ~component();
        
        virtual component& operator = (const component&);
        friend std::ostream& operator << (std::ostream&, const component&);
};
inline bool operator==(const component& lhs, const component& rhs){return ((lhs.node == rhs.node)&&(lhs.dof == rhs.dof)); }
    
} //namespace simcoon