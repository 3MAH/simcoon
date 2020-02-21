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

///@file equation.hpp
///@brief Characteristics of an equation between nodes

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/component.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace simcoon{
    
//======================================
class equation
//======================================
{
    private:
            
    protected:
            
    public :

        std::vector<component> components;
    
        equation(); 	//default constructor
        equation(const std::vector<component> &);	//constructor with parameters
        
        equation(const equation &);	//Copy constructor
        virtual ~equation();
        
        virtual equation& operator = (const equation&);
        friend std::ostream& operator << (std::ostream&, const equation&);
};
    
} //namespace simcoon
