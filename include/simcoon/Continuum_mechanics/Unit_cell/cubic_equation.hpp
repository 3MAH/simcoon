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

///@file cubic_equation.hpp
///@brief Class that gathers periodic equations between nodes in a cubic Representative Volume Element

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace simcoon{
    
//======================================
class cubic_equation
//======================================
{
    private:
            
    protected:
            
    public :
    
    std::vector<Node> CD_nodes;
    std::vector<int> CD;
    std::vector<std::string> CD_set_name;
    std::vector<int> list_dofs;
    
        std::vector<equation> Corner_listXmYpZm;
        std::vector<equation> Corner_listXpYmZm;
        std::vector<equation> Corner_listXpYpZm;
        std::vector<equation> Corner_listXmYmZp;
        std::vector<equation> Corner_listXmYpZp;
        std::vector<equation> Corner_listXpYmZp;
        std::vector<equation> Corner_listXpYpZp;
    
        std::vector<std::vector<equation> > Edge_listXpYm;
        std::vector<std::vector<equation> > Edge_listXpYp;
        std::vector<std::vector<equation> > Edge_listXmYp;
        std::vector<std::vector<equation> > Edge_listXpZm;
        std::vector<std::vector<equation> > Edge_listXpZp;
        std::vector<std::vector<equation> > Edge_listXmZp;
        std::vector<std::vector<equation> > Edge_listYpZm;
        std::vector<std::vector<equation> > Edge_listYpZp;
        std::vector<std::vector<equation> > Edge_listYmZp;
        std::vector<std::vector<equation> > Face_listXp;
        std::vector<std::vector<equation> > Face_listYp;
        std::vector<std::vector<equation> > Face_listZp;

        cubic_equation(); 	//default constructor
        cubic_equation(const cubic_mesh &, const cubic_mesh &, const unsigned int &, const unsigned int &); 	//Parameter constructor
        
        cubic_equation(const cubic_equation &);	//Copy constructor
        ~cubic_equation();
    
        void construct(const cubic_mesh &, const cubic_mesh &, const unsigned int &, const unsigned int &);
        
        //    virtual cubic_equation& operator = (const cubic_equation&);
        friend std::ostream& operator << (std::ostream&, const cubic_equation&);
};
    
} //namespace simcoon
