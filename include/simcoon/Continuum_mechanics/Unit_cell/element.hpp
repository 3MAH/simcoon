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

///@file element.hpp
///@brief Declaration of an element object (in the Finite Element sense)
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace simcoon{
    
//======================================
class Element
//======================================
{
    
    private:
    
    protected:
    
    public :
    
        std::string type;
        unsigned int number;
        std::vector<Node> nodes;
    
        Element(); 	//default constructor
        Element(const std::string &, const unsigned int &, const std::vector<Node> &);	//Constructor with parameters
        Element(const Element&);	//Copy constructor
        virtual ~Element();
        
        void initialize(const std::string &, const unsigned int &, const std::vector<Node> &); // Construct from the default constructor

        virtual Element& operator = (const Element&);
    
        friend void write_aba_format(std::ostream&, const Element&);
        friend std::ostream& operator << (std::ostream&, const Element&);
};

inline bool operator==(const Element& lhs, const Element& rhs){return (lhs.number == rhs.number); }
inline bool operator!=(const Element& lhs, const Element& rhs){return !operator==(lhs,rhs);}
inline bool operator< (const Element& lhs, const Element& rhs){return (lhs.number < rhs.number); }
inline bool operator> (const Element& lhs, const Element& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const Element& lhs, const Element& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const Element& lhs, const Element& rhs){return !operator< (lhs,rhs);}
    
} //namespace simcoon

