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

///@file node.hpp
///@brief Declaration of a node object (in the Finite Element sense)
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

namespace simcoon{
    
//======================================
class Node
//======================================
{
    
    private:
    
    protected:
    
    public :
    
        unsigned int number;
        Point coords;
    
        Node(); 	//default constructor
        Node(const unsigned int &, const Point &);	//Constructor with parameters
        Node(const Node&);	//Copy constructor
        virtual ~Node();
        
        void initialize(const unsigned int &, const Point &); // Construct from the default constructor
        
        virtual Node& operator = (const Node&);
    
        friend std::ostream& operator << (std::ostream&, const Node&);
};

inline bool operator==(const Node& lhs, const Node& rhs){return (lhs.number == rhs.number); }
inline bool operator!=(const Node& lhs, const Node& rhs){return !operator==(lhs,rhs);}
inline bool operator< (const Node& lhs, const Node& rhs){return (lhs.number < rhs.number); }
inline bool operator> (const Node& lhs, const Node& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const Node& lhs, const Node& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const Node& lhs, const Node& rhs){return !operator< (lhs,rhs);}
    
} //namespace simcoon

