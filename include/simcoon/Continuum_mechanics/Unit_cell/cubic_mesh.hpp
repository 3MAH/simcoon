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

///@file cubic_mesh.hpp
///@brief Characteristics of a cubic mesh
///@version 1.0

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
class cubic_mesh
//======================================
{
    private:
            
    protected:
            
    public :
            
        bool is_perio;

        static std::string set_name_all;
        static std::vector<std::string> set_name_corners;
        static std::vector<std::string> set_name_edges;
        static std::vector<std::string> set_name_faces;

        std::vector<Node> Node_list;
        std::vector<Point> Point_list;

        std::string Node_list_name;

        Node center_node;

        Kernel::Iso_cuboid_3 cuboid; //For the bounding box
        double volume; //The characeristic size of the box
        Kernel::FT Dx;
        Kernel::FT Dy;
        Kernel::FT Dz;
        double size_box;

        std::shared_ptr<Node> Corner_listXmYmZm;
        std::shared_ptr<Node> Corner_listXmYpZm;
        std::shared_ptr<Node> Corner_listXpYmZm;
        std::shared_ptr<Node> Corner_listXpYpZm;
        std::shared_ptr<Node> Corner_listXmYmZp;
        std::shared_ptr<Node> Corner_listXmYpZp;
        std::shared_ptr<Node> Corner_listXpYmZp;
        std::shared_ptr<Node> Corner_listXpYpZp;

        std::shared_ptr<std::vector<Node> > Edge_listXmYm;
        std::shared_ptr<std::vector<Node> > Edge_listXpYm;
        std::shared_ptr<std::vector<Node> > Edge_listXpYp;
        std::shared_ptr<std::vector<Node> > Edge_listXmYp;
        std::shared_ptr<std::vector<Node> > Edge_listXmZm;
        std::shared_ptr<std::vector<Node> > Edge_listXpZm;
        std::shared_ptr<std::vector<Node> > Edge_listXpZp;
        std::shared_ptr<std::vector<Node> > Edge_listXmZp;
        std::shared_ptr<std::vector<Node> > Edge_listYmZm;
        std::shared_ptr<std::vector<Node> > Edge_listYpZm;
        std::shared_ptr<std::vector<Node> > Edge_listYpZp;
        std::shared_ptr<std::vector<Node> > Edge_listYmZp;

        std::shared_ptr<std::vector<Node> > Face_listXm;
        std::shared_ptr<std::vector<Node> > Face_listYm;
        std::shared_ptr<std::vector<Node> > Face_listZm;
        std::shared_ptr<std::vector<Node> > Face_listXp;
        std::shared_ptr<std::vector<Node> > Face_listYp;
        std::shared_ptr<std::vector<Node> > Face_listZp;

        std::vector< std::shared_ptr<Node> > list_of_corners;
        std::vector< std::shared_ptr<std::vector<Node> > > list_of_edges;
        std::vector< std::shared_ptr<std::vector<Node> > > list_of_faces;

        cubic_mesh(); 	//default constructor
        cubic_mesh(const std::vector<Node> &, const std::string &);	//constructor with parameters

        cubic_mesh(const cubic_mesh &);	//Copy constructor
        ~cubic_mesh();

        void initialize(const std::vector<Node> &, const std::string &); // Construct from the default constructor
        void construct();                                // Construct the lists with make_shared

        //    void check_duplicates(const double &);
        void get_domain();
        void find_pairs();
        void construct_lists();

        //    virtual cubic_mesh& operator = (const cubic_mesh&);
        friend std::ostream& operator << (std::ostream&, const cubic_mesh&);
};
    
} //namespace simcoon