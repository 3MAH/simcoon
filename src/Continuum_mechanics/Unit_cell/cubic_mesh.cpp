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

///@file cubic_mesh.cpp
///@brief Characteristics of a mesh for a Representative Volume Element.
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/squared_distance_3.h> //for 3D functions

#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for cubic_mesh===================================

//=====Public methods for cubic_mesh============================================

/*!
 \brief default constructor
 */

std::string cubic_mesh::set_name_all = "AllNodes";
std::vector<std::string> cubic_mesh::set_name_corners = {"CornerXmYmZm", "CornerXpYmZm", "CornerXpYpZm", "CornerXmYpZm", "CornerXmYmZp", "CornerXpYmZp", "CornerXpYpZp", "CornerXmYpZp"};
std::vector<std::string> cubic_mesh::set_name_edges = {"EdgeXmYm", "EdgeXpYm", "EdgeXpYp", "EdgeXmYp", "EdgeXmZm", "EdgeXpZm", "EdgeXpZp", "EdgeXmZp", "EdgeYmZm", "EdgeYpZm", "EdgeYpZp", "EdgeYmZp"};
std::vector<std::string> cubic_mesh::set_name_faces = {"FaceXm", "FaceXp", "FaceYm", "FaceYp", "FaceZm", "FaceZp"};
    
    
//-------------------------------------------------------------
cubic_mesh::cubic_mesh()
//-------------------------------------------------------------
{
    is_perio = false;
    Node_list_name = "";
    volume = 0.;
    Dx = 0.;
    Dy = 0.;
    Dz = 0.;
    size_box = 0.;
}
    
/*!
 \brief Constructor
 \param Node_list : list of nodes of the cubic mesh
 \param Node_list_name : Output file (.inp) of the normalized list of nodes and sets
 \n\n
 */

//-------------------------------------------------------------
cubic_mesh::cubic_mesh(const std::vector<Node> &mNode_list, const std::string &mNode_list_name)
//-------------------------------------------------------------
{
    is_perio = false;
    Node_list = mNode_list;
    Node_list_name = mNode_list_name;
    
    volume = 0.;
    Dx = 0.;
    Dy = 0.;
    Dz = 0.;
    size_box = 0.;
    
    construct();
}

/*!
 \brief Copy constructor
 \param cm cubic_mesh object to duplicate
 */
    
//------------------------------------------------------
cubic_mesh::cubic_mesh(const cubic_mesh& cm)
//------------------------------------------------------
{
    is_perio = cm.is_perio;
    Node_list = cm.Node_list;
    Node_list_name = cm.Node_list_name;
    
    cuboid = cm.cuboid;
    volume = cm.volume;
    Dx = cm.Dx;
    Dy = cm.Dy;
    Dz = cm.Dz;
    size_box = cm.size_box;
    
    center_node = cm.center_node;
    
    Edge_listXmYm = cm.Edge_listXmYm;
    Edge_listXpYm = cm.Edge_listXpYm;
    Edge_listXpYp = cm.Edge_listXpYp;
    Edge_listXmYp = cm.Edge_listXmYp;
    Edge_listXmZm = cm.Edge_listXmZm;
    Edge_listXpZm = cm.Edge_listXpZm;
    Edge_listXpZp = cm.Edge_listXpZp;
    Edge_listXmZp = cm.Edge_listXmZp;
    Edge_listYmZm = cm.Edge_listYmZm;
    Edge_listYpZm = cm.Edge_listYpZm;
    Edge_listYpZp = cm.Edge_listYpZp;
    Edge_listYmZp = cm.Edge_listYmZp;
    
    Face_listXm = cm.Face_listXm;
    Face_listYm = cm.Face_listYm;
    Face_listZm = cm.Face_listZm;
    Face_listXp = cm.Face_listXp;
    Face_listYp = cm.Face_listYp;
    Face_listZp = cm.Face_listZp;
    
    Corner_listXmYmZm = cm.Corner_listXmYmZm;
    Corner_listXmYpZm = cm.Corner_listXmYpZm;
    Corner_listXpYmZm = cm.Corner_listXpYmZm;
    Corner_listXpYpZm = cm.Corner_listXpYpZm;
    Corner_listXmYmZp = cm.Corner_listXmYmZp;
    Corner_listXmYpZp = cm.Corner_listXmYpZp;
    Corner_listXpYmZp = cm.Corner_listXpYmZp;
    Corner_listXpYpZp = cm.Corner_listXpYpZp;
    
    Faces = cm.Faces;
}

/*!
 \brief Destructor
 
 Deletes cubic_mesh, the vectors and matrix.
 */

//-------------------------------------
cubic_mesh::~cubic_mesh() {}
//-------------------------------------
    
//-------------------------------------------------------------
void cubic_mesh::initialize(const std::vector<Node> &mNode_list, const std::string &mNode_list_name)
//-------------------------------------------------------------
{
    Node_list = mNode_list;
    Node_list_name = mNode_list_name;
    
    construct();
}
    
//-------------------------------------------------------------
void cubic_mesh::construct()
//-------------------------------------------------------------
{
    Corner_listXmYmZm = std::make_shared<Node>();
    Corner_listXmYpZm = std::make_shared<Node>();
    Corner_listXpYmZm = std::make_shared<Node>();
    Corner_listXpYpZm = std::make_shared<Node>();
    Corner_listXmYmZp = std::make_shared<Node>();
    Corner_listXmYpZp = std::make_shared<Node>();
    Corner_listXpYmZp = std::make_shared<Node>();
    Corner_listXpYpZp = std::make_shared<Node>();
    
    Edge_listXmYm = std::make_shared<std::vector<Node> >();
    Edge_listXpYm = std::make_shared<std::vector<Node> >();
    Edge_listXpYp = std::make_shared<std::vector<Node> >();
    Edge_listXmYp = std::make_shared<std::vector<Node> >();
    Edge_listXmZm = std::make_shared<std::vector<Node> >();
    Edge_listXpZm = std::make_shared<std::vector<Node> >();
    Edge_listXpZp = std::make_shared<std::vector<Node> >();
    Edge_listXmZp = std::make_shared<std::vector<Node> >();
    Edge_listYmZm = std::make_shared<std::vector<Node> >();
    Edge_listYpZm = std::make_shared<std::vector<Node> >();
    Edge_listYpZp = std::make_shared<std::vector<Node> >();
    Edge_listYmZp = std::make_shared<std::vector<Node> >();
    
    Face_listXm = std::make_shared<std::vector<Node> >();
    Face_listYm = std::make_shared<std::vector<Node> >();
    Face_listZm = std::make_shared<std::vector<Node> >();
    Face_listXp = std::make_shared<std::vector<Node> >();
    Face_listYp = std::make_shared<std::vector<Node> >();
    Face_listZp = std::make_shared<std::vector<Node> >();
}
    
//-------------------------------------------------------------
void cubic_mesh::get_domain()
//-------------------------------------------------------------
{
    for (auto n : Node_list) {
        Point_list.push_back(n.coords);
    }
    
    cuboid = CGAL::bounding_box(Point_list.begin(), Point_list.end());
    
    Point cuboid_center(0.5*(cuboid.xmax()+cuboid.xmin()),0.5*(cuboid.ymax()+cuboid.ymin()),0.5*(cuboid.zmax()+cuboid.zmin()));

    Dx = cuboid.xmax()-cuboid.xmin();
    Dy = cuboid.ymax()-cuboid.ymin();
    Dz = cuboid.zmax()-cuboid.zmin();

    volume = CGAL::to_double(cuboid.volume());
    size_box = pow(volume,1./3.);
    center_node = closest_node(Node_list, cuboid_center);
    
    //Defines the points, faces and edges
    *Corner_listXmYmZm = find_corner(Node_list, cuboid[0], size_box);
    *Corner_listXmYpZm = find_corner(Node_list, cuboid[3], size_box);
    *Corner_listXpYmZm = find_corner(Node_list, cuboid[1], size_box);
    *Corner_listXpYpZm = find_corner(Node_list, cuboid[2], size_box);
    *Corner_listXmYmZp = find_corner(Node_list, cuboid[5], size_box);
    *Corner_listXmYpZp = find_corner(Node_list, cuboid[4], size_box);
    *Corner_listXpYmZp = find_corner(Node_list, cuboid[6], size_box);
    *Corner_listXpYpZp = find_corner(Node_list, cuboid[7], size_box);
    
/*    Line EdgeXmYm(Corner_listXmYmZm->coords,Corner_listXmYmZp->coords);
    Line EdgeXpYm(Corner_listXpYmZm->coords,Corner_listXpYmZp->coords);
    Line EdgeXpYp(Corner_listXpYpZm->coords,Corner_listXpYpZp->coords);
    Line EdgeXmYp(Corner_listXmYpZm->coords,Corner_listXmYpZp->coords);
    Line EdgeXmZm(Corner_listXmYmZm->coords,Corner_listXmYpZm->coords);
    Line EdgeXpZm(Corner_listXpYmZm->coords,Corner_listXpYpZm->coords);
    Line EdgeXpZp(Corner_listXpYmZp->coords,Corner_listXpYpZp->coords);
    Line EdgeXmZp(Corner_listXmYmZp->coords,Corner_listXmYpZp->coords);
    Line EdgeYmZm(Corner_listXmYmZm->coords,Corner_listXpYmZm->coords);
    Line EdgeYpZm(Corner_listXmYpZm->coords,Corner_listXpYpZm->coords);
    Line EdgeYpZp(Corner_listXmYpZp->coords,Corner_listXpYpZp->coords);
    Line EdgeYmZp(Corner_listXmYmZp->coords,Corner_listXpYmZp->coords);*/

    Line EdgeXmYm(cuboid[0],cuboid[5]);
    Line EdgeXpYm(cuboid[1],cuboid[6]);
    Line EdgeXpYp(cuboid[2],cuboid[7]);
    Line EdgeXmYp(cuboid[3],cuboid[4]);
    Line EdgeXmZm(cuboid[0],cuboid[3]);
    Line EdgeXpZm(cuboid[1],cuboid[2]);
    Line EdgeXpZp(cuboid[6],cuboid[7]);
    Line EdgeXmZp(cuboid[5],cuboid[4]);
    Line EdgeYmZm(cuboid[0],cuboid[1]);
    Line EdgeYpZm(cuboid[3],cuboid[2]);
    Line EdgeYpZp(cuboid[4],cuboid[7]);
    Line EdgeYmZp(cuboid[5],cuboid[6]);
    
    Edges.push_back(EdgeXmYm);
    Edges.push_back(EdgeXpYm);
    Edges.push_back(EdgeXpYp);
    Edges.push_back(EdgeXmYp);
    Edges.push_back(EdgeXmZm);
    Edges.push_back(EdgeXpZm);
    Edges.push_back(EdgeXpZp);
    Edges.push_back(EdgeXmZp);
    Edges.push_back(EdgeYmZm);
    Edges.push_back(EdgeYpZm);
    Edges.push_back(EdgeYpZp);
    Edges.push_back(EdgeYmZp);
    
    *Edge_listXmYm = find_edge(Node_list, EdgeXmYm, size_box);
    *Edge_listXpYm = find_edge(Node_list, EdgeXpYm, size_box);
    *Edge_listXpYp = find_edge(Node_list, EdgeXpYp, size_box);
    *Edge_listXmYp = find_edge(Node_list, EdgeXmYp, size_box);
    *Edge_listXmZm = find_edge(Node_list, EdgeXmZm, size_box);
    *Edge_listXpZm = find_edge(Node_list, EdgeXpZm, size_box);
    *Edge_listXpZp = find_edge(Node_list, EdgeXpZp, size_box);
    *Edge_listXmZp = find_edge(Node_list, EdgeXmZp, size_box);
    *Edge_listYmZm = find_edge(Node_list, EdgeYmZm, size_box);
    *Edge_listYpZm = find_edge(Node_list, EdgeYpZm, size_box);
    *Edge_listYpZp = find_edge(Node_list, EdgeYpZp, size_box);
    *Edge_listYmZp = find_edge(Node_list, EdgeYmZp, size_box);

/*    Plane FaceXm(Corner_listXmYmZm->coords, Corner_listXmYmZp->coords, Corner_listXmYpZm->coords);
    Plane FaceXp(Corner_listXpYmZm->coords, Corner_listXpYmZp->coords, Corner_listXpYpZm->coords);
    Plane FaceYm(Corner_listXmYmZm->coords, Corner_listXmYmZp->coords, Corner_listXpYmZm->coords);
    Plane FaceYp(Corner_listXmYpZm->coords, Corner_listXmYpZp->coords, Corner_listXpYpZm->coords);
    Plane FaceZm(Corner_listXmYmZm->coords, Corner_listXpYmZm->coords, Corner_listXmYpZm->coords);
    Plane FaceZp(Corner_listXmYmZp->coords, Corner_listXpYmZp->coords, Corner_listXpYpZp->coords);*/
    
    Plane FaceXm(cuboid[0], cuboid[3], cuboid[5]);
    Plane FaceXp(cuboid[1], cuboid[2], cuboid[6]);
    Plane FaceYm(cuboid[0], cuboid[1], cuboid[5]);
    Plane FaceYp(cuboid[3], cuboid[2], cuboid[4]);
    Plane FaceZm(cuboid[0], cuboid[1], cuboid[3]);
    Plane FaceZp(cuboid[5], cuboid[6], cuboid[4]);
        
    Faces.push_back(FaceXm);
    Faces.push_back(FaceXp);
    Faces.push_back(FaceYm);
    Faces.push_back(FaceYp);
    Faces.push_back(FaceZm);
    Faces.push_back(FaceZp);
            
    *Face_listXm = find_face(Node_list, FaceXm, size_box);
    *Face_listXp = find_face(Node_list, FaceXp, size_box);
    *Face_listYm = find_face(Node_list, FaceYm, size_box);
    *Face_listYp = find_face(Node_list, FaceYp, size_box);
    *Face_listZm = find_face(Node_list, FaceZm, size_box);
    *Face_listZp = find_face(Node_list, FaceZp, size_box);
    
    //Now exclude the edges of the faces
    //Xm
    vector_difference_nodes(*Face_listXm, *Edge_listXmYm);
    vector_difference_nodes(*Face_listXm, *Edge_listXmYp);
    vector_difference_nodes(*Face_listXm, *Edge_listXmZm);
    vector_difference_nodes(*Face_listXm, *Edge_listXmZp);
    
    //Xp
    vector_difference_nodes(*Face_listXp, *Edge_listXpYm);
    vector_difference_nodes(*Face_listXp, *Edge_listXpYp);
    vector_difference_nodes(*Face_listXp, *Edge_listXpZm);
    vector_difference_nodes(*Face_listXp, *Edge_listXpZp);
    
    //Ym
    vector_difference_nodes(*Face_listYm, *Edge_listXmYm);
    vector_difference_nodes(*Face_listYm, *Edge_listXpYm);
    vector_difference_nodes(*Face_listYm, *Edge_listYmZm);
    vector_difference_nodes(*Face_listYm, *Edge_listYmZp);

    //Yp
    vector_difference_nodes(*Face_listYp, *Edge_listXmYp);
    vector_difference_nodes(*Face_listYp, *Edge_listXpYp);
    vector_difference_nodes(*Face_listYp, *Edge_listYpZm);
    vector_difference_nodes(*Face_listYp, *Edge_listYpZp);

    //Zm
    vector_difference_nodes(*Face_listZm, *Edge_listXmZm);
    vector_difference_nodes(*Face_listZm, *Edge_listXpZm);
    vector_difference_nodes(*Face_listZm, *Edge_listYmZm);
    vector_difference_nodes(*Face_listZm, *Edge_listYpZm);
    
    //Zp
    vector_difference_nodes(*Face_listZp, *Edge_listXmZp);
    vector_difference_nodes(*Face_listZp, *Edge_listXpZp);
    vector_difference_nodes(*Face_listZp, *Edge_listYmZp);
    vector_difference_nodes(*Face_listZp, *Edge_listYpZp);
    
    //Now exclude the corners from the edges
    //XmYm
    vector_erase_node(*Edge_listXmYm, *Corner_listXmYmZm);
    vector_erase_node(*Edge_listXmYm, *Corner_listXmYmZp);

    //XmYp
    vector_erase_node(*Edge_listXmYp, *Corner_listXmYpZm);
    vector_erase_node(*Edge_listXmYp, *Corner_listXmYpZp);

    //XpYm
    vector_erase_node(*Edge_listXpYm, *Corner_listXpYmZm);
    vector_erase_node(*Edge_listXpYm, *Corner_listXpYmZp);

    //XpYp
    vector_erase_node(*Edge_listXpYp, *Corner_listXpYpZm);
    vector_erase_node(*Edge_listXpYp, *Corner_listXpYpZp);
    
    //XmZm
    vector_erase_node(*Edge_listXmZm, *Corner_listXmYmZm);
    vector_erase_node(*Edge_listXmZm, *Corner_listXmYpZm);

    //XmZp
    vector_erase_node(*Edge_listXmZp, *Corner_listXmYmZp);
    vector_erase_node(*Edge_listXmZp, *Corner_listXmYpZp);

    //XpZm
    vector_erase_node(*Edge_listXpZm, *Corner_listXpYmZm);
    vector_erase_node(*Edge_listXpZm, *Corner_listXpYpZm);
    
    //XpZp
    vector_erase_node(*Edge_listXpZp, *Corner_listXpYmZp);
    vector_erase_node(*Edge_listXpZp, *Corner_listXpYpZp);
    
    //YmZm
    vector_erase_node(*Edge_listYmZm, *Corner_listXmYmZm);
    vector_erase_node(*Edge_listYmZm, *Corner_listXpYmZm);
    
    //YmZp
    vector_erase_node(*Edge_listYmZp, *Corner_listXmYmZp);
    vector_erase_node(*Edge_listYmZp, *Corner_listXpYmZp);
    
    //YpZm
    vector_erase_node(*Edge_listYpZm, *Corner_listXmYpZm);
    vector_erase_node(*Edge_listYpZm, *Corner_listXpYpZm);
    
    //YpZp
    vector_erase_node(*Edge_listYpZp, *Corner_listXmYpZp);
    vector_erase_node(*Edge_listYpZp, *Corner_listXpYpZp);
}
    
//-------------------------------------------------------------
void cubic_mesh::find_pairs(const double &min_dist, const double &dist_replace)
//-------------------------------------------------------------
{
    bool perio_test = true;
    is_perio = true;
    cout << "faces Xm Xp" << endl;
    perio_test = find_face_pair(*Face_listXm, *Face_listXp, Faces[0], Faces[1], size_box, min_dist, dist_replace);
    if (perio_test == false)
        is_perio = false;
    cout << "faces Ym Yp" << endl;
    perio_test = find_face_pair(*Face_listYm, *Face_listYp, Faces[2], Faces[3], size_box, min_dist, dist_replace);
    if (perio_test == false)
        is_perio = false;
    cout << "faces Zm Zp" << endl;
    perio_test = find_face_pair(*Face_listZm, *Face_listZp, Faces[4], Faces[5], size_box, min_dist, dist_replace);
   if (perio_test == false)
       is_perio = false;
    cout << "edges XY" << endl;
    perio_test = find_edge_pair(*Edge_listXmYm, *Edge_listXpYm, *Edge_listXpYp, *Edge_listXmYp, Edges[0], Edges[1], Edges[2], Edges[3], size_box, min_dist, dist_replace);
    if (perio_test == false)
        is_perio = false;
    cout << "edges XZ" << endl;
    perio_test = find_edge_pair(*Edge_listXmZm, *Edge_listXpZm, *Edge_listXpZp, *Edge_listXmZp, Edges[4], Edges[5], Edges[6], Edges[7], size_box, min_dist, dist_replace);
    if (perio_test == false)
        is_perio = false;
    cout << "edges YZ" << endl;
    perio_test = find_edge_pair(*Edge_listYmZm, *Edge_listYpZm, *Edge_listYpZp, *Edge_listYmZp, Edges[8], Edges[9], Edges[10], Edges[11], size_box, min_dist, dist_replace);
    if (perio_test == false)
        is_perio = false;
}

//-------------------------------------------------------------
void cubic_mesh::construct_lists()
//-------------------------------------------------------------
{
    if (Corner_listXmYmZm->number != 0)
        list_of_corners.push_back(Corner_listXmYmZm);
    if (Corner_listXpYmZm->number != 0)
        list_of_corners.push_back(Corner_listXpYmZm);
    if (Corner_listXpYpZm->number != 0)
        list_of_corners.push_back(Corner_listXpYpZm);
    if (Corner_listXmYpZm->number != 0)
        list_of_corners.push_back(Corner_listXmYpZm);
    if (Corner_listXmYmZp->number != 0)
        list_of_corners.push_back(Corner_listXmYmZp);
    if (Corner_listXpYmZp->number != 0)
        list_of_corners.push_back(Corner_listXpYmZp);
    if (Corner_listXpYpZp->number != 0)
        list_of_corners.push_back(Corner_listXpYpZp);
    if (Corner_listXmYpZp->number != 0)
        list_of_corners.push_back(Corner_listXmYpZp);

    if (Edge_listXmYm->size() > 0)
        list_of_edges.push_back(Edge_listXmYm);
    if (Edge_listXpYm->size() > 0)
        list_of_edges.push_back(Edge_listXpYm);
    if (Edge_listXpYp->size() > 0)
        list_of_edges.push_back(Edge_listXpYp);
    if (Edge_listXmYp->size() > 0)
        list_of_edges.push_back(Edge_listXmYp);
    if (Edge_listXmZm->size() > 0)
        list_of_edges.push_back(Edge_listXmZm);
    if (Edge_listXpZm->size() > 0)
        list_of_edges.push_back(Edge_listXpZm);
    if (Edge_listXpZp->size() > 0)
        list_of_edges.push_back(Edge_listXpZp);
    if (Edge_listXmZp->size() > 0)
        list_of_edges.push_back(Edge_listXmZp);
    if (Edge_listYmZm->size() > 0)
        list_of_edges.push_back(Edge_listYmZm);
    if (Edge_listYpZm->size() > 0)
        list_of_edges.push_back(Edge_listYpZm);
    if (Edge_listYpZp->size() > 0)
        list_of_edges.push_back(Edge_listYpZp);
    if (Edge_listYmZp->size() > 0)
        list_of_edges.push_back(Edge_listYmZp);
    
    if (Face_listXm->size() > 0)
        list_of_faces.push_back(Face_listXm);
    if (Face_listXp->size() > 0)
        list_of_faces.push_back(Face_listXp);
    if (Face_listYm->size() > 0)
        list_of_faces.push_back(Face_listYm);
    if (Face_listYp->size() > 0)
        list_of_faces.push_back(Face_listYp);
    if (Face_listZm->size() > 0)
        list_of_faces.push_back(Face_listZm);
    if (Face_listZp->size() > 0)
        list_of_faces.push_back(Face_listZp);
}
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const cubic_mesh& cm)
//--------------------------------------------------------------------------
{
    s << "Display info on the cubic mesh:\n";
    s << "Node list file name: " << cm.Node_list_name << "\n";
    
    s << "Bounding box: " << cm.cuboid << std::endl;
    s << "volume box: " << cm.volume << std::endl;
    s << "size_box box: " << cm.size_box << std::endl;

    s << "Corner list: \n";
    for (auto n : cm.list_of_corners) {
        s << n->number << ", ";
    }
    
    s << "\n\nEdges list: ";
    s << "\nEdge_listXmYm:\t";
    for (auto n : *cm.Edge_listXmYm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXpYm:\t";
    for (auto n : *cm.Edge_listXpYm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXpYp:\t";
    for (auto n : *cm.Edge_listXpYp) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXmYp:\t";
    for (auto n : *cm.Edge_listXmYp) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXmZm:\t";
    for (auto n : *cm.Edge_listXmZm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXpZm:\t";
    for (auto n : *cm.Edge_listXpZm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXpZp:\t";
    for (auto n : *cm.Edge_listXpZp) {
        s << n.number << ", ";
    }
    s << "\nEdge_listXmZp:\t";
    for (auto n : *cm.Edge_listXmZp) {
        s << n.number << ", ";
    }
    s << "\nEdge_listYmZm:\t";
    for (auto n : *cm.Edge_listYmZm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listYpZm:\t";
    for (auto n : *cm.Edge_listYpZm) {
        s << n.number << ", ";
    }
    s << "\nEdge_listYpZp:\t";
    for (auto n : *cm.Edge_listYpZp) {
        s << n.number << ", ";
    }
    s << "\nEdge_listYmZp:\t";
    for (auto n : *cm.Edge_listYmZp) {
        s << n.number << ", ";
    }
    
    s << "\n\nfaces: \n";
    s << "\nFace_listXm:\t";
    for (auto n : *cm.Face_listXm) {
        s << n.number << ", ";
    }
    s << "\nFace_listXp:\t";
    for (auto n : *cm.Face_listXp) {
        s << n.number << ", ";
    }
    s << "\nFace_listYm:\t";
    for (auto n : *cm.Face_listYm) {
        s << n.number << ", ";
    }
    s << "\nFace_listYp:\t";
    for (auto n : *cm.Face_listYp) {
        s << n.number << ", ";
    }
    s << "\nFace_listZm:\t";
    for (auto n : *cm.Face_listZm) {
        s << n.number << ", ";
    }
    s << "\nFace_listZp:\t";
    for (auto n : *cm.Face_listZp) {
        s << n.number << ", ";
    }
    return s;
}
    
} //namespace simcoon
