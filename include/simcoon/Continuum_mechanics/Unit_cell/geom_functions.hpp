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

///@file geom_functions.hpp
///@brief specific geometric functions
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_equation.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector_3;

namespace simcoon{
    
// Computes the closest node of a prescribed point from a set (vector) of nodes

std::vector<Node> copy_list_nodes(const std::vector<Node> &);
    
Node closest_node(const std::vector<Node> &, const Point &);

void vector_difference_nodes(std::vector<Node> &, const std::vector<Node> &);
    
void vector_erase_node(std::vector<Node> &, const Node &);
    
Node find_corner(const std::vector<Node> &, const Point &, const double &, const double & = 1.E-6);

std::vector<Node> find_edge(const std::vector<Node> &, const Line &, const double &, const double & = 1.E-6);
    
bool find_edge_pair(std::vector<Node> &, std::vector<Node> &, std::vector<Node> &, std::vector<Node> &, const double &, const double & = 1.E-6, const double & = 1.E-4);
    
std::vector<Node> find_face(const std::vector<Node> &, const Plane &, const double &, const double & = 1.E-6);
    
bool find_face_pair(std::vector<Node> &, std::vector<Node> &, const double &, const double & = 1.E-6, const double & = 1.E-4);

Node duplicate_node(const Node &, unsigned int &);
    
std::vector<Node> duplicate_list_nodes(const std::vector<Node> &, unsigned int &);
    
void translate_node(Node &, const Kernel::FT &, const int &, unsigned int &);

void translate_node(Node &, const Vector_3 &, unsigned int &);
    
void translate_nodes(std::vector<Node> &, const Kernel::FT &, const int &, unsigned int &);
    
void translate_nodes(std::vector<Node> &, const Vector_3 &, unsigned int &);
    
cubic_mesh perio_RVE(cubic_mesh &RVE, unsigned int &);
    
void set_faceXp(std::vector<Node> &, const cubic_mesh &);
void set_faceYp(std::vector<Node> &, const cubic_mesh &);
void set_faceZp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXpYm(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXpYp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXmYp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXpZm(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXpZp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeXmZp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeYpZm(std::vector<Node> &, const cubic_mesh &);
void set_EdgeYpZp(std::vector<Node> &, const cubic_mesh &);
void set_EdgeYmZp(std::vector<Node> &, const cubic_mesh &);
    
void set_CornerXmYpZm(std::vector<Node> &, const cubic_mesh &);
void set_CornerXpYmZm(std::vector<Node> &, const cubic_mesh &);
void set_CornerXpYpZm(std::vector<Node> &, const cubic_mesh &);
void set_CornerXmYmZpm(std::vector<Node> &, const cubic_mesh &);
void set_CornerXmYpZp(std::vector<Node> &, const cubic_mesh &);
void set_CornerXpYmZp(std::vector<Node> &, const cubic_mesh &);
void set_CornerXpYpZp(std::vector<Node> &, const cubic_mesh &);
    
unsigned int index_from_dof(const unsigned int &, const unsigned int &);

void replace_perio_eq(equation &, const cubic_equation &, const arma::mat &, const cubic_mesh &, const std::string &, const unsigned int &, const unsigned int &);

std::vector<equation> MPC_equations_non_perio(const cubic_mesh &, const cubic_mesh &, const cubic_equation &, const unsigned int &, const unsigned int &);
    
} //namespace simcoon

