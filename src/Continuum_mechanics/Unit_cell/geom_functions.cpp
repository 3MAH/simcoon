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

///@file geom_functions.cpp
///@brief Geometric fonctions to apply Periodic Boundary Condition to an Element
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <CGAL/Vector_3.h>
#include <CGAL/intersections.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Direction_3 Direction;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Aff_transformation_3<Kernel> Aff_transformation_3;

using namespace std;
using namespace arma;

namespace simcoon{
    
std::vector<Node> copy_list_nodes(const std::vector<Node> &nodes) {

    std::vector<Node> copied_list;
    for(auto n : nodes) {
        copied_list.push_back(n);
    }
    
    if (copied_list.size() == 0) {
        cout << "There is no node on the list to copy" << endl;
    }
    return copied_list;
}
    
Node closest_node(const std::vector<Node> &nodes, const Point &point) {
    
    double min_distance = squared_distance(nodes[0].coords, point);
    double distance = 0.;
    Node clos_Node = nodes[0];
    for(auto n : nodes) {
        distance = squared_distance(n.coords, point);
        if (distance < min_distance) {
            clos_Node = n;
            min_distance = distance;
        }
    }
    return clos_Node;
}
    
void vector_difference_nodes(std::vector<Node> &set, const std::vector<Node> &dif) {
    
    std::vector<Node> v(set.size());
    std::vector<Node>::iterator it;
    
    std::vector<Node> dif_sort = dif;
    
    std::sort (set.begin(),set.end());
    std::sort (dif_sort.begin(),dif_sort.end());
    
    it=std::set_difference (set.begin(), set.end(), dif_sort.begin(), dif_sort.end(), v.begin());
    v.resize(it-v.begin());
    
    //std::cout << "The difference has " << (v.size()) << " elements:\n";
    //for (it=v.begin(); it!=v.end(); ++it)
    //    std::cout << ' ' << *it;
    //std::cout << '\n';
    
    set = v;
}
    
void vector_erase_node(std::vector<Node> &set, const Node &node) {
    
    set.erase(std::remove(set.begin(), set.end(), node), set.end());
}
    
Node find_corner(const std::vector<Node> &nodes, const Point &corner, const double &sizebox, const double &min_dist) {
    
    Point zero(0.,0.,0.);
    Node clos_Node(0, zero);
    unsigned int count = 0;
    for(auto n : nodes) {
        if (squared_distance(n.coords, corner) < min_dist/sizebox) {
            clos_Node = n;
            count++;
        }
    }
    if (count == 0) {
        cout << "There is no node close to the selected corner of coordinates : " << corner.x() << ",\t" << corner.y() << ",\t" << corner.z() << endl;
    }
    if (count > 1) {
        cout << "There are duplicate nodes close to the selected corner (" << corner.x() << ",\t" << corner.y() << ",\t" << corner.z() << ", the last one in the list has been selected" << endl;
    }
    return clos_Node;
        
}
    
std::vector<Node> find_edge(const std::vector<Node> &nodes, const Line &edge, const double &sizebox, const double &min_dist) {
    
    std::vector<Node> on_edge;
    unsigned int count = 0;
    for(auto n : nodes) {
        if (squared_distance(n.coords, edge) < min_dist/sizebox) {
            on_edge.push_back(n);
            count++;
        }
    }
    if (count == 0) {
        cout << "There is no node close to the selected edge" << endl;
    }
    return on_edge;
}

/*std::vector<Node> grid_edge(const Point &min, const Point &max, cont int &n_p, const bool &exclude) {
    
    std::vector<Node> grid;
    
    double pfactor = 0.;
    int pinc = 1;
    int z = 1;
    int pcol = 0;
    
    int n_samples = n_p;
    vec doe = zeros(n_samples);
    
    if(exclude == true) {
        ///Determination of parameters_equally_spaced
        for(int j=0; j<n_param; j++) {
            
            pcol = pow(spop,j);
            pinc=1;
            z=1;
            
            for(int i=0; i<n_samples; i++) {
                
                pfactor = (double) pinc/(spop+1);
                
                doe(i,j) = params[j].min_value + pfactor*(params[j].max_value-params[j].min_value);
                
                if (z==pcol) {
                    if (pinc==spop) {
                        pinc = 0;
                    }
                    pinc++;
                    z=0;
                }	
                z++;
            }
        }
        return doe;
    else{
        ///Determination of parameters_equally_spaced
        for(int j=0; j<n_param; j++) {
            
            pcol = pow(spop,j);
            pinc=0;
            z=1;
            
            for(int i=0; i<n_samples; i++) {
                
                pfactor = (double) pinc/(spop-1);
                
                doe(i,j) = params[j].min_value + pfactor*(params[j].max_value-params[j].min_value);
                
                if (z==pcol) {
                    if (pinc==spop-1) {
                        pinc = -1;
                    }
                    pinc++;
                    z=0;
                }
                z++;
            }
        }
        return doe;
    }
}*/
    
bool find_edge_pair(std::vector<Node> &edge1, std::vector<Node> &edge2, std::vector<Node> &edge3, std::vector<Node> &edge4, const double &sizebox, const double &min_dist, const double &dist_replace) {
    
    Col<int> len(4);
    len(0) = edge1.size();
    len(1) = edge2.size();
    len(2) = edge3.size();
    len(3) = edge4.size();
        
    Line edge(edge1[0].coords, edge1.back().coords);
        
    Line edge2_inter(edge2[0].coords, edge2.back().coords);
    Line edge3_inter(edge2[0].coords, edge2.back().coords);
    Line edge4_inter(edge2[0].coords, edge2.back().coords);
        
    if(abs(max(len) - min(len)) != 0)
        return false;

    std::vector<Node> edge2_temp;
    std::vector<Node> edge3_temp;
    std::vector<Node> edge4_temp;
    unsigned int count=0;
    
    for(unsigned int i=0; i<edge1.size(); i++) {
        Plane temp = edge.perpendicular_plane(edge1[i].coords);
        
        count = 0;
        for(auto n : edge2) {
            if (squared_distance(n.coords, temp) < min_dist/sizebox) {
                edge2_temp.push_back(n);
                count++;
            }
        }
        if (count > 1 ) {
            cout << "There are duplicate nodes on the edge";
            return false;
        }
        if (count == 0) {
            auto result = intersection(temp, edge2_inter);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    Node temp_point = closest_node(edge2, *p);
                    if (squared_distance(temp_point.coords, *p) < dist_replace) {
                        temp_point.coords = *p;
                        edge2_temp.push_back(temp_point);
                    }
                    else {
                        cout << "The closest point is too far for a replace" << endl;
                        return false;
                    }
                }
                else {
                    cout << "The intersection between edge and constructing plane is not a Point" << endl;
                    return false;

                }
            }
            else {
                cout << "There are no pair node on this part of the edge, and intersection between edge and constructing plane is null" << endl;
                return false;
            }
        }
        
        count = 0;
        for(auto n : edge3) {
            if (squared_distance(n.coords, temp) < min_dist/sizebox) {
                edge3_temp.push_back(n);
                count++;
            }
        }
        if (count > 1 ) {
            cout << "There are duplicate nodeson the edge";
            return false;
        }
        if (count == 0) {
            auto result = intersection(temp, edge3_inter);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    Node temp_point = closest_node(edge3, *p);
                    if (squared_distance(temp_point.coords, *p) < dist_replace) {
                        temp_point.coords = *p;
                        edge3_temp.push_back(temp_point);
                    }
                    else {
                        cout << "The closest point is too far for a replace" << endl;
                        return false;
                    }
                }
                else {
                    cout << "The intersection between edge and constructing plane is not a Point" << endl;
                    return false;
                    
                }
            }
            else {
                cout << "There are no pair node on this part of the edge, and intersection between edge and constructing plane is null" << endl;
                return false;
            }
        }
        
        count = 0;
        for(auto n : edge4) {
            if (squared_distance(n.coords, temp) < min_dist/sizebox) {
                edge4_temp.push_back(n);
                count++;
            }
        }
        if (count > 1 ) {
            cout << "There are duplicate nodes on the edge";
            return false;
        }
        if (count == 0) {
            auto result = intersection(temp, edge4_inter);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    Node temp_point = closest_node(edge4, *p);
                    if (squared_distance(temp_point.coords, *p) < dist_replace) {
                        temp_point.coords = *p;
                        edge4_temp.push_back(temp_point);
                    }
                    else {
                        cout << "The closest point is too far for a replace" << endl;
                        return false;
                    }
                }
                else {
                    cout << "The intersection between edge and constructing plane is not a Point" << endl;
                    return false;
                    
                }
            }
            else {
                cout << "There are no pair node on this part of the edge, and intersection between edge and constructing plane is null" << endl;
                return false;
            }
        }
        
    }

    edge2 = edge2_temp;
    edge3 = edge3_temp;
    edge4 = edge4_temp;
    return true;
}
    
std::vector<Node> find_face(const std::vector<Node> &nodes, const Plane &face, const double &sizebox, const double &min_dist) {
    
    std::vector<Node> on_face;
    unsigned int count = 0;
    for(auto n : nodes) {
        if (squared_distance(n.coords, face) < min_dist/sizebox) {
            on_face.push_back(n);
            count++;
        }
    }
    if (count == 0) {
        cout << "There is no node close to the selected face" << endl;
    }
    return on_face;
}

bool find_face_pair(std::vector<Node> &face1, std::vector<Node> &face2, const double &sizebox, const double &min_dist, const double &dist_replace)
{
    Plane face(face1[0].coords, face1[1].coords, face1.back().coords);
    Plane face2_inter(face2[0].coords, face2[1].coords, face2.back().coords);
    
    if(face1.size() - face2.size() != 0)
        return false;
    
    std::vector<Node> face2_temp;
    unsigned int count=0;
    
    for(unsigned int i=0; i<face1.size(); i++) {
        Line temp = face.perpendicular_line(face1[i].coords);
        count = 0;
        for(auto n : face2) {
            if (squared_distance(n.coords, temp) < min_dist/sizebox) {
                face2_temp.push_back(n);
                count++;
            }
        }
        if (count > 1 ) {
            cout << "There are duplicate nodes on the face";
            return false;
        }
        if (count == 0) {
            auto result = intersection(temp, face2_inter);
            if (result) {
                if (const Point *p = boost::get<Point>(&*result)) {
                    Node temp_point = closest_node(face2, *p);
                    if (squared_distance(temp_point.coords, *p) < dist_replace) {
                        temp_point.coords = *p;
                        face2_temp.push_back(temp_point);
                    }
                    else {
                        cout << "The closest point is too far for a replace" << endl;
                        return false;
                    }
                }
                else {
                    cout << "The intersection between edge and constructing plane is not a Point" << endl;
                    return false;
                    
                }
            }
            else {
                cout << "There are no pair node on this part of the edge, and intersection between edge and constructing plane is null" << endl;
                return false;
            }
        }
    }
    
    face2 = face2_temp;
    return true;
}

Node duplicate_node(const Node &node, unsigned int &nb_nodes) {

    assert(nb_nodes > 0);
    Node duplicated = node;
    nb_nodes++;
    if (nb_nodes !=0) {
        duplicated.number = nb_nodes;
    }
    return duplicated;
}
    
std::vector<Node> duplicate_list_nodes(const std::vector<Node> &nodes, unsigned int &nb_nodes) {
    
    assert(nb_nodes > 0);
    std::vector<Node> duplicated_list;
    for(auto n : nodes) {
        duplicated_list.push_back(n);
        nb_nodes++;
        n.number = nb_nodes;
    }
    
    if (duplicated_list.size() == 0) {
        cout << "There is no node on the list to duplicate";
    }
    
    return duplicated_list;
}
    
void translate_node(Node &node, const Kernel::FT &Dp, const int &axis, unsigned int &nb_nodes) {
    
    Vector_3 trans_vec;//(0,0,0);
    
    switch (axis) {
        case 1: {
            trans_vec = Vector_3(Dp,0,0);
            break;
        }
        case 2: {
            trans_vec = Vector_3(0,Dp,0);
            break;
        }
        case 3: {
            trans_vec = Vector_3(0,0,Dp);
            break;
        }
        default: {
            cout << "Error: in translate_nodes function: Please enter a valid axis : 1, 2 or 3" << endl;
            exit(0);
        }
    }
    
    Aff_transformation_3 translate(CGAL::TRANSLATION, trans_vec);
    Point translated = node.coords.transform(translate);
    node.coords = translated;
    if (nb_nodes != 0) {
        nb_nodes++;
        node.number = nb_nodes;
    }
}
    
void translate_node(Node &node, const Vector_3 &trans_vec, unsigned int &nb_nodes) {
    
    Aff_transformation_3 translate(CGAL::TRANSLATION, trans_vec);
    Point translated = node.coords.transform(translate);
    node.coords = translated;
    if (nb_nodes != 0) {
        nb_nodes++;
        node.number = nb_nodes;
    }
}
    
void translate_nodes(std::vector<Node> &nodes, const Kernel::FT &Dp, const int &axis, unsigned int &nb_nodes) {
        
    Vector_3 trans_vec;//(0,0,0);
    
    switch (axis) {
            case 1: {
                trans_vec = Vector_3(Dp,0,0);
                break;
            }
            case 2: {
                trans_vec = Vector_3(0,Dp,0);
                break;
            }
            case 3: {
                trans_vec = Vector_3(0,0,Dp);
                break;
            }
            default: {
                cout << "Error: in translate_nodes function: Please enter a valid axis : 1, 2 or 3" << endl;
                exit(0);
            }
    }
    
    Aff_transformation_3 translate(CGAL::TRANSLATION, trans_vec);
    Point translated;
    for(unsigned int i=0; i<nodes.size(); i++) {
        translated = nodes[i].coords.transform(translate);
        nodes[i].coords = translated;
        if (nb_nodes != 0) {
            nb_nodes++;
            nodes[i].number = nb_nodes;
        }
    }
}
    
void translate_nodes(std::vector<Node> &nodes, const Vector_3 &trans_vec, unsigned int &nb_nodes) {
    
    Aff_transformation_3 translate(CGAL::TRANSLATION, trans_vec);
    Point translated;
    for(unsigned int i=0; i<nodes.size(); i++) {
        translated = nodes[i].coords.transform(translate);
        nodes[i].coords = translated;
        if (nb_nodes != 0) {
            nb_nodes++;
            nodes[i].number = nb_nodes;
        }
    }
}
    
cubic_mesh perio_RVE(cubic_mesh &RVE, unsigned int &nb_nodes) {
    
    cubic_mesh perio_mesh;
    
    perio_mesh.is_perio = RVE.is_perio;
    perio_mesh.Node_list_name = RVE.Node_list_name;
    
    perio_mesh.cuboid = RVE.cuboid;
    perio_mesh.volume = RVE.volume;
    perio_mesh.Dx = RVE.Dx;
    perio_mesh.Dy = RVE.Dy;
    perio_mesh.Dz = RVE.Dz;
    perio_mesh.size_box = RVE.size_box;
    
    perio_mesh.center_node = RVE.center_node;
    
    perio_mesh.construct();
    
    Vector_3 trans_vec;
    
    std::vector<Node> temp_list;
    
    *perio_mesh.Face_listXm = copy_list_nodes(*RVE.Face_listXm);
    *perio_mesh.Face_listYm = copy_list_nodes(*RVE.Face_listYm);
    *perio_mesh.Face_listZm = copy_list_nodes(*RVE.Face_listZm);
    
    *perio_mesh.Face_listXp = copy_list_nodes(*RVE.Face_listXm);
    trans_vec = Vector_3(perio_mesh.Dx,0.0,0.0);
    translate_nodes(*perio_mesh.Face_listXp, trans_vec, nb_nodes);

    *perio_mesh.Face_listYp = copy_list_nodes(*RVE.Face_listYm);
    trans_vec = Vector_3(0.0,perio_mesh.Dy,0.0);
    translate_nodes(*perio_mesh.Face_listYp, trans_vec, nb_nodes);

    *perio_mesh.Face_listZp = copy_list_nodes(*RVE.Face_listZm);
    trans_vec = Vector_3(0.0,0.0,perio_mesh.Dz);
    translate_nodes(*perio_mesh.Face_listZp, trans_vec, nb_nodes);

    //Edge 1 : colinear to Z
    *perio_mesh.Edge_listXmYm = copy_list_nodes(*RVE.Edge_listXmYm);

    *perio_mesh.Edge_listXpYm = copy_list_nodes(*RVE.Edge_listXmYm);
    trans_vec = Vector_3(perio_mesh.Dx,0.0,0.0);
    translate_nodes(*perio_mesh.Edge_listXpYm, trans_vec, nb_nodes);
    
    *perio_mesh.Edge_listXpYp = copy_list_nodes(*RVE.Edge_listXmYm);
    trans_vec = Vector_3(perio_mesh.Dx,perio_mesh.Dy,0.0);
    translate_nodes(*perio_mesh.Edge_listXpYp, trans_vec, nb_nodes);
    
    *perio_mesh.Edge_listXmYp = copy_list_nodes(*RVE.Edge_listXmYm);
    trans_vec = Vector_3(0.0,perio_mesh.Dy,0.0);
    translate_nodes(*perio_mesh.Edge_listXmYp, trans_vec, nb_nodes);
    
    //Edge 2 : colinear to Y
    *perio_mesh.Edge_listXmZm = copy_list_nodes(*RVE.Edge_listXmZm);
    
    *perio_mesh.Edge_listXpZm = copy_list_nodes(*RVE.Edge_listXmZm);
    trans_vec = Vector_3(perio_mesh.Dx,0.0,0.0);
    translate_nodes(*perio_mesh.Edge_listXpZm, trans_vec, nb_nodes);

    *perio_mesh.Edge_listXpZp = copy_list_nodes(*RVE.Edge_listXmZm);
    trans_vec = Vector_3(perio_mesh.Dx,0.0,perio_mesh.Dz);
    translate_nodes(*perio_mesh.Edge_listXpZp, trans_vec, nb_nodes);
    
    *perio_mesh.Edge_listXmZp = copy_list_nodes(*RVE.Edge_listXmZm);
    trans_vec = Vector_3(0.0,0.0,perio_mesh.Dz);
    translate_nodes(*perio_mesh.Edge_listXmZp, trans_vec, nb_nodes);

    //Edge 3 : colinear to Y
    *perio_mesh.Edge_listYmZm = copy_list_nodes(*RVE.Edge_listYmZm);
    
    *perio_mesh.Edge_listYpZm = copy_list_nodes(*RVE.Edge_listYmZm);
    trans_vec = Vector_3(0.0,perio_mesh.Dy,0.0);
    translate_nodes(*perio_mesh.Edge_listYpZm, trans_vec, nb_nodes);
    
    *perio_mesh.Edge_listYpZp = copy_list_nodes(*RVE.Edge_listYmZm);
    trans_vec = Vector_3(0.0,perio_mesh.Dy,perio_mesh.Dz);
    translate_nodes(*perio_mesh.Edge_listYpZp, trans_vec, nb_nodes);

    *perio_mesh.Edge_listYmZp = copy_list_nodes(*RVE.Edge_listYmZm);
    trans_vec = Vector_3(0.0,0.0,perio_mesh.Dz);
    translate_nodes(*perio_mesh.Edge_listYmZp, trans_vec, nb_nodes);
    
    //Corners
    *perio_mesh.Corner_listXmYmZm = *RVE.Corner_listXmYmZm;
    *perio_mesh.Corner_listXmYpZm = *RVE.Corner_listXmYpZm;
    *perio_mesh.Corner_listXpYmZm = *RVE.Corner_listXpYmZm;
    *perio_mesh.Corner_listXpYpZm = *RVE.Corner_listXpYpZm;
    *perio_mesh.Corner_listXmYmZp = *RVE.Corner_listXmYmZp;
    *perio_mesh.Corner_listXmYpZp = *RVE.Corner_listXmYpZp;
    *perio_mesh.Corner_listXpYmZp = *RVE.Corner_listXpYmZp;
    *perio_mesh.Corner_listXpYpZp = *RVE.Corner_listXpYpZp;

/*    *perio_mesh.Corner_listXmYmZm = *RVE.Corner_listXmYmZm;
    *perio_mesh.Corner_listXmYpZm = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(0.0,perio_mesh.Dy,0.);
    translate_node(*perio_mesh.Corner_listXmYpZm, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXpYmZm = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(perio_mesh.Dx,0.0,0.);
    translate_node(*perio_mesh.Corner_listXpYmZm, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXpYpZm = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(perio_mesh.Dx,perio_mesh.Dx,0.);
    translate_node(*perio_mesh.Corner_listXpYpZm, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXmYmZp = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(0.0,0.0,perio_mesh.Dz);
    translate_node(*perio_mesh.Corner_listXmYmZp, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXmYpZp = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(0.0,perio_mesh.Dy,perio_mesh.Dz);
    translate_node(*perio_mesh.Corner_listXmYpZp, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXpYmZp = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(perio_mesh.Dx,0.0,perio_mesh.Dz);
    translate_node(*perio_mesh.Corner_listXpYmZp, trans_vec, nb_nodes);
    *perio_mesh.Corner_listXpYpZp = *RVE.Corner_listXmYmZm;
    trans_vec = Vector_3(perio_mesh.Dx,perio_mesh.Dy,perio_mesh.Dz);
    translate_node(*perio_mesh.Corner_listXpYpZp, trans_vec, nb_nodes);*/
    
/*    *perio_mesh.Corner_listXmYpZm = duplicate_node(*RVE.Corner_listXmYpZm, nb_nodes);
    *perio_mesh.Corner_listXpYmZm = duplicate_node(*RVE.Corner_listXpYmZm, nb_nodes);
    *perio_mesh.Corner_listXpYpZm = duplicate_node(*RVE.Corner_listXpYpZm, nb_nodes);
    *perio_mesh.Corner_listXmYmZp = duplicate_node(*RVE.Corner_listXmYmZp, nb_nodes);
    *perio_mesh.Corner_listXmYpZp = duplicate_node(*RVE.Corner_listXmYpZp, nb_nodes);
    *perio_mesh.Corner_listXpYmZp = duplicate_node(*RVE.Corner_listXpYmZp, nb_nodes);
    *perio_mesh.Corner_listXpYpZp = duplicate_node(*RVE.Corner_listXpYpZp, nb_nodes);*/
    
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listXm), std::end(*perio_mesh.Face_listXm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listXp), std::end(*perio_mesh.Face_listXp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listYm), std::end(*perio_mesh.Face_listYm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listYp), std::end(*perio_mesh.Face_listYp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listZm), std::end(*perio_mesh.Face_listZm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Face_listZp), std::end(*perio_mesh.Face_listZp));

    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXmYm), std::end(*perio_mesh.Edge_listXmYm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXpYm), std::end(*perio_mesh.Edge_listXpYm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXpYp), std::end(*perio_mesh.Edge_listXpYp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXmYp), std::end(*perio_mesh.Edge_listXmYp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXmZm), std::end(*perio_mesh.Edge_listXmZm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXpZm), std::end(*perio_mesh.Edge_listXpZm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXpZp), std::end(*perio_mesh.Edge_listXpZp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listXmZp), std::end(*perio_mesh.Edge_listXmZp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listYmZm), std::end(*perio_mesh.Edge_listYmZm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listYpZm), std::end(*perio_mesh.Edge_listYpZm));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listYpZp), std::end(*perio_mesh.Edge_listYpZp));
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), std::begin(*perio_mesh.Edge_listYmZp), std::end(*perio_mesh.Edge_listYmZp));
    
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXmYmZm);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXmYpZm);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXpYmZm);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXpYpZm);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXmYmZp);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXmYpZp);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXpYmZp);
    perio_mesh.Node_list.insert(std::end(perio_mesh.Node_list), *perio_mesh.Corner_listXpYpZp);
    
    return perio_mesh;
}

void set_faceXm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listXm), std::end(*cm_perio.Face_listXm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYm), std::end(*cm_perio.Edge_listXmYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYp), std::end(*cm_perio.Edge_listXmYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZm), std::end(*cm_perio.Edge_listXmZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZp), std::end(*cm_perio.Edge_listXmZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
}
    
void set_faceYm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listYm), std::end(*cm_perio.Face_listYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYm), std::end(*cm_perio.Edge_listXmYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYm), std::end(*cm_perio.Edge_listXpYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZm), std::end(*cm_perio.Edge_listYmZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZp), std::end(*cm_perio.Edge_listYmZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
}

void set_faceZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listZm), std::end(*cm_perio.Face_listZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZm), std::end(*cm_perio.Edge_listXmZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZm), std::end(*cm_perio.Edge_listXpZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZm), std::end(*cm_perio.Edge_listYmZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZm), std::end(*cm_perio.Edge_listYpZm));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
}
    
void set_EdgeXmYm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYm), std::end(*cm_perio.Edge_listXmYm));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
}

void set_EdgeXmZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZm), std::end(*cm_perio.Edge_listXmZm));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
}

void set_EdgeYmZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZm), std::end(*cm_perio.Edge_listYmZm));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
}

void set_faceXp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listXp), std::end(*cm_perio.Face_listXp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYm), std::end(*cm_perio.Edge_listXpYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYp), std::end(*cm_perio.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZm), std::end(*cm_perio.Edge_listXpZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZp), std::end(*cm_perio.Edge_listXpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

void set_faceYp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listYp), std::end(*cm_perio.Face_listYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYp), std::end(*cm_perio.Edge_listXmYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYp), std::end(*cm_perio.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZm), std::end(*cm_perio.Edge_listYpZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZp), std::end(*cm_perio.Edge_listYpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

void set_faceZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listZp), std::end(*cm_perio.Face_listZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZp), std::end(*cm_perio.Edge_listXmZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZp), std::end(*cm_perio.Edge_listXpZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZp), std::end(*cm_perio.Edge_listYmZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZp), std::end(*cm_perio.Edge_listYpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

void set_EdgeXpYm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYm), std::end(*cm_perio.Edge_listXpYm));
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
}

void set_EdgeXpYp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYp), std::end(*cm_perio.Edge_listXpYp));
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

void set_EdgeXmYp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYp), std::end(*cm_perio.Edge_listXmYp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
}
    
void set_EdgeXpZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZm), std::end(*cm_perio.Edge_listXpZm));
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
}

void set_EdgeXpZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZp), std::end(*cm_perio.Edge_listXpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

void set_EdgeXmZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZp), std::end(*cm_perio.Edge_listXmZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
}

void set_EdgeYpZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZm), std::end(*cm_perio.Edge_listYpZm));
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
}

void set_EdgeYpZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZp), std::end(*cm_perio.Edge_listYpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}
    
void set_EdgeYmZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZp), std::end(*cm_perio.Edge_listYmZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
}

void set_CornerXmYpZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
}

void set_CornerXpYmZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
}

void set_CornerXpYpZm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
}

void set_CornerXmYmZpm(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
}

void set_CornerXmYpZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
}

void set_CornerXpYmZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
}

void set_CornerXpYpZp(std::vector<Node> &set, const cubic_mesh &cm_perio) {
    set.clear();
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
}

/*void fluctuation(equation &eq, std::vector<Node> &set_face, const cubic_mesh& cm, const cubic_equation &cubic_eq, const std::string &name_set, const unsigned int &dof) {

    double xmin = cm.cuboid.xmin();
    double ymin = cm.cuboid.ymin();
    double zmin = cm.cuboid.zmin();
    
    mat Dxyz = {{-1.*cm.Dx,-0.5*cm.Dx,-0.5*cm.Dx},{-0.5*cm.Dy,-1.*cm.Dy,-0.5*cm.Dy},{-0.5*cm.Dz,-0.5*cm.Dz,-1.*cm.Dz}};
    umat CD_num = {{0,3,4},{3,1,5},{4,5,2}};
    vec perio_disp = zeros(3);
        
    perio_disp(0) += (CGAL::to_double(n.coords.x()) - xmin)/cm.Dx;
    perio_disp(1) += (CGAL::to_double(n.coords.y()) - ymin)/cm.Dy;
    perio_disp(2) += (CGAL::to_double(n.coords.z()) - zmin)/cm.Dz;
    
    
    if(name_set == cm.set_name_faces[0]) {
        for(auto n:set_face)
        //FaceXm
        for (unsigned int i=0; i<3; i++) {
            temp_component.node = cubic_eq.CD_nodes[CD_num(0,i)];
            temp_component.dof = 1;
            temp_component.coef = perio_disp(i);
            eq.component.push_back(temp_component);
        }
    }
    else if(name_set == cm.set_name_faces[2]) {
        //FaceYp
        for (unsigned int i=0; i<3; i++) {
            temp_component.node = cubic_eq.CD_nodes[CD_num(1,i)];
            temp_component.dof = 1;
            temp_component.coef = perio_disp(i);
            eq.component.push_back(temp_component);
        }
    }
    else if(name_set == cm.set_name_faces[3]) {
        //FaceZp
        for (unsigned int i=0; i<3; i++) {
            temp_component.node = cubic_eq.CD_nodes[CD_num(2,i)];
            temp_component.dof = 1;
            temp_component.coef = perio_disp(i);
            eq.component.push_back(temp_component);
        }
    }
}*/
    
unsigned int index_from_dof(const unsigned int &dof, const unsigned int &loading_type) {
    
    UNUSED(loading_type);
    
    if (dof < 4) {
        return dof-1;
    }
    else if(dof == 11) {
        return 0;
    }
    else {
        cout << "in Continuum_mechanics/Unit/cell/geom_functions.cpp : error: dof is not recognized" << endl;
        return 0;
    }
}
    
void replace_perio_eq(equation &eq, const cubic_equation &cubic_eq, const mat &weight, const cubic_mesh &cm, const std::string &name_set, const unsigned int &loading_type, const unsigned int &dof) {
    
    UNUSED(weight);
    std::vector<equation> set_eq;
    //Iterator on the node
    std::vector<component>::iterator it;
//    double weight_temp = 0.;
    equation coefs_temps;
    unsigned int N=0;
    uvec is_unique = ones<uvec>(N);
    
    unsigned int j_dof = index_from_dof(dof, loading_type);
    
    if(name_set == cm.set_name_faces[1]) {
        //FaceXp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Face_listXp[j_dof]), std::end(cubic_eq.Face_listXp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpYm[j_dof]), std::end(cubic_eq.Edge_listXpYm[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpYp[j_dof]), std::end(cubic_eq.Edge_listXpYp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpZm[j_dof]), std::end(cubic_eq.Edge_listXpZm[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpZp[j_dof]), std::end(cubic_eq.Edge_listXpZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_faces[3]) {
        //FaceYp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Face_listYp[j_dof]), std::end(cubic_eq.Face_listYp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXmYp[j_dof]), std::end(cubic_eq.Edge_listXmYp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpYp[j_dof]), std::end(cubic_eq.Edge_listXpYp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYpZm[j_dof]), std::end(cubic_eq.Edge_listYpZm[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYpZp[j_dof]), std::end(cubic_eq.Edge_listYpZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_faces[5]) {
        //FaceZp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Face_listZp[j_dof]), std::end(cubic_eq.Face_listZp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXmZp[j_dof]), std::end(cubic_eq.Edge_listXmZp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpZp[j_dof]), std::end(cubic_eq.Edge_listXpZp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYmZp[j_dof]), std::end(cubic_eq.Edge_listYmZp[j_dof]));
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYpZp[j_dof]), std::end(cubic_eq.Edge_listYpZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[1]) {
        //EdgeXpYm
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpYm[j_dof]), std::end(cubic_eq.Edge_listXpYm[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[2]) {
        //EdgeXpYp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpYp[j_dof]), std::end(cubic_eq.Edge_listXpYp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[3]) {
        //EdgeXmYp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXmYp[j_dof]), std::end(cubic_eq.Edge_listXmYp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[5]) {
        //EdgeXpZm
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpZm[j_dof]), std::end(cubic_eq.Edge_listXpZm[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZm[j_dof]);
    }
    else if(name_set == cm.set_name_edges[6]) {
        //EdgeXpZp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXpZp[j_dof]), std::end(cubic_eq.Edge_listXpZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[7]) {
        //EdgeXmZp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listXmZp[j_dof]), std::end(cubic_eq.Edge_listXmZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[9]) {
        //EdgeYpZm
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYpZm[j_dof]), std::end(cubic_eq.Edge_listYpZm[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZm[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZm[j_dof]);
    }
    else if(name_set == cm.set_name_edges[10]) {
        //EdgeYpZp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYpZp[j_dof]), std::end(cubic_eq.Edge_listYpZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYpZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYpZp[j_dof]);
    }
    else if(name_set == cm.set_name_edges[11]) {
        //EdgeYmZp
        set_eq.clear();
        set_eq.insert(std::end(set_eq), std::begin(cubic_eq.Edge_listYmZp[j_dof]), std::end(cubic_eq.Edge_listYmZp[j_dof]));
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXmYmZp[j_dof]);
        set_eq.insert(std::end(set_eq), cubic_eq.Corner_listXpYmZp[j_dof]);
    }
    
    equation temp_eq;
    bool replace = false;
    temp_eq.components.push_back(eq.components[0]);
    for (it = eq.components.begin()+1; it != eq.components.end(); ++it) {
        for(auto n = set_eq.begin(); ((n != set_eq.end())&&(replace != true)); ++n) {
            if(it->node == n->components[0].node) {

                coefs_temps.components = n->components;
                for(unsigned int k=0; k<coefs_temps.components.size(); k++) {
                    coefs_temps.components[k].coef *= -1.0*it->coef;
                }
                replace = true;
                temp_eq.components.insert(temp_eq.components.end(), coefs_temps.components.begin()+1, coefs_temps.components.end());
            }
        }
        if(replace == false) {
            temp_eq.components.push_back(*it);
        }
        else {
            replace = false;
        }
    }
    eq = temp_eq;
    
    N = eq.components.size();
    is_unique = ones<uvec>(N);
    for (unsigned int k=0; k<N-1; k++) {
        for (unsigned int l=k+1; l<N; l++) {
            if ((eq.components[k] == eq.components[l])&&(is_unique(l) == 1)) {
                is_unique(l) = 0;
                eq.components[k].coef += eq.components[l].coef;
                eq.components[l].coef = 0.;
            }
        }
    }
    for (it = eq.components.begin(); it != eq.components.end();) {
        if (fabs(it->coef) < 1.E-6)
            it = eq.components.erase(it);
        else
            ++it;
    }
}
    
    

} //namespace simcoon
