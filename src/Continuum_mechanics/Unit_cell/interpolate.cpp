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

///@file interpolate.cpp
///@brief To interpolate 1d or 2d fields
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <assert.h>
#include <armadillo>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


using namespace std;
using namespace arma;

namespace simcoon{

void find_neighbours_set(vector<vector<Node> > &neighbours, vector<mat> &neighbours_dist, const vector<Node> &inpoints, const vector<Node> &interpoints, const unsigned int &N){

    unsigned int n_inter = interpoints.size();
    double min_dist = 1.E-6;
    
    std::list<Point> points;
    for(auto n:inpoints) {
        points.push_back(n.coords);
    }
    neighbours.resize(n_inter);
    neighbours_dist.resize(n_inter);
    for (unsigned int i=0; i<n_inter; i++) {
        neighbours[i].resize(N+2);
        neighbours_dist[i] = zeros(N+1,4);
    }
    
    Tree tree(points.begin(), points.end());
        // Initialize the search structure, and search all N points
    for (unsigned int i=0; i<n_inter; i++) {
        
        Point query = interpoints[i].coords;
        Neighbor_search search(tree, query, N+1);
        // report the N nearest neighbors and their distance
        // This should sort all N points by increasing distance from origin
        neighbours[i][0] = interpoints[i];
        unsigned int j=1;
        for(Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
            for(auto n : inpoints) {
                if (squared_distance(n.coords, it->first) < min_dist) {
                    neighbours[i][j] = n;
                    if (it->second > min_dist) {
                        neighbours_dist[i](j-1,0) = sqrt(it->second);
                    }
                    else
                        neighbours_dist[i](j-1,0) = 0.;
                }
            }
            j++;
        }
        
        for(unsigned int j=0; j<N+1; j++) {
            neighbours_dist[i](j,1) = fabs(CGAL::to_double(neighbours[i][j+1].coords.x()) - CGAL::to_double(neighbours[i][0].coords.x()));
            neighbours_dist[i](j,2) = fabs(CGAL::to_double(neighbours[i][j+1].coords.y()) - CGAL::to_double(neighbours[i][0].coords.y()));
            neighbours_dist[i](j,3) = fabs(CGAL::to_double(neighbours[i][j+1].coords.z()) - CGAL::to_double(neighbours[i][0].coords.z()));
        }
        
        if((N==2)||(neighbours_dist[i](N,0) - neighbours_dist[i](N-1,0) > min_dist)){
            neighbours[i].pop_back();
            neighbours_dist[i].shed_row(N);
        }
        
        for(auto n:neighbours[i]) {
            cout << "neighbours[" << i << "] = \n" << n << endl;
        }
    }
}

void set_weights(vector<mat> &weight, const vector<mat> &neigh_dist, const unsigned int &N_in, const double &min_dis) {
   
    double sum_weight = 0.;
    double p = 1.;
    //define the weight
    for(unsigned int i=0;i<N_in; i++) {
        for(unsigned int k=0;k<4; k++) {
            
            sum_weight = 0.;
            if(neigh_dist[i](0,0) < min_dis) {
                weight[i](0,k) = 1.;
                for (unsigned int j=1; j<neigh_dist[i].n_rows; j++) {
                    weight[i](j,k) = 0.;
                }
            }
            else {
                for (unsigned int j=0; j<neigh_dist[i].n_rows; j++) {
                    if (neigh_dist[i](j,k) > min_dis) {
                        weight[i](j,k) = 1./pow(neigh_dist[i](j,k),p);
                    }
                    else {
                        weight[i](j,k) = 1.;
                    }
                    sum_weight += weight[i](j,k);
                }
                for (unsigned int j=0; j<neigh_dist[i].n_rows; j++) {
                    weight[i](j,k) = weight[i](j,k)/sum_weight;
                }
            }
        }
    }

}
                                       
                                       
equation set_equation(const std::vector<Node> &neigh, const arma::mat &weight, const unsigned int &face, const unsigned int &dof) {

    Mat<int> face_dof = {{0,2,3},{1,0,3},{1,2,0}};
    equation eq;
    eq.components.resize(neigh.size());
    for (unsigned int j=0; j<neigh.size(); j++) {
        eq.components[j].node = neigh[j];
        eq.components[j].dof = dof;
        if (j==0) {
            eq.components[j].coef = 1.0;
        }
        else {
            if (face == 0) {
                eq.components[j].coef = -1.0*weight(j-1,0);
            }
            else {
                eq.components[j].coef = -1.0*weight(j-1,face_dof(face-1,dof-1));
            }
        }
    }
    return eq;
}
    
} //namespace simcoon
