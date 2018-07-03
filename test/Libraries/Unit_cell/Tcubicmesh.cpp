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

///@file Tcubicmesh.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "aba_cmesh"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iterator>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;


BOOST_AUTO_TEST_CASE( aba_construct_mesh )
{
    
    string path_data = "data";
    string path_run = "run_aba";
//    string pointsfile = "Points.inp";
    string pointsfile = "Twill_Weave.txt";

    string PBC_file_name = "PBC_file_name.inp";
    string CDN_file_name = "CDN_file_name.inp";
    string inputfile = path_data + "/" + pointsfile;
    string buffer;

    unsigned int nb_nodes = 0;
    unsigned int loading_type = 1;
    
    std::vector<Node> nodes;
    read_mesh(nodes, path_data, pointsfile);
    nb_nodes = nodes.size();
    cubic_mesh cm(nodes, pointsfile);
    cm.get_domain();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes);
    cm_perio.construct_lists();
    
    write_PBC(cm, nb_nodes, path_data, PBC_file_name);
    write_NonPerio_CDN(cm, cm_perio, loading_type, path_data, CDN_file_name);
}

BOOST_AUTO_TEST_CASE( aba_construct_perio_mesh ) {

    string path_data = "data_Unit_cell";
    string pointsfile = "node_nperio.txt";
    string points_out = "node_perio0.txt";

    string src_file = path_data + "/" + pointsfile;
    string dst_file = path_data + "/" + points_out;
    boost::filesystem::copy_file(src_file,dst_file,boost::filesystem::copy_option::overwrite_if_exists);
    
    string buffer;
    
    unsigned int nb_nodes = 0;
    
    std::vector<Node> nodes;
    read_mesh(nodes, path_data, pointsfile);
    nb_nodes = nodes.size();
    cubic_mesh cm(nodes, pointsfile);
    cm.get_domain();
    cm.construct_lists();
    
    cubic_mesh cm_perio = perio_RVE(cm, nb_nodes);
    cm_perio.construct_lists();
    
    append_perio_nodes(cm_perio, path_data, points_out);
    
}
