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

///@file Tphase_characteristics.cpp
///@brief Test for Phase characteristics objects
///@version 1.0

#include <gtest/gtest.h>
#include <fstream>
#include <iterator>
#include <armadillo>


#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read_json.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tphase_characteristics, read_write)
{
    string umat_name;
    string inputfile;
    string outputfile;
    string path_data = "data";
    string path_inputfile;
    string path_outputfile;

    vec props = {2,0};

    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;

    //Phases
    phase_characteristics rve_phase;

    rve_phase.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_phase.construct(0,1); //The rve is supposed to be mechanical only here
    inputfile = "phases" + to_string(int(rve_phase.sptr_matprops->props(1))) + ".json";
    outputfile = "phases1.json";

    read_phase_json(rve_phase, path_data, inputfile);
    write_phase_json(rve_phase, path_data, outputfile);
    path_inputfile = path_data + "/" + inputfile;
    path_outputfile = path_data + "/" + outputfile;

    std::ifstream ifs1_phase(path_inputfile);
    std::ifstream ifs2_phase(path_outputfile);

    std::istream_iterator<char> b1_phase(ifs1_phase), e1_phase;
    std::istream_iterator<char> b2_phase(ifs2_phase), e2_phase;

    ASSERT_TRUE(std::equal(b1_phase, e1_phase, b2_phase));

    //Layers
    phase_characteristics rve_layer;
    rve_layer.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_layer.construct(1,1); //The rve is supposed to be mechanical only here

    inputfile = "layers" + to_string(int(rve_layer.sptr_matprops->props(1))) + ".json";
    outputfile = "layers1.json";

    read_layer_json(rve_layer, path_data, inputfile);
    write_layer_json(rve_layer, path_data, outputfile);

    path_inputfile = path_data + "/" + inputfile;
    path_outputfile = path_data + "/" + outputfile;

    std::ifstream ifs1_layer(path_inputfile);
    std::ifstream ifs2_layer(path_outputfile);

    std::istream_iterator<char> b1_layer(ifs1_layer), e1_layer;
    std::istream_iterator<char> b2_layer(ifs2_layer), e2_layer;

    ASSERT_TRUE(std::equal(b1_layer, e1_layer, b2_layer));

    //Ellipsoid
    phase_characteristics rve_ellipsoid;
    rve_ellipsoid.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_ellipsoid.construct(2,1); //The rve is supposed to be mechanical only here

    inputfile = "ellipsoids" + to_string(int(rve_ellipsoid.sptr_matprops->props(1))) + ".json";
    outputfile = "ellipsoids1.json";

    read_ellipsoid_json(rve_ellipsoid, path_data, inputfile);
    write_ellipsoid_json(rve_ellipsoid, path_data, outputfile);

    path_inputfile = path_data + "/" + inputfile;
    path_outputfile = path_data + "/" + outputfile;

    std::ifstream ifs1_ellipsoid(path_inputfile);
    std::ifstream ifs2_ellipsoid(path_outputfile);

    std::istream_iterator<char> b1_ellipsoid(ifs1_ellipsoid), e1_ellipsoid;
    std::istream_iterator<char> b2_ellipsoid(ifs2_ellipsoid), e2_ellipsoid;

    ASSERT_TRUE(std::equal(b1_ellipsoid, e1_ellipsoid, b2_ellipsoid));

    //Cylinder
    phase_characteristics rve_cylinder;
    rve_cylinder.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_cylinder.construct(3,1); //The rve is supposed to be mechanical only here

    inputfile = "cylinders" + to_string(int(rve_cylinder.sptr_matprops->props(1))) + ".json";
    outputfile = "cylinders1.json";

    read_cylinder_json(rve_cylinder, path_data, inputfile);
    write_cylinder_json(rve_cylinder, path_data, outputfile);

    path_inputfile = path_data + "/" + inputfile;
    path_outputfile = path_data + "/" + outputfile;

    std::ifstream ifs1_cylinder(path_inputfile);
    std::ifstream ifs2_cylinder(path_outputfile);

    std::istream_iterator<char> b1_cylinder(ifs1_cylinder), e1_cylinder;
    std::istream_iterator<char> b2_cylinder(ifs2_cylinder), e2_cylinder;

    ASSERT_TRUE(std::equal(b1_cylinder, e1_cylinder, b2_cylinder));

}
