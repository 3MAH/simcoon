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

///@file Tmaterials.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "aba_material"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iterator>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( aba_write_materials )
{
    string path_data = "data";
    string inputfile_1 = "Nmat_1.inp";
    string inputfile_2 = "Nmat_2.inp";
    string inputfile_3 = "Nmat_3.inp";
    string inputfile_4 = "Nmat_4.inp";
    
    //Phases
    phase_characteristics rve_ellipsoid;
    string umat_name_macro = "MIMTN";
    vec props_rve = {2,0};
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    double density = 5.;
    double conductivity = 12.;
    
    rve_ellipsoid.sptr_matprops->update(0, umat_name_macro, 1, psi_rve, theta_rve, phi_rve, props_rve.n_elem, props_rve);
    rve_ellipsoid.construct(2,1); //The rve is supposed to be mechanical only here
    string inputfile = "Nellipsoids" + to_string(int(rve_ellipsoid.sptr_matprops->props(1))) + ".dat";
    read_ellipsoid(rve_ellipsoid, path_data, inputfile);
    
    int mid = 100000 + rve_ellipsoid.sptr_matprops->props(1)*1000 + rve_ellipsoid.sub_phases[0].sptr_matprops->number;
    int nstatev = rve_ellipsoid.sub_phases[0].sptr_sv_global->nstatev;
    int number = rve_ellipsoid.sub_phases[0].sptr_matprops->number;
    
    std::string umat_name = rve_ellipsoid.sub_phases[0].sptr_matprops->umat_name;
    int save = rve_ellipsoid.sub_phases[0].sptr_matprops->save;
    
    double psi_mat = rve_ellipsoid.sub_phases[0].sptr_matprops->psi_mat;
    double theta_mat = rve_ellipsoid.sub_phases[0].sptr_matprops->theta_mat;
    double phi_mat = rve_ellipsoid.sub_phases[0].sptr_matprops->phi_mat;
    
    int nprops = rve_ellipsoid.sub_phases[0].sptr_matprops->nprops;
    vec props = rve_ellipsoid.sub_phases[0].sptr_matprops->props;
    
    //Phases
    aba_material rve_mat_1;

    //aba_material rve_mat_2(3);
    //aba_material rve_mat_21(3,true, 200);
    
    aba_material rve_mat_3(number, mid, umat_name, save, psi_mat, theta_mat, phi_mat, density, conductivity, nprops, nstatev, props);
    
    aba_material rve_mat_4(rve_mat_3);
    rve_mat_1.update(*rve_ellipsoid.sub_phases[0].sptr_matprops, density, conductivity, *rve_ellipsoid.sub_phases[0].sptr_sv_global, mid);

    aba_material rve_mat_2 = rve_mat_1;
    
    unsigned int loading_type = 1;
    
    rve_mat_1.write(loading_type, path_data, inputfile_1);
    rve_mat_2.write(loading_type, path_data, inputfile_2);
    rve_mat_3.write(loading_type, path_data, inputfile_3);
    rve_mat_4.write(loading_type, path_data, inputfile_4);
    
    string path_inputfile_1 = path_data + "/" + inputfile_1;
    string path_inputfile_2 = path_data + "/" + inputfile_2;
    string path_inputfile_3 = path_data + "/" + inputfile_3;
    string path_inputfile_4 = path_data + "/" + inputfile_4;
    
    std::ifstream ifs1_phase(path_inputfile_1);
    std::ifstream ifs11_phase(path_inputfile_1);
    std::ifstream ifs111_phase(path_inputfile_1);
    std::ifstream ifs2_phase(path_inputfile_2);
    std::ifstream ifs3_phase(path_inputfile_3);
    std::ifstream ifs4_phase(path_inputfile_4);
    
    std::istream_iterator<char> b1_phase(ifs1_phase), e1_phase;
    std::istream_iterator<char> b11_phase(ifs11_phase), e11_phase;
    std::istream_iterator<char> b111_phase(ifs111_phase), e111_phase;
    std::istream_iterator<char> b2_phase(ifs2_phase), e2_phase;
    std::istream_iterator<char> b3_phase(ifs3_phase), e3_phase;
    std::istream_iterator<char> b4_phase(ifs4_phase), e4_phase;

    std::vector<aba_material> aba_mats;
    aba_mats.push_back(rve_mat_1);
    aba_mats.push_back(rve_mat_2);
    aba_mats.push_back(rve_mat_3);
    aba_mats.push_back(rve_mat_4);
    update_materials(aba_mats, 2, 0, 2, path_data);
    write_materials(aba_mats, loading_type, path_data, "Nmat_0.inp");
    
    section_characteristics section_rve;
    update_sections(section_rve, 2, 0, 2, path_data);
    write_sections(section_rve, loading_type, path_data, "Nmatsec_0.inp");
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_phase, e1_phase, b2_phase, e2_phase);
    BOOST_CHECK_EQUAL_COLLECTIONS(b11_phase, e11_phase, b3_phase, e3_phase);
    BOOST_CHECK_EQUAL_COLLECTIONS(b111_phase, e111_phase, b4_phase, e4_phase);
    
}

BOOST_AUTO_TEST_CASE( read_write_aba )
{
    string umat_name;
    string inputfile;
    string outputfile;
    string path_data = "data";
    string path_inputfile;
    string path_outputfile;

    unsigned int loading_type = 1;
    
    //Sections
    std::vector<section_characteristics> sections;
    
    inputfile = "Nsections0.dat";
    outputfile = "Nsections1.dat";
    
    read_sections(sections, loading_type, path_data, inputfile);
    write_sections(sections, loading_type, path_data, outputfile);
    path_inputfile = path_data + "/" + inputfile;
    path_outputfile = path_data + "/" + outputfile;
    
/*    std::ifstream ifs1_phase(path_inputfile);
    std::ifstream ifs2_phase(path_outputfile);
    
    std::istream_iterator<char> b1_phase(ifs1_phase), e1_phase;
    std::istream_iterator<char> b2_phase(ifs2_phase), e2_phase;
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_phase, e1_phase, b2_phase, e2_phase);*/
}

