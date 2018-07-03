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

///@file Tsteps.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "aba_step"
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iterator>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_meca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( aba_write_steps )
{
    string path_data = "data";
    string pathfile = "path.txt";
    
    ///Usefull UMAT variables
    std::vector<block> blocks;  //loading blocks
    
    double T_init = 0.;
    bool nlgeom = false;
    int loading_type = 1;
    int max_temp = 10;
    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);
    
    //Phases
    std::vector<std::shared_ptr<step> > aba_steps;
    
    update_steps(aba_steps, blocks, nlgeom, loading_type, max_temp);
    write_steps(aba_steps, loading_type, T_init, path_data, "Nstep_0.inp");
    
/*    string path_inputfile_1 = path_data + "/" + inputfile_1;
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
    write_materials(aba_mats, path_data, "Nmat_0.inp");
    
    BOOST_CHECK_EQUAL_COLLECTIONS(b1_phase, e1_phase, b2_phase, e2_phase);
    BOOST_CHECK_EQUAL_COLLECTIONS(b11_phase, e11_phase, b3_phase, e3_phase);
    BOOST_CHECK_EQUAL_COLLECTIONS(b111_phase, e111_phase, b4_phase, e4_phase);*/
    
}
