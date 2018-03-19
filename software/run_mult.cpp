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

///@file Identification_LM.cpp
///@brief main: to identify // SMA model based on several uniaxial tests // Micromechanical parameters in a multiscale model of a composite material
///@author Yves Chemisky, Fodil Meraghni, Boris Piotrowski, Nicolas Despringre
///@version 0.9
///@date 01-18-2016

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/script.hpp>

#define stationnarity 1.E-12		/// Stationnary stopping criteria (no more evolution of the cost function)

using namespace std;
using namespace arma;
using namespace simcoon;

int main() {
    
	///Allow non-repetitive pseudo-random number generation
	srand(time(0));
    
	int TOOL = 1;   ///Which code is going to compute numerical files
    ofstream result;    ///Output stream, with parameters values and cost function
    
    int n_param = 0;
    int n_consts=2;
    int nfiles = 4; //number of files for the identification
    int id0 = 0;
    
    vector<parameters> params(n_param);  //vector of parameters
    vector<constants> consts(n_consts);  //vector of constants

    //Read the parameters
    read_constants(n_consts, consts, nfiles);
    
    string path_data = "data/";
    string path_keys = "data/key/";
    //Copy of the parameters to the keys folder

    copy_constants(consts, path_data, path_keys);
    generation genrun(nfiles,n_param, id0);

    string data_num_name = "NUM";
    string data_num_folder = "num_data";
    string data_num_ext = "dat";
    
    launch_solver(genrun, nfiles, params, consts, data_num_folder, data_num_name, data_num_ext);


    return 0;
}
