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

///@file read.cpp
///@brief To read parameters, variables and input values for functions
///@version 1.0

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <simcoon/Continuum_Mechanics/Functions/read.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void read_func_N(vec &params, vec &variables, string &N_file, const string &path_data, const string &inputfile) {
        
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramfunc;
    
    unsigned int nparams = 0;
    unsigned int nvariables = 0;
    
    paramfunc.open(path_inputfile, ios::in);
    if(paramfunc) {
        paramfunc >> buffer >> nparams >> buffer >> nvariables;
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " in the folder :" << path_data << endl;
    }
    paramfunc.close();
    params = zeros(nparams);
    variables = zeros(nvariables);
    
    paramfunc.open(path_inputfile, ios::in);
    paramfunc >> buffer >> buffer >> buffer >> buffer;
    paramfunc >> buffer;
    for(unsigned int i=0; i<nparams; i++) {
        paramfunc >> buffer >> params(i);
    }
    paramfunc >> buffer;
    for(unsigned int i=0; i<nvariables; i++) {
        paramfunc >> buffer >> variables(i);
    }
    
    paramfunc >> buffer >> N_file;
    paramfunc.close();
}

} //namespace simcoon
