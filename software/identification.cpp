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

#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/identification.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

int main() {
    
    ofstream result;    ///Output stream, with parameters values and cost function
    
    int n_param;
    int n_consts;
    int nfiles = 0; //number of files for the identification
    
    //Parameters for the optimization software
	int ngen;
	int aleaspace;	
	int apop;	
	int spop;
	int ngboys;
	int maxpop;
    int station_nb;
    double probaMut;
    double pertu;
    double c;	///Lagrange penalty parameters
    double p0;
	double lambdaLM;
    //Read the identification control
    
    string path_data = "data";
    string path_keys = "keys";
    string path_results = "results";
    string materialfile = "material.dat";
    string outputfile = "id_params.txt";
    string simulfile = "simul.txt";
    
    string file_essentials = "ident_essentials.inp";
    string file_control = "ident_control.inp";

    string simul_type = "SOLVE";

    ident_essentials(n_param, n_consts, nfiles, path_data, file_essentials);
    ident_control(ngen, aleaspace, apop, spop, ngboys, maxpop, station_nb, probaMut, pertu, c, p0, lambdaLM, path_data, file_control);
    run_identification(simul_type,n_param, n_consts, nfiles, ngen, aleaspace, apop, spop, ngboys, maxpop, station_nb, path_data, path_keys, path_results, materialfile, outputfile, simulfile, probaMut, pertu, c, p0, lambdaLM);

}
