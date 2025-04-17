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

///@file Tidentification.cpp
///@brief Test for Identification algorithm
///@version 1.0

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/identification.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Solver/read.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(Tidentification, identification_layers)
{
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
    double station_lim;
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
    ident_control(ngen, aleaspace, apop, spop, ngboys, maxpop, station_nb, station_lim, probaMut, pertu, c, p0, lambdaLM, path_data, file_control);
    run_identification(simul_type,n_param, n_consts, nfiles, ngen, aleaspace, apop, spop, ngboys, maxpop, station_nb, station_lim, path_data, path_keys, path_results, materialfile, outputfile, simulfile, probaMut, pertu, c, p0, lambdaLM);
    
    string umat_name;
    unsigned int nprops = 0;
    unsigned int nstatev = 0;
    vec props;
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    
    phase_characteristics rve_layer_0;
    phase_characteristics rve_layer_1;
    rve_layer_0.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_layer_0.construct(1,1); //The rve is supposed to be mechanical only here
    rve_layer_1.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
    rve_layer_1.construct(1,1); //The rve is supposed to be mechanical only here

    string path_comparison = "comparison";
    
    read_layer(rve_layer_0, path_comparison, "Nlayers0.dat");
    read_layer(rve_layer_0, path_results, "Nlayers0.dat");

    for(unsigned int i=0; i<nprops; i++) {
        if (fabs(rve_layer_0.sptr_matprops->props(i)) > simcoon::iota) {
            EXPECT_LT( pow(pow(rve_layer_0.sptr_matprops->props(i),2.) - pow(rve_layer_1.sptr_matprops->props(i),2.),0.5)/fabs(rve_layer_0.sptr_matprops->props(i)),1.E-6);
        }
        else {
            EXPECT_LT( pow(pow(rve_layer_0.sptr_matprops->props(i),2.) - pow(rve_layer_1.sptr_matprops->props(i),2.),0.5),1.E-6);
        }
    }
}
