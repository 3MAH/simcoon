/* This file is part of simcoon private.
 
 Only part of simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file func_N.cpp
///@brief Functions that computes accumulative laws
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/func_N.hpp>

//Add the header of the cumulative function you want

using namespace std;
using namespace arma;

namespace simcoon{
    
//This function returns a file with the values computed according to an implemented function (to define)
void func_N(const vec &params, const vec &variables, const string& N_file, const string& outputfile, const string& path_data, const string& path_results) {
    
    std::string buffer;
    std::string path_inputfile = path_data + "/" + N_file;
    std::string path_outputfile = path_results + "/" + outputfile;
    std::ifstream cum_N;
    unsigned int nN = 0;
    
    cum_N.open(path_inputfile, ios::in);
    if(cum_N) {
        while (!cum_N.eof())
        {
            getline (cum_N,buffer);
            if (buffer != "") {
                nN++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << N_file << " that details the cumulative function characteristics in the folder: " << path_data << endl;
        return;
    }
    cum_N.close();
    
    vec N = zeros(nN);
    cum_N.open(path_inputfile, ios::in);
    for(unsigned int i=0; i<nN; i++) {
        cum_N >> N(i) >> buffer;
    }
    cum_N.close();
    
    //In here you are asked to introduced the function you want to solve
    UNUSED(params);
    UNUSED(variables);
    vec y;
//    vec y = p_cumulative(N, variables(0), variables(1), params); Insert here the fonction you want
    
    //write in the file
    std::ofstream cum_N_out;
    cum_N_out.open(path_outputfile, ios::out);
    
    for(unsigned int i=0; i<N.n_elem; i++) {
        cum_N_out << N(i) << "\t" << y(i) << "\n";
    }
    cum_N_out.close();
}

} //namespace simcoon
