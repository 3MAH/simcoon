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

///@file read
///@brief read and construct for complex objects construction from input files
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/opti_data.hpp>
#include <simcoon/Simulation/Identification/read.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void read_parameters(const int &n_param, vector<parameters> &params) {
    
    ifstream paraminit;
    string buffer;
    
    ///@brief Properties of the parameters reading, use "parameters.dat" to specify the parameters of a model
    paraminit.open("data/parameters.inp", ios::in);
    if(paraminit) {
        paraminit >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
        for(int i=0; i<n_param; i++) {
            paraminit >> buffer >> buffer >> buffer >> buffer >> params[i].ninput_files;
            params[i].input_files.resize(params[i].ninput_files);
            for(int j=0; j<params[i].ninput_files; j++) {
                paraminit >> buffer;
            }
        }
    }
    else {
        cout << "Error: cannot open parameters.inp file \n";
        exit(0);
    }
    paraminit.close();
    
    paraminit.open("data/parameters.inp", ios::in);
    if(paraminit) {
        paraminit >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
        for(int i=0; i<n_param; i++) {
            paraminit >> params[i].number >> params[i].min_value >> params[i].max_value >> params[i].key >> buffer;
            for(int j=0; j<params[i].ninput_files; j++) {
                paraminit >> params[i].input_files[j];
            }
        }
    }
    paraminit.close();
}
    
void read_constants(const int &n_consts, vector<constants> &consts, const int &nfiles) {
    
    ifstream paraminit;
    string buffer;
    
    ///@brief Properties of the parameters reading, use "parameters.dat" to specify the parameters of a model
    paraminit.open("data/constants.inp", ios::in);
    if(paraminit) {
        paraminit >> buffer >> buffer >> buffer >> buffer >> buffer;
        for(int i=0; i<n_consts; i++) {
            paraminit >> buffer >> buffer;
            for(int j=0; j<nfiles; j++) {
                paraminit >> buffer;
            }
            paraminit >> consts[i].ninput_files;
            consts[i].resize(consts[i].ninput_files, nfiles);
            for(int j=0; j<consts[i].ninput_files; j++) {
                paraminit >> buffer;
            }
        }
    }
    else {
        cout << "Error: cannot open parameters.inp file \n";
        exit(0);
    }
    paraminit.close();
    
    paraminit.open("data/constants.inp", ios::in);
    if(paraminit) {
        paraminit >> buffer >> buffer >> buffer >> buffer >> buffer;
        for(int i=0; i<n_consts; i++) {
            paraminit >> consts[i].number >> consts[i].key;
            for(int j=0; j<nfiles; j++) {
                paraminit >> consts[i].input_values(j);
            }
            paraminit>> buffer;
            for(int j=0; j<consts[i].ninput_files; j++) {
                paraminit >> consts[i].input_files[j];
            }
        }
    }
    paraminit.close();
}
    
void read_data_exp(const int &nfiles, vector<opti_data> &datas) {
    
    ifstream paraminit;
    string buffer;
    paraminit.open("data/files_exp.inp", ios::in);
    if(!paraminit) {
        cout << "Error: cannot open files_exp.inp file\n";
        exit(0);
    }
    datas.resize(nfiles);
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> datas[i].name;
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> datas[i].ncolumns;
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> datas[i].ninfo;
        datas[i].constructc_data();
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        for(int j=0; j<datas[i].ninfo;j++) {
            paraminit >> datas[i].c_data(j);
            assert(datas[i].c_data(j)>=0);
            assert(datas[i].c_data(j)<= datas[i].ncolumns);
        }
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> datas[i].skiplines;
    }
    paraminit.close();
}
    
void read_data_weights(const int &nfiles, Col<int> &weight_types, vec &weights_file, vector<vec> &weights_cols, vector<opti_data> &weights, const vector<opti_data> &data_exp) {
    
    ifstream paraminit;
    string buffer;
    paraminit.open("data/files_weights.inp", ios::in);
    if(!paraminit) {
        cout << "Error: cannot open files_weights.inp file\n";
        exit(0);
    }

    weights.resize(nfiles);
    for (int i=0; i<nfiles; i++) {
        weights[i].name = data_exp[i].name;
        weights[i].number = data_exp[i].number;
        weights[i].ndata = data_exp[i].ndata;
        weights[i].ninfo = data_exp[i].ninfo;
        weights[i].ncolumns = data_exp[i].ncolumns;
        weights[i].skiplines = data_exp[i].skiplines;
        weights[i].constructc_data();
    }
    
    //Load info for the weight type 1 : Weight for each data file
    paraminit >> buffer >> buffer >> weight_types(0);
    if (weight_types(0) == 0) {
        paraminit >> buffer;
    }
    else if(weight_types(0) == 1) {
        paraminit >> buffer;
        weights_file.resize(nfiles);
        for(int i = 0; i<nfiles; i++) {
            paraminit >> weights_file(i);
        }
    }
    else {
        cout << "Please enter 0 or 1 for the weight type 1 : Weight for each data point\n";
        exit(0);
    }
    
    //Load info for the weight type 2 : Weight for each data columns
    paraminit >> buffer >> buffer >> weight_types(1);
    
    if (weight_types(1) == 0) {
        paraminit >> buffer;
    }
    else if (weight_types(1) == 1) {
        paraminit >> buffer;
    }
    else if((weight_types(1) == 2)||(weight_types(1) == 3)) {
        paraminit >> buffer;
        for(int i = 0; i<nfiles; i++) {
            weights_cols[i].resize(weights[i].ninfo);
            for (int j=0; j<weights[i].ninfo; j++) {
                paraminit >> weights_cols[i](j);
            }
        }
    }
    else {
        cout << "Please enter 0 or 1 or 2 or 3 for the weight type 2 : Weight for each data point\n";
        exit(0);
    }
    
    //Load info for the weight type 3 : Weight for each data point
    paraminit >> buffer >> buffer >> weight_types(2);
    if (weight_types(2) == 0) {
        paraminit >> buffer;
    }
    else if(weight_types(2) == 1) {
        paraminit >> buffer;        
        for(int i=0; i<nfiles;i++) {
            for(int j=0; j<weights[i].ninfo;j++) {
                paraminit >> weights[i].c_data(j);
                assert(weights[i].c_data(j)>0);
                assert(weights[i].c_data(j)<= weights[i].ncolumns);
            }
        }
    }
    else {
        cout << "Please enter 0 or 1 for the weight type 3 : Weight for each data point\n";
        exit(0);
    }
    paraminit.close();
}

void read_data_num(const int &nfiles, const vector<opti_data> &data_exp, vector<opti_data> &data_num) {
    
    ifstream paraminit;
    string buffer;
    paraminit.open("data/files_num.inp", ios::in);
    if(!paraminit) {
        cout << "Error: cannot open files_num.inp file\n";
        exit(0);
    }
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> data_num[i].ncolumns;
        data_num[i].ninfo = data_exp[i].ninfo;
        data_num[i].constructc_data();
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        for(int j=0; j<data_num[i].ninfo; j++) {
            paraminit >> data_num[i].c_data(j);
            assert(data_num[i].c_data(j)>0);
            assert(data_num[i].c_data(j)<=data_num[i].ncolumns);
        }
    }
    
    paraminit >> buffer;
    for(int i=0; i<nfiles;i++) {
        paraminit >> data_num[i].skiplines;
    }
    paraminit.close();
}

void ident_essentials(int &n_param, int &n_consts, int &n_files, const string &path, const string &filename) {

    string pathfile = path + "/" + filename;
    ifstream param_essentials;
    string buffer;

    param_essentials.open(pathfile, ios::in);
    if(!param_essentials) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
        
    ///Get the control values for the genetic algorithm
    param_essentials >> buffer >> n_param;
    param_essentials >> buffer >> n_consts;
    param_essentials >> buffer >> n_files;
    
    param_essentials.close();
    
}
    
void ident_control(int &ngen, int &aleaspace, int &apop, int &spop, int &ngboys, int &maxpop, int &station_nb, double &station_lim, double &probaMut, double &pertu, double &c, double &p0, double &lambdaLM, const string &path, const string &filename) {
    
    string pathfile = path + "/" + filename;
    ifstream param_control;
    string buffer;
        
    param_control.open(pathfile, ios::in);
    if(!param_control) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
        
    ///Get the control values for the genetic algorithm
    param_control >> buffer >> ngen;
    param_control >> buffer >> aleaspace;
    ///Get the state of the initial population : 0 = equidistant individuals, 1 = random individuals, 2 = previously computed population, 3 = equidistant individuals with boundary ones
    if((aleaspace==0)||(aleaspace==1))
        param_control >> buffer >> spop;
    else if((aleaspace==2)||(aleaspace==3))
        param_control >> buffer >> apop;
    else {
        cout << "Please select if the initial space is filled with random or equidistant values\n";
        exit(0);
    }
    
    param_control >> buffer >> ngboys;
    param_control >> buffer >> maxpop;
    
    param_control >> buffer >> station_nb;
    param_control >> buffer >> station_lim;
    param_control >> buffer >> probaMut;
    param_control >> buffer >> pertu;
    
    param_control >> buffer >> c >> p0;
    param_control >> buffer >> lambdaLM;
    
    param_control.close();
}
    
void read_gen(int &apop, mat &samples, const int &n_param) {
    
    ifstream paraminit;
    string buffer;
    paraminit.open("data/gen0.inp", ios::in);
    if(!paraminit) {
        cout << "Error: cannot open data/gen0.inp file\n";
        exit(0);
    }
    
    apop=0;
    //Read the number of lines (-1 since their is a header), to get the pop number
    while (!paraminit.eof())
    {
        getline (paraminit,buffer);
        if (buffer != "") {
            apop++;
        }
    }
    paraminit.close();
    apop--;
    samples.resize(apop, n_param);
    
    paraminit.open("data/gen0.inp", ios::in);
    paraminit >> buffer >> buffer >> buffer;
    for(int j=0; j<n_param; j++) {
        paraminit >> buffer;
    }
    
    for(int i=0;i<apop;i++) {
        paraminit >> buffer >> buffer >> buffer;
        for(int j=0; j<n_param; j++) {
            paraminit >> samples(i,j);
        }
    }
    paraminit.close();
}
    
} //namespace simcoon
