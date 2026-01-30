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
///@brief Solver utility functions for mixed boundary conditions
///@version 2.0

#include <fstream>
#include <sstream>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Simulation/Solver/output.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void read_matprops(string &umat_name, unsigned int &nprops, vec &props, unsigned int &nstatev,
                   double &psi_rve, double &theta_rve, double &phi_rve,
                   const string &path_data, const string &materialfile) {

    string filename = path_data + "/" + materialfile;
    ifstream propsdata(filename);

    if (!propsdata) {
        throw runtime_error("Could not open material file: " + filename);
    }

    string buffer;

    // Skip header line
    getline(propsdata, buffer);

    // Read material properties
    propsdata >> umat_name >> nprops >> nstatev >> psi_rve >> theta_rve >> phi_rve;

    // Convert angles from degrees to radians
    psi_rve *= (M_PI / 180.0);
    theta_rve *= (M_PI / 180.0);
    phi_rve *= (M_PI / 180.0);

    props = zeros(nprops);
    for (unsigned int i = 0; i < nprops; i++) {
        propsdata >> props(i);
    }

    propsdata.close();
}

Col<int> subdiag2vec() {

    Col<int> temp;
    temp = zeros<Col<int> >(6);
    temp(0) = 0;
	temp(1) = 3;
	temp(2) = 1;
	temp(3) = 4;
	temp(4) = 5;
	temp(5) = 2;

    return temp;
}

/// Function that fills the matrix K for mixed strain/stress boundary conditions
void Lt_2_K(const mat &Lt, mat &K, const Col<int> &cBC_meca, const double &lambda)
{
	K = zeros(6,6);

    for (int i=0; i<6; i++) {
        if (cBC_meca(i)) {
            K.row(i) = Lt.row(i);
        }
        else
            K(i,i) = lambda;
    }
}

/// Function that fills the matrix K for mixed strain/stress/thermal boundary conditions
void Lth_2_K(const mat &dSdE, mat &dSdT, mat &dQdE, mat &dQdT, mat &K, const Col<int> &cBC_meca, const int &cBC_T, const double &lambda)
{
	K = zeros(7,7);

    K.submat(0, 0, 5, 5) = dSdE;
    K.submat(0, 6, 5, 6) = dSdT;
    K.submat(6, 0, 6, 5) = dQdE;
    K.submat(6, 6, 6, 6) = dQdT;

    for (int i=0; i<6; i++) {
        if (cBC_meca(i) == 0) {
            K.row(i) = 0.*K.row(i);
            K(i,i) = lambda;
        }
    }
    if (cBC_T == 0) {
            K.row(6) = 0.*K.row(6);
            K(6,6) = lambda;
    }
}

void check_path_output(const std::vector<block> &blocks, const solver_output &so) {

    /// Reading blocks
    for(unsigned int i = 0 ; i < blocks.size() ; i++) {
        
        switch(blocks[i].type) {
            case 1: {
                
                for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                    shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                    
                    if (sptr_meca->mode == 3) {
                        if((so.o_type(i) == 2)||((so.o_type(i) == 1)&&(so.o_nfreq(i) != 1))) {
                            cout << "The output nfreq is not compatible with the number of increments of the step)" << endl;
                            break;
                        }
                    }
                    else {
                        
                        if(so.o_type(i) == 1) {
                            if(sptr_meca->ninc%so.o_nfreq(i) > 0) {
                                cout << "The output nfreq is not compatible with the number of increments of the step)" << endl;
                                break;
                            }
                        }
                        else if(so.o_type(i) == 2) {
                            if((fmod(1, so.o_tfreq(i)) > sim_limit)||(fmod(so.o_tfreq(i), sptr_meca->Dn_inc) > 0.)) {
                                cout << "The output tfreq is not compatible with the time of increments of the step)" << endl;
                                break;
                            }
                        }
                        
                    }
                    
                }
                break;
            }
            case 2: {
                
                for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                    shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                    
                    if (sptr_thermomeca->mode == 3) {
                        if((so.o_type(i) == 2)||(sptr_thermomeca->ninc%so.o_nfreq(i))) {
                            cout << "The output nfreq is not compatible with the number of increments of the step)" << endl;
                            break;
                        }
                    }
                    else {
                        if(so.o_type(i) == 1) {
                            if(sptr_thermomeca->ninc%so.o_nfreq(i) > 0) {
                                cout << "The output nfreq is not compatible with the number of increments of the step)" << endl;
                                break;
                            }
                        }
                        else if(so.o_type(i) == 2) {
                            if((fmod(1, so.o_tfreq(i)) > sim_limit)||(fmod(so.o_tfreq(i), sptr_thermomeca->Dn_inc) > 0.)) {
                                cout << "The output tfreq is not compatible with the time of increments of the step)" << endl;
                                break;
                            }
                        }
                    }
                }
                break;
            }
            default: {
                cout << "The block type is incorrect. Please enter a valid block type (1) : Mechanical (2) Thermomechanical" << endl;
                break;
            }
        }
        
    }
    
}

} //namespace simcoon
