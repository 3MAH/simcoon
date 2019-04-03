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
///@brief To read from material.dat and path.dat
///@version 1.0

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

/// Function that fills the matrix Tdsde for mix strain/stress conditions
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

/// Function that fills the matrix Tdsde for mix strain/stress conditions
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

void solver_essentials(int &solver_type, const string &path, const string &filename) {
        
    string pathfile = path + "/" + filename;
    ifstream solver_essentials;
    string buffer;
    
    solver_essentials.open(pathfile, ios::in);
    if(!solver_essentials) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
    
    ///Get the control values for the solver
    solver_essentials >> buffer >> solver_type;
    solver_essentials.close();
}

void solver_control(double &div_tnew_dt_solver, double &mul_tnew_dt_solver, int &miniter_solver, int &maxiter_solver, int &inforce_solver, double &precision_solver, double &lambda_solver, const string &path, const string &filename) {
    
    string pathfile = path + "/" + filename;
    ifstream solver_control;
    string buffer;
    
    solver_control.open(pathfile, ios::in);
    if(!solver_control) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
    
    ///Get the control values for the solver
    solver_control >> buffer >> div_tnew_dt_solver;
    solver_control >> buffer >> mul_tnew_dt_solver;
    solver_control >> buffer >> miniter_solver;
    solver_control >> buffer >> maxiter_solver;
    solver_control >> buffer >> inforce_solver;
    solver_control >> buffer >> precision_solver;
    solver_control >> buffer >> lambda_solver;
    solver_control.close();

}
    
void read_matprops(string &umat_name, unsigned int &nprops, vec &props, unsigned int &nstatev, double &psi_rve, double &theta_rve, double &phi_rve, const string &path_data, const string &materialfile) {

    ///Material properties reading, use "material.dat" to specify parameters values
	string buffer;
	ifstream propsmat;
    string path_materialfile = path_data + "/" + materialfile;
	propsmat.open(path_materialfile, ios::in);
	if(propsmat) {
        
		string buffer;
		propsmat >> buffer >> buffer >> umat_name >> buffer >> nprops >> buffer >> nstatev;
	}
	else {
		cout << "Error: cannot open the file " << materialfile << " in the folder :" << path_data << endl;
	}
	
	char *cmname = new char [umat_name.length()];
	strcpy (cmname, umat_name.c_str());
    
	propsmat.close();
    	
	props = zeros(nprops);
    
	propsmat.open(path_materialfile, ios::in);
	if(propsmat) {
		string buffer;
		propsmat >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> psi_rve >> buffer >> theta_rve >> buffer >> phi_rve >> buffer;
        
		for(unsigned int i=0;i<nprops;i++)
			propsmat >> buffer >> props(i);
	}
	else {
		cout << "Error: cannot open the file " << materialfile << " in the folder :" << path_data << endl;
        return;
	}
    
    psi_rve*=(sim_pi/180.);
    theta_rve*=(sim_pi/180.);
    phi_rve*=(sim_pi/180.);
    
	propsmat.close();
    
}

void read_output(solver_output &so, const int &nblock, const int &nstatev, const string &path_data, const string &outputfile) {
        
    string buffer;
    string file;
    string path_outputfile = path_data + "/" + outputfile;
    
    ifstream cyclic_output;
	cyclic_output.open(path_outputfile, ios::in);
	if(cyclic_output)
	{
        cyclic_output >> buffer;
        cyclic_output >> buffer >> so.o_nb_meca;
        so.o_meca.zeros(so.o_nb_meca);
        for (int i=0; i<so.o_nb_meca; i++) {
            cyclic_output >> so.o_meca(i);
        }
        cyclic_output >> buffer >> so.o_nb_T;
        
        ///Selection of the wanted umat statev, use "cyclic.dat" to specify wanted internal variables
        cyclic_output >> buffer >> buffer;
        if ((buffer == "all") || (buffer == "All") || (buffer == "ALL")){
            so.o_wanted_statev.zeros(1);
            so.o_nw_statev = -1;
            so.o_wanted_statev(0) = -1;
        }
        else if(atoi(buffer.c_str()) != 0){
            so.o_nw_statev = atoi(buffer.c_str());
            so.o_wanted_statev.zeros(so.o_nw_statev);
            so.o_range_statev.zeros(so.o_nw_statev);
            for (int i = 0; i < so.o_nw_statev; i++){
                cyclic_output >> buffer >> buffer;
                if ((buffer == "from") || (buffer == "From") || (buffer == "FROM")){
                    cyclic_output >> so.o_wanted_statev(i) >> buffer >> so.o_range_statev(i);
                }
                else{
                    so.o_wanted_statev(i) = atoi(buffer.c_str());
                    so.o_range_statev(i) = so.o_wanted_statev(i);
                }
                
                if(so.o_range_statev(i) > nstatev -1) {
                    cout << "Error : The range of outputed statev is greater than the actual number of statev!\n";
                    cout << "Check output file and/or material input file\n" << endl;
                    
                    return;
                }
            }
        }
        else {
            so.o_nw_statev = 0;
        }
        
        cyclic_output >> buffer >> buffer >> buffer;
        for(int i = 0 ; i < nblock ; i++){
            cyclic_output >> buffer >> so.o_type(i);
            if(so.o_type(i) == 1)
                cyclic_output >> so.o_nfreq(i);
            else if(so.o_type(i) == 2)
                cyclic_output >> so.o_tfreq(i);
            else
                cyclic_output >> buffer;
        }
        cyclic_output.close();
    }
    else {
//        cout << "The file data/output.dat is not present, so default output is selected\n";
        so.o_nb_meca = 6;
        so.o_meca.zeros(so.o_nb_meca);
        so.o_meca = {0,1,2,3,4,5};
        so.o_nb_T = 1;
        so.o_nw_statev = 0;
        
        for(int i = 0 ; i < nblock ; i++){
            so.o_type(i) = 1;
            so.o_nfreq(i) = 1;
        }
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
                        if((so.o_type(i) == 2)||(sptr_meca->ninc%so.o_nfreq(i))) {
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
    
void read_path(std::vector<block> &blocks, double &T, const string &path_data, const string &pathfile) {
    
	/// Reading the loading path file, Path.txt
    string buffer;
    string pathfile_inc;
    int conver;
    char bufferchar;
    unsigned int nblock;
    Col<int> Equiv = subdiag2vec();
    
    std::string path_inputfile = path_data + "/" + pathfile;
    std::ifstream path;
	path.open(path_inputfile, ios::in);
	if(!path)
	{
		cout << "Error: cannot open the file " << pathfile << " in the folder :" << path_data << "\n";
	}

	///temperature is initialized
	path >> buffer >> T >> buffer >> nblock;
    blocks.resize(nblock);
    
    /// Reading path_file
    for(unsigned int i = 0 ; i < nblock ; i++){
        
        path >> buffer >> blocks[i].number >> buffer >> blocks[i].type >> buffer >> blocks[i].control_type >> buffer >> blocks[i].ncycle >> buffer >> blocks[i].nstep;

        if (blocks[i].number != i+1) {
            cout << "The number of blocks could not be found. Please verify the blocks order in the path file";
        }
        
        blocks[i].generate();
        
        switch(blocks[i].type) {
            case 1: {

                for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                    path >> buffer >> blocks[i].steps[j]->mode;
                    blocks[i].steps[j]->number = j+1;
                    
                    if ((blocks[i].steps[j]->mode == 1)||(blocks[i].steps[j]->mode == 2)) {
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        unsigned int size_meca = sptr_meca->BC_meca.n_elem;
                        
                        path >> buffer >> sptr_meca->Dn_init >> buffer >> sptr_meca->Dn_mini >> buffer >> sptr_meca->Dn_inc >> buffer >> sptr_meca->BC_Time >> buffer;
                        
                        if (sptr_meca->control_type <= 3) {
                            for(unsigned int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 83){
                                    sptr_meca->cBC_meca(Equiv(k)) = 1;
                                    path >> sptr_meca->BC_meca(Equiv(k));
                                }
                                else if (conver == 69){
                                    sptr_meca->cBC_meca(Equiv(k)) = 0;
                                    path >> sptr_meca->BC_meca(Equiv(k));
                                }
                            }
                        }
                        if (sptr_meca->control_type >= 4) {
                            for(unsigned int k = 0 ; k < size_meca ; k++) {
                                sptr_meca->cBC_meca(k) = 0;
                                path >> sptr_meca->BC_meca(k);
                            }
                        }
                        
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            path >> sptr_meca->BC_T;
                        }
                        else
                            cout << "Error, This is a mechanical step, only temperature boundary condition is allowed here\n";

                    }
                    else if (blocks[i].steps[j]->mode == 3) {
                        
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        unsigned int size_meca = sptr_meca->BC_meca.n_elem;
                        
                        path >> buffer >> pathfile_inc >> buffer >> sptr_meca->Dn_init >> buffer >> sptr_meca->Dn_mini >> buffer;
                        sptr_meca->file = path_data + "/" + pathfile_inc;
                        
                        if (sptr_meca->control_type <= 3) {
                            for(unsigned int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 83){
                                    sptr_meca->cBC_meca(Equiv(k)) = 1;
                                }
                                else if (conver == 69) {
                                    sptr_meca->cBC_meca(Equiv(k)) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_meca->cBC_meca(Equiv(k)) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        else if (sptr_meca->control_type == 4) {
                            for(unsigned int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 70){
                                    sptr_meca->cBC_meca(k) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_meca->cBC_meca(k) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        else if (sptr_meca->control_type == 5) {
                            for(unsigned int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 85){
                                    sptr_meca->cBC_meca(k) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_meca->cBC_meca(k) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            sptr_meca->cBC_T = 0;                       //This is the classical temperature imposed in the file
                        }
                        else if (conver == 48) {                        //This is a special case where the temperature is constant
                            sptr_meca->cBC_T = 2;
                        }
                    }
                    else {
                        cout << "Please enter a suitable block mode (1 for linear, 2 for sinusoidal, 3 for user-input)";
                    }
                }
                break;
            }
            case 2: {
                
                for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                    path >> buffer >> blocks[i].steps[j]->mode;
                    blocks[i].steps[j]->number = j+1;
                    
                    if ((blocks[i].steps[j]->mode == 1)||(blocks[i].steps[j]->mode == 2)) {
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        unsigned int size_meca = sptr_thermomeca->BC_meca.n_elem;
                        
                        path >> buffer >> sptr_thermomeca->Dn_init >> buffer >> sptr_thermomeca->Dn_mini >> buffer >> sptr_thermomeca->Dn_inc >> buffer >> sptr_thermomeca->BC_Time >> buffer;
                    
                        if (sptr_thermomeca->control_type <= 3) {
                            for(int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 83){
                                    sptr_thermomeca->cBC_meca(Equiv(k)) = 1;
                                    path >> sptr_thermomeca->BC_meca(Equiv(k));
                                }
                                else if (conver == 69){
                                    sptr_thermomeca->cBC_meca(Equiv(k)) = 0;
                                    path >> sptr_thermomeca->BC_meca(Equiv(k));
                                }
                            }
                        }
                        if (sptr_thermomeca->control_type >= 4) {
                            for(int k = 0 ; k < size_meca ; k++) {
                                sptr_thermomeca->cBC_meca(k) = 0;
                                path >> sptr_thermomeca->BC_meca(k);
                            }
                        }
                        
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 81){
                            sptr_thermomeca->cBC_T = 1;
                            path >> sptr_thermomeca->BC_T;
                        }
                        else if (conver == 84){
                            sptr_thermomeca->cBC_T = 0;
                            path >> sptr_thermomeca->BC_T;
                        }
                        else if (conver == 67) {
                            sptr_thermomeca->cBC_T = 3;
                            path >> sptr_thermomeca->BC_T;
                        }
                        
                    }
                    else if (blocks[i].steps[j]->mode == 3) {
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        unsigned int size_meca = sptr_thermomeca->BC_meca.n_elem;
                        
                        path >> buffer >> pathfile_inc >> buffer >> sptr_thermomeca->Dn_init >> buffer >> sptr_thermomeca->Dn_mini >> buffer;
                        sptr_thermomeca->file = path_data + "/" + pathfile_inc;
                        
                        if (sptr_thermomeca->control_type <= 3) {
                            for(int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 83){
                                    sptr_thermomeca->cBC_meca(Equiv(k)) = 1;
                                }
                                else if (conver == 69) {
                                    sptr_thermomeca->cBC_meca(Equiv(k)) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_thermomeca->cBC_meca(Equiv(k)) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        else if (sptr_thermomeca->control_type == 4) {
                            for(int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 70){
                                    sptr_thermomeca->cBC_meca(k) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_thermomeca->cBC_meca(k) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        else if (sptr_thermomeca->control_type == 5) {
                            for(int k = 0 ; k < size_meca ; k++) {
                                path >> bufferchar;
                                conver = bufferchar;
                                if (conver == 85){
                                    sptr_thermomeca->cBC_meca(k) = 0;
                                }
                                else if (conver == 48) {
                                    sptr_thermomeca->cBC_meca(k) = 2;      // this is a special stress-controlled one (it is different since it does not read data from the path_inc tabular file)
                                }
                            }
                        }
                        
                        path >> buffer >> bufferchar;
                        conver = bufferchar;
                        if (conver == 84){
                            sptr_thermomeca->cBC_T = 0;                       //This is the classical temperature imposed in the file
                        }
                        else if (conver == 81){
                            sptr_thermomeca->cBC_T = 1;                       //This is the classical heat flux quantitt imposed in the file
                        }
                        else if (conver == 48) {                        //This is a special case where the temperature is constant
                            sptr_thermomeca->cBC_T = 2;
                        }
                        else if (conver == 67) {                        //This is a special case where the convexion is assumed
                            sptr_thermomeca->cBC_T = 3;
                            path >> sptr_thermomeca->BC_T;                              //Get the tau
                        }
                        
                        
                    }
                    else {
                        cout << "Please enter a suitable block mode (1 for linear, 2 for sinusoidal, 3 for user-input)";
                    }
                    
                }
                break;
            }
            default: {
                cout << "Please enter a valid block type (1 for mechanical, 2 for thermomechanical)\n";
                break;
            }
            
        }
    }
    path.close();
    
}

} //namespace simcoon
