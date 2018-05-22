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

///@file script.cpp
///@brief Scripts that allows to run identification algorithms based on Smart+ Control functions
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <algorithm>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/generation.hpp>
#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/optimize.hpp>
#include <simcoon/Simulation/Identification/script.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/solver.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Phase/write.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_mechanics/Material/PDF.hpp>
#include <simcoon/Continuum_mechanics/Material/read.hpp>
#include <simcoon/Continuum_mechanics/Material/ODF2Nphases.hpp>
#include <simcoon/Continuum_mechanics/Material/PDF2Nphases.hpp>
#include <simcoon/Continuum_mechanics/Functions/func_N.hpp>
#include <simcoon/Continuum_mechanics/Functions/read.hpp>

using namespace std;
using namespace arma;
using namespace arma;

namespace simcoon{
    
//This function will copy the parameters files
void copy_parameters(const vector<parameters> &params, const string &src_path, const string &dst_path) {

    string src_files;
    string dst_files;
    
    for (auto pa : params) {
        for(auto ifiles : pa.input_files) {
            src_files = src_path + "/" + ifiles;
            dst_files = dst_path + "/" + ifiles;
            boost::filesystem::copy_file(src_files,dst_files,boost::filesystem::copy_option::overwrite_if_exists);
        }
    }
}

//This function will copy the parameters files
void copy_constants(const vector<constants> &consts, const string &src_path, const string &dst_path) {
    
    string src_files;
    string dst_files;
    
    for (auto co : consts) {
        for(auto ifiles : co.input_files) {
            src_files = src_path + "/" + ifiles;
            dst_files = dst_path + "/" + ifiles;
            boost::filesystem::copy_file(src_files,dst_files,boost::filesystem::copy_option::overwrite_if_exists);
        }
    }
}
    
//This function will replace the keys by the parameters
void apply_parameters(const vector<parameters> &params, const string &dst_path) {
  
    string mod_files;
    string buffer;
    
    ifstream in_files;
    ofstream ou_files;
    
    for (auto pa : params) {
        for(auto ifiles : pa.input_files) {
            mod_files = dst_path + "/" + ifiles;

            in_files.open(mod_files, ios::in);
            
            std::vector<string> str;
            while (!in_files.eof())
            {
                getline(in_files,buffer);
                str.push_back(buffer);
            }
            in_files.close();
            
            ou_files.open(mod_files);
            for (auto s : str) {
                boost::replace_all(s, pa.key, to_string(pa.value));
                ou_files << s << "\n";
            }
            ou_files.close();
        }
    }

}

//This function will replace the keys by the parameters
void apply_constants(const vector<constants> &consts, const string &dst_path) {
    
    string mod_files;
    string buffer;
    
    ifstream in_files;
    ofstream ou_files;
    
    for (auto co : consts) {
        for(auto ifiles : co.input_files) {
            mod_files = dst_path + "/" + ifiles;
            
            in_files.open(mod_files, ios::in);
            
            std::vector<string> str;
            while (!in_files.eof())
            {
                getline(in_files,buffer);
                str.push_back(buffer);
            }
            in_files.close();
            
            ou_files.open(mod_files);
            for (auto s : str) {
                boost::replace_all(s, co.key, to_string(co.value));
                ou_files << s << "\n";
            }
            ou_files.close();
        }
    }
    
}
    
void launch_solver(const individual &ind, const int &nfiles, vector<parameters> &params, vector<constants> &consts, const string &path_results, const string &name, const string &path_data, const string &path_keys, const string &materialfile)
{
	string outputfile;
    string simulfile;
	string pathfile;
    
    string name_ext = name.substr(name.length()-4,name.length());
    string name_root = name.substr(0,name.length()-4); //to remove the extension
    
	//#pragma omp parallel for private(sstm, path)
    for (int i = 0; i<nfiles; i++) {
        ///Creating the right path & output filenames
        
        outputfile = name_root + "_" + to_string(ind.id) + "_" + to_string(i+1) + name_ext;
        pathfile = "path_id_" + to_string(i+1) + ".txt";
        
        string umat_name;
        unsigned int nprops = 0;
        unsigned int nstatev = 0;
        vec props;
        
        double psi_rve = 0.;
        double theta_rve = 0.;
        double phi_rve = 0.;
                
        //Replace the constants
        for (unsigned int k=0; k<consts.size(); k++) {
            consts[k].value = consts[k].input_values(i);
        }
        //Replace the parameters
        for (unsigned int k=0; k<params.size(); k++) {
            params[k].value = ind.p(k);
        }
        
        copy_constants(consts, path_keys, path_data);
        copy_parameters(params, path_keys, path_data);
        
        apply_constants(consts, path_data);
        apply_parameters(params, path_data);
        
        int solver_type = 0;
        
        double div_tnew_dt_solver = 0.;
        double mul_tnew_dt_solver = 0.;
        int miniter_solver = 0;
        int maxiter_solver = 0;
        int inforce_solver = 0;
        double precision_solver = 0.;
        double lambda_solver = 0.;
        
        solver_essentials(solver_type, path_data);
        solver_control(div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data);
        
        //Then read the material properties
        read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
        ///Launching the solver with relevant parameters
        solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, path_results, pathfile, outputfile);
        
        //Get the simulation files according to the proper name
        outputfile = path_results + "/" + name_root + "_" + to_string(ind.id) + "_" + to_string(i+1) + "_global-0" + name_ext;
        simulfile = path_results + "/" + name_root + "_" + to_string(ind.id)  + "_" + to_string(i+1) + name_ext;
        
        boost::filesystem::copy_file(outputfile,simulfile,boost::filesystem::copy_option::overwrite_if_exists);
    }
}
    
void launch_odf(const individual &ind, vector<parameters> &params, const string &path_results, const string &name, const string &path_data, const string &path_keys, const string &materialfile)
{

    string inputfile;
    string outputfile;
    string simulfile;
    
    string name_ext = name.substr(name.length()-4,name.length());
    string name_root = name.substr(0,name.length()-4); //to remove the extension
    
    //Replace the parameters
    for (unsigned int k=0; k<params.size(); k++) {
        params[k].value = ind.p(k);
    }
    
    copy_parameters(params, path_keys, path_data);
    apply_parameters(params, path_data);
    
    string umat_name;
    unsigned int nprops;
    unsigned int nstatev = 0;
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props;
    
    //Then read the material properties
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    phase_characteristics rve_init;
    rve_init.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, nprops, props);
    // The vector of props should be = {nphases_out,nscale,geom_type,npeak};
    //int nphases_in = int(props(0));
    int nphases_out = int(props(1));
    int nscale_in = int(props(2));
    int nscale_out = int(props(3));
    int geom_type = int(props(4));
    int npeak = int(props(5));
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            rve_init.construct(0,1); //The rve is supposed to be mechanical only here
            inputfile = "Nphases" + to_string(nscale_in) + ".dat";
            read_phase(rve_init, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            rve_init.construct(1,1); //The rve is supposed to be mechanical only here
            inputfile = "Nlayers" + to_string(nscale_in) + ".dat";
            read_layer(rve_init, path_data, inputfile);
            break;
        }
        case 2: {
            rve_init.construct(2,1); //The rve is supposed to be mechanical only here
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(nscale_in) + ".dat";
            read_ellipsoid(rve_init, path_data, inputfile);
            break;
        }
        case 3: {
            rve_init.construct(3,1); //The rve is supposed to be mechanical only here
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(nscale_in) + ".dat";
            read_cylinder(rve_init, path_data, inputfile);
            break;
        }
    }
    
    double angle_min = 0.;
    double angle_max = 180.;
    string peakfile = "Npeaks" + to_string(npeak) + ".dat";
    
    ODF odf_rve(0, false, angle_min, angle_max);
    read_peak(odf_rve, path_data, peakfile);
    
    phase_characteristics rve = discretize_ODF(rve_init, odf_rve, 1, nphases_out,0);
    
    if(rve.shape_type == 0) {
        outputfile = "Nphases" + to_string(nscale_out) + ".dat";
        write_phase(rve, path_data, outputfile);
    }
    if(rve.shape_type == 1) {
        outputfile = "Nlayers" + to_string(nscale_out) + ".dat";
        write_layer(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 2) {
        outputfile = "Nellipsoids" + to_string(nscale_out) + ".dat";
        write_ellipsoid(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 3) {
        outputfile = "Ncylinders" + to_string(nscale_out) + ".dat";
        write_cylinder(rve, path_data, outputfile);
    }
    
    //Get the simulation files according to the proper name
    simulfile = path_results + "/" + name_root + "_" + to_string(ind.id)  +"_" + to_string(1) + name_ext;
    outputfile = path_data + "/" + outputfile;
    boost::filesystem::copy_file(outputfile,simulfile,boost::filesystem::copy_option::overwrite_if_exists);
}


void launch_pdf(const individual &ind, vector<parameters> &params, const string &path_results, const string &name, const string &path_data, const string &path_keys, const string &materialfile)
{

    string inputfile;
    string outputfile;
    string simulfile;
    
    string name_ext = name.substr(name.length()-4,name.length());
    string name_root = name.substr(0,name.length()-4); //to remove the extension
    
    //Replace the parameters
    for (unsigned int k=0; k<params.size(); k++) {
        params[k].value = ind.p(k);
    }
    
    copy_parameters(params, path_keys, path_data);
    apply_parameters(params, path_data);
    
    string umat_name;
    unsigned int nprops;
    unsigned int nstatev = 0;    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props;
    
    //Then read the material properties
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    phase_characteristics rve_init;
    rve_init.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, nprops, props);
    
    
    // The vector of props should be = {nphases_out,nscale,geom_type,npeak};
    //int nphases_in = int(props(0));
    int nphases_out = int(props(1));		//called "nphases_rve" (props(6)) in main (software/PDF.cpp)
    int nscale_in = int(props(2));		//called "num_file_in" (props(1)) in main (software/PDF.cpp)
    int nscale_out = int(props(3));		//called "num_file_out" (props(5)) in main (software/PDF.cpp)
    int geom_type = int(props(4));		//related to umat_name in main (software/PDF.cpp)
    int npeak = int(props(5));		//called "num_file_peaks" (props(11)) in main (software/PDF.cpp)
    double parameter_min = int(props(6));		//as (props(8)) in main (software/PDF.cpp)
    double parameter_max = int(props(7));		//as (props(9)) in main (software/PDF.cpp)
    double num_Parameter = int(props(8));		//as (props(10)) in main (software/PDF.cpp)
    
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            rve_init.construct(0,1); //The rve is supposed to be mechanical only here
            inputfile = "Nphases" + to_string(nscale_in) + ".dat";
            read_phase(rve_init, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            rve_init.construct(1,1); //The rve is supposed to be mechanical only here
            inputfile = "Nlayers" + to_string(nscale_in) + ".dat";
            read_layer(rve_init, path_data, inputfile);
            break;
        }
        case 2: {
            rve_init.construct(2,1); //The rve is supposed to be mechanical only here
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(nscale_in) + ".dat";
            read_ellipsoid(rve_init, path_data, inputfile);
            break;
        }
        case 3: {
            rve_init.construct(3,1); //The rve is supposed to be mechanical only here
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(nscale_in) + ".dat";
            read_cylinder(rve_init, path_data, inputfile);
            break;
        }
    }
    

    string peakfile = "Npeaks" + to_string(npeak) + ".dat";
    
    PDF pdf_rve(num_Parameter, parameter_min, parameter_max);    ////HERE Parameter index is set to 0 as default
    read_peak(pdf_rve, path_data, peakfile);
    
    phase_characteristics rve = discretize_PDF(rve_init, pdf_rve, 1, nphases_out);
    
    if(rve.shape_type == 0) {
        outputfile = "Nphases" + to_string(nscale_out) + ".dat";
        write_phase(rve, path_data, outputfile);
    }
    if(rve.shape_type == 1) {
        outputfile = "Nlayers" + to_string(nscale_out) + ".dat";
        write_layer(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 2) {
        outputfile = "Nellipsoids" + to_string(nscale_out) + ".dat";
        write_ellipsoid(rve, path_data, outputfile);
    }
    else if(rve.shape_type == 3) {
        outputfile = "Ncylinders" + to_string(nscale_out) + ".dat";
        write_cylinder(rve, path_data, outputfile);
    }
    
    //Get the simulation files according to the proper name
    simulfile = path_results + "/" + name_root + "_" + to_string(ind.id)  +"_" + to_string(1) + name_ext;
    outputfile = path_data + "/" + outputfile;
    boost::filesystem::copy_file(outputfile,simulfile,boost::filesystem::copy_option::overwrite_if_exists);
}
    
void launch_func_N(const individual &ind, const int &nfiles, vector<parameters> &params, vector<constants> &consts, const string &path_results, const string &name, const string &path_data, const string &path_keys, const string &materialfile)
{

    string outputfile;
    string simulfile;
    string pathfile;
    
    string name_ext = name.substr(name.length()-4,name.length());
    string name_root = name.substr(0,name.length()-4); //to remove the extension
    
    for (int i = 0; i<nfiles; i++) {
        ///Creating the right path & output filenames
        
        outputfile = name_root + "_" + to_string(ind.id) + "_" + to_string(i+1) + name_ext;
        pathfile = "path_id_" + to_string(i+1) + ".txt";
        
        //Replace the constants
        for (unsigned int k=0; k<consts.size(); k++) {
            consts[k].value = consts[k].input_values(i);
        }
        //Replace the parameters
        for (unsigned int k=0; k<params.size(); k++) {
            params[k].value = ind.p(k);
        }
        
        copy_constants(consts, path_keys, path_data);
        copy_parameters(params, path_keys, path_data);
        
        apply_constants(consts, path_data);
        apply_parameters(params, path_data);
        
        vec props;
        vec variables;
        string N_file;
        read_func_N(props, variables, N_file, path_data, materialfile);
        func_N(props, variables, pathfile, outputfile, path_data, path_results);
        
        //Get the simulation files according to the proper name
        //simulfile = path_results + "/" + name_root + "_" + to_string(ind.id)  + "_" + to_string(i+1) + name_ext;
        //boost::filesystem::copy_file(outputfile,simulfile,boost::filesystem::copy_option::overwrite_if_exists);
    }
}
    
void run_simulation(const string &simul_type, const individual &ind, const int &nfiles, vector<parameters> &params, vector<constants> &consts, vector<opti_data> &data_num, const string &folder, const string &name, const string &path_data, const string &path_keys, const string &inputdatafile) {
    
    //In the simulation run, make sure that we remove all the temporary files
    boost::filesystem::path path_to_remove(folder);
    for (boost::filesystem::directory_iterator end_dir_it, it(path_to_remove); it!=end_dir_it; ++it) {
        uintmax_t temp = boost::filesystem::remove_all(it->path());
    }
    
    std::map<std::string, int> list_simul;
    list_simul = {{"SCRIPT",0},{"SOLVE",1},{"ODF",2},{"PDF",3},{"FUNCN",4}};
    
    switch (list_simul[simul_type]) {
            
        case 0: {
            //to finish
            break;
        }
        case 1: {
            launch_solver(ind, nfiles, params, consts, folder, name, path_data, path_keys, inputdatafile);
            break;
        }
        case 2: {
            launch_odf(ind, params, folder, name, path_data, path_keys, inputdatafile);
            break;
        }
        case 3: {
            launch_pdf(ind, params, folder, name, path_data, path_keys, inputdatafile);
            break;
        }
        case 4: {
            launch_func_N(ind, nfiles, params, consts, folder, name, path_data, path_keys, inputdatafile);
            break;
        }
        default: {
            cout << "\n\nError in run_simulation : The specified solver (" << simul_type << ") does not exist.\n";
            return;
        }
    }
    
    for (int i = 0; i<nfiles; i++) {
        
        string simulfile;
        
        string name_ext = name.substr(name.length()-4,name.length());
        string name_root = name.substr(0,name.length()-4); //to remove the extension
        simulfile = name_root + + "_" + to_string(ind.id)  +"_" + to_string(i+1) + name_ext;
        
        data_num[i].name = simulfile;
        data_num[i].import(folder);
    }
    
}
    
double calc_cost(const vec &vexp, vec &vnum, const vec &W, const vector<opti_data> &data_num, const vector<opti_data> &data_exp, const int &nfiles, const int &sizev) {

    vnum = calcV(data_num, data_exp, nfiles, sizev);    
    return calcC(vexp, vnum, W);
}
     
mat calc_sensi(const individual &gboy, generation &n_gboy, const string &simul_type, const int &nfiles, const int &n_param, vector<parameters> &params, vector<constants> &consts, vec &vnum0, vector<opti_data> &data_num, vector<opti_data> &data_exp, const string &folder, const string &name, const string &path_data, const string &path_keys, const int &sizev, const vec &Dp_n, const string &materialfile) {
    
    //delta
    vec delta = 0.01*ones(n_param);
    
    mat S = zeros(sizev,n_param);
    //genrun part of the gradient
    
    run_simulation(simul_type, gboy, nfiles, params, consts, data_num, folder, name, path_data, path_keys, materialfile);
    vnum0 = calcV(data_num, data_exp, nfiles, sizev);
    
    for(int j=0; j<n_param; j++) {
        n_gboy.pop[j].p = gboy.p;
        if (fabs(Dp_n(j)) > 0.) {
            delta(j) *= Dp_n(j);
//            delta(j) *= (0.1*gboy.p(j));
            n_gboy.pop[j].p(j) += delta(j);
        }
        else {
            delta(j) *= (0.1*gboy.p(j));
            n_gboy.pop[j].p(j) += delta(j);
        }
    }
    
    for(int j=0; j<n_param; j++) {
        //run the simulation
        run_simulation(simul_type, n_gboy.pop[j], nfiles, params, consts, data_num, folder, name, path_data, path_keys, materialfile);
        vec vnum = calcV(data_num, data_exp, nfiles, sizev);
        
        calcS(S, vnum, vnum0, j, delta);
    }
    return S;
}
    
    
    
    
} //namespace simcoon
