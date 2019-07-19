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

///@file identification.cpp
///@brief The main identification function
///@version 1.0

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/random.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/optimize.hpp>
#include <simcoon/Simulation/Identification/generation.hpp>
#include <simcoon/Simulation/Identification/methods.hpp>
#include <simcoon/Simulation/Identification/opti_data.hpp>
#include <simcoon/Simulation/Identification/doe.hpp>
#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/script.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
        
void run_identification(const std::string &simul_type, const int &n_param, const int &n_consts, const int &nfiles, const int &ngen, const int &aleaspace, int &apop, int &spop, const int &ngboys, const int &maxpop, const int &stationnarity_nb, const std::string &path_data, const std::string &path_keys, const std::string &path_results, const std::string &materialfile, const std::string &outputfile, const std::string &data_num_name, const double &probaMut, const double &pertu, const double &c, const double &p0, const double &lambdaLM) {

    std::string data_num_ext = data_num_name.substr(data_num_name.length()-4,data_num_name.length());
    std::string data_num_name_root = data_num_name.substr(0,data_num_name.length()-4); //to remove the extension
    
    cout << boost::filesystem::current_path().string() << endl;
    
    //Check if the required directories exist:
    if(!boost::filesystem::is_directory(path_data)) {
        cout << "error: the folder for the data, " << path_data << ", is not present" << endl;
        return;
    }
    if(!boost::filesystem::is_directory(path_keys)) {
        cout << "error: the folder for the keys, " << path_keys << ", is not present" << endl;
        return;
    }
    if(!boost::filesystem::is_directory(path_results)) {
        cout << "The folder for the results, " << path_results << ", is not present and has been created" << endl;
        boost::filesystem::create_directory(path_results);
    }
    
    //Check consistency of data
    if((aleaspace==0)||(aleaspace==1)) {
        if(maxpop > spop*n_param) {
            cout << "Please increase the mesh grid for the first generation (Space population) or reduce the max number population per subgeneration\n";
            exit(0);
        }
    }
    else if((aleaspace==2)||(aleaspace==3)) {
        if(maxpop > apop) {
            cout << "Please increase the Space population or reduce the max number population per subgeneration\n";
            exit(0);
        }
    }
    
    if(ngboys > maxpop) {
        cout << "Please increase the the max number population per subgeneration or reduce the number of gboys\n";
        exit(0);
    }
    
    ///Allow non-repetitive pseudo-random number generation
    srand(time(0));
    ofstream result;    ///Output stream, with parameters values and cost function

    //Define the parameters
    vector<parameters> params(n_param);  //vector of parameters
    vector<constants> consts(n_consts);  //vector of constants
    vec Dp = zeros(n_param);
    vec p = zeros(n_param);
    
    //Read the parameters and constants
    read_parameters(n_param, params);
    read_constants(n_consts, consts, nfiles);
    
    int idnumber=1;
    int id0=0;
    
    //Get the experimental data file
    string data_exp_folder="exp_data";
    if(!boost::filesystem::is_directory(data_exp_folder)) {
        cout << "The folder for the experimental data, " << data_exp_folder << ", is not present" << endl;
        return;
    }
    
    //Get the experimental data and build the exp vector, and get the size of vectors
    int sizev = 0;
    vector<opti_data> data_exp(nfiles);
    read_data_exp(nfiles, data_exp);
    for(int i=0; i<nfiles; i++) {
        data_exp[i].import(data_exp_folder);
        sizev += data_exp[i].ndata * data_exp[i].ninfo;
    }
    vec vexp = calcV(data_exp, data_exp, nfiles, sizev);
    
    //Get the weight data and build the weight vector
    vector<opti_data> data_weight(nfiles);
    Col<int> weight_types(3);
    vec weight_files = zeros(nfiles);
    vector<vec> weight_cols(nfiles);
    read_data_weights(nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);
    for(int i=0; i<nfiles; i++) {
        data_weight[i].import(data_exp_folder);
    }
    vec W = calcW(sizev, nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);
    
    //Get the data structures for the num data
    vector<opti_data> data_num(nfiles);
    read_data_num(nfiles, data_exp, data_num);
    vec vnum = zeros(sizev);   //num vector
    
    //Data structure has been created. Next is the generation of structures to compute cost function and associated derivatives
    mat S(sizev,n_param);
    Col<int> pb_col;
    pb_col.zeros(n_param);
    
    result.open(outputfile,  ios::out);
    result << "g" << "\t";
    result << "nindividual" << "\t";
    result << "cost" << "\t";
    for(int i=0; i<n_param; i++) {
        result << "p(" << i << ")" << "\t";
    }
    result << "\n";
    result.close();

    ///Optimization process
    ///Creation of the generation table
    vector<generation> gen(ngen+1);
    vector<generation> gboys(ngen+1);

    generation geninit;
    int g=0;

    gen[g].construct(maxpop, n_param, id0, lambdaLM);
    if(ngboys) {
        gboys[g].construct(ngboys, n_param, id0, lambdaLM);
    }
    gen_initialize(geninit, spop, apop, idnumber, aleaspace, n_param, params, lambdaLM);
    
    string data_num_folder = "num_data";
    if(!boost::filesystem::is_directory(data_num_folder)) {
        cout << "The folder for the numerical data, " << data_num_folder << ", is not present and has been created" << endl;
        boost::filesystem::create_directory(data_num_folder);
    }
        
    /// Run the simulations corresponding to each individual
    /// The simulation input files should be ready!
    for(int i=0; i<geninit.size(); i++) {
        
        run_simulation(simul_type, geninit.pop[i], nfiles, params, consts, data_num, data_num_folder, data_num_name, path_data, path_keys, materialfile);
        
        //Calculation of the cost function
        geninit.pop[i].cout = calc_cost(vexp, vnum, W, data_num, data_exp, nfiles, sizev);
    }
    
    //Classification of bests
    for(int i=0; i<maxpop; i++) {
        gen[0].pop[i]=geninit.pop[i];
    }
    gen[0].classify();

    result.open(outputfile,  ios::out | ios::app);
    for(int i=0; i<maxpop; i++) {
        result << 0 << "\t" << gen[0].pop[i].id << "\t" << gen[0].pop[i].cout << "\t";
        for(int j=0; j<n_param;j++) {
            result << gen[0].pop[i].p(j);
            if(j==n_param-1)
                result << "\n";
            else
                result << "\t";
        }
    }
    result.close();

    cout << "\nCost function (Best set of parameters)  = " << gen[0].pop[0].cout << "\n";

    ///Next step : Classify the best one, and compute the next generation!
    // Here we can choose the type of optimizer we want
    ///Creation of the generation table
    generation gensons(maxpop, n_param, id0);
    generation n_gboys(n_param, n_param, id0);
    
    //get the first gboys to be optimized via gradient_based
    for(int i=0; i<ngboys; i++) {
        gboys[0].pop[i] = gen[0].pop[i];
    }
    
    double costnm1 = 0.;
    double stationnarity = 1.E-12;		/// Stationnary stopping criteria (no more evolution of the cost function)
    
    std::vector<double> cost_gb_cost_n(ngboys);
    std::vector<vec> Dp_gb_n(ngboys);
    
    for(int i=0; i<ngboys; i++) {
        Dp_gb_n[i] = zeros(n_param);
    }
    //bool bad_des = false;
    int compt_des = 0;
    
    while((g<ngen)&&(compt_des < stationnarity_nb)) {
//    while(g<ngen) {
        
        costnm1 = gen[g].pop[0].cout;
        
        /// Run the simulations corresponding to each individual
        /// The simulation input files should be ready!
        if (maxpop > 1) {
            
            genetic(gen[g], gensons, idnumber, probaMut, pertu, params);
            ///prepare the individuals to run
            
            for(int i=0; i<gensons.size(); i++) {
                run_simulation(simul_type, gensons.pop[i], nfiles, params, consts, data_num, data_num_folder, data_num_name, path_data, path_keys, materialfile);
                //Calculation of the cost function
                gensons.pop[i].cout = calc_cost(vexp, vnum, W, data_num, data_exp, nfiles, sizev);
            }
            
        }
        for (int i=0; i<ngboys; i++) {
            
            cost_gb_cost_n[i] = gen[g].pop[i].cout;
            
            S = calc_sensi(gboys[g].pop[i], n_gboys, simul_type, nfiles, n_param, params, consts, vnum, data_num, data_exp, data_num_folder, data_num_name, path_data, path_keys, sizev, Dp_gb_n[i], materialfile);
            gboys[g].pop[i].cout = calcC(vexp, vnum, W);
            p = gboys[g].pop[i].p;
            ///Compute the parameters increment
            Dp = calcDp(S, vexp, vnum, W, p, params, gboys[g].pop[i].lambda, c, p0, n_param, pb_col);
            p += Dp;
            
            Dp_gb_n[i] = Dp;

            for(int j=0; j < n_param; j++) {
                if(p(j) > params[j].max_value)
                    p(j) = params[j].max_value;
                if(p(j) < params[j].min_value)
                    p(j) = params[j].min_value;
            }
            gboys[g].pop[i].p = p;
            
            if(gboys[g].pop[i].cout > cost_gb_cost_n[i]) {
                //bad_des = true;
                gboys[g].pop[i].lambda *= 3;
                gboys[g].pop[i].p = gen[g].pop[i].p;
                gboys[g].pop[i].cout = cost_gb_cost_n[i];
            }
            else if(gboys[g].pop[i].cout < cost_gb_cost_n[i]) {
                gboys[g].pop[i].lambda *= 0.5;
            }
                
        }
        
        ///Find the bests
        g++;
        find_best(gen[g], gboys[g], gen[g-1], gboys[g-1], gensons, maxpop, n_param, id0);
        write_results(result, outputfile, gen[g], g, maxpop, n_param);
        
        if(fabs(costnm1 - gen[g].pop[0].cout) < stationnarity) {
            compt_des++;
        }
        else
            compt_des = 0;
        
        cout << "Cost function (Best set of parameters) = " << gen[g].pop[0].cout << "\n";
        
        //Replace the parameters
        for (unsigned int k=0; k<params.size(); k++) {
            params[k].value = gen[g].pop[0].p(k);
        }
        
        //In the simulation run, make sure that we remove all the temporary files
        boost::filesystem::path path_to_remove(data_num_folder);
        for (boost::filesystem::directory_iterator end_dir_it, it(path_to_remove); it!=end_dir_it; ++it) {
            boost::filesystem::remove_all(it->path());
        }
        
        //Run the identified simulation and store results in the results folder
        run_simulation(simul_type, gen[g].pop[0], nfiles, params, consts, data_num, path_results, data_num_name, path_data, path_keys, materialfile);
        
        for (int i = 0; i<nfiles; i++) {
            string simulfile = path_results + "/" + data_num_name_root + "_" + to_string(gen[g].pop[0].id)  + "_" + to_string(i+1) + data_num_ext;
            string finalfile = path_results + "/" + data_num_name_root + "_" + to_string(i+1) + data_num_ext;
            boost::filesystem::rename(simulfile, finalfile);
        }
        
        copy_parameters(params, path_keys, path_results);
        apply_parameters(params, path_results);
    }
    
}

} //namespace simcoon
