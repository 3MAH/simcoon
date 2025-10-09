#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>
#include <assert.h>

#include <simcoon/Simulation/Identification/identification.hpp>
#include <simcoon/Simulation/Identification/constants.hpp>
#include <simcoon/Simulation/Identification/parameters.hpp>
#include <simcoon/Simulation/Identification/read.hpp>
#include <simcoon/Simulation/Identification/script.hpp>
#include <simcoon/Simulation/Identification/optimize.hpp>
#include <simcoon/python_wrappers/Libraries/Identification/identification.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {
    
//This function computes the identifcation of materials parameters for one/multiple homogeneous mixed thermomechanical loading experiment
void identification(const std::string &simul_type_py, const int &n_param, const int &n_consts, const int &nfiles, const int &ngen, const int &aleaspace, const int &pop_py, const int &ngboys, const int &maxpop, const std::string &path_data_py, const std::string &path_keys_py, const std::string &path_results_py, const std::string &materialfile_py, const std::string &outputfile_py) {
    
    int apop = 0;
    int spop = 0;
    
    if(aleaspace == 2)
        apop = pop_py;
    else if(aleaspace < 2)
        spop = pop_py;

    int station_nb = 6;
    double station_lim = 1.E-12;
    simcoon::run_identification(simul_type_py,n_param, n_consts, nfiles, ngen, aleaspace, apop, spop, ngboys, maxpop, station_nb, station_lim, path_data_py, path_keys_py, path_results_py, materialfile_py, outputfile_py);
}

py::list read_constants_py(const int &nconstants, const int &nfiles) {
    std::vector<simcoon::constants> consts(nconstants);
    py::list list_to_return = py::cast(consts);
    return list_to_return;
}
    
py::list read_parameters_py(const int &nparams) {
    std::vector<simcoon::parameters> params(nparams);
    py::list list_to_return = py::cast(params);
    return list_to_return;
}
    
//This function will copy the constant files
void copy_constants_py(const py::list &consts_py, const string &src_path, const string &dst_path) {

    std::vector<simcoon::constants> consts = consts_py.cast<std::vector<simcoon::constants>>();
    simcoon::copy_constants(consts, src_path, dst_path);
}

//This function will copy the parameters files
void copy_parameters_py(const py::list &params_py, const string &src_path, const string &dst_path) {
    
    std::vector<simcoon::parameters> params = params_py.cast<std::vector<simcoon::parameters>>();
    simcoon::copy_parameters(params, src_path, dst_path);
}
    
void apply_constants_py(const py::list &consts_py, const string &dst_path) {
    
    std::vector<simcoon::constants> consts = consts_py.cast<std::vector<simcoon::constants>>();
    simcoon::apply_constants(consts, dst_path);
}

void apply_parameters_py(const py::list &params_py, const string &dst_path) {
    
    std::vector<simcoon::parameters> params = params_py.cast<std::vector<simcoon::parameters>>();
    simcoon::apply_parameters(params, dst_path);
}
    
double calc_cost(const int &nfiles, const py::list &data_num_names_list) {

    //Get the data structures
    std::vector<simcoon::opti_data> data_exp(nfiles);
    std::vector<simcoon::opti_data> data_weight(nfiles);
    std::vector<simcoon::opti_data> data_num(nfiles);
    
    Col<int> weight_types(3);
    vec weight_files = zeros(nfiles);
    vector<vec> weight_cols(nfiles);
    
    simcoon::read_data_exp(nfiles, data_exp);
    simcoon::read_data_weights(nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);
    simcoon::read_data_num(nfiles, data_exp, data_num);
    
    /// Get the data vectors
    ///Import of the experimental data
    string data_exp_folder="exp_data";
    string data_num_folder="num_data";
    
    int sizev = 0;
    for(int i=0; i<nfiles;i++) {

        py::object item = data_num_names_list[i];
        std::string data_num_item = item.cast<std::string>();

        data_exp[i].import(data_exp_folder);
        data_weight[i].import(data_exp_folder);
        sizev += data_exp[i].ndata * data_exp[i].ninfo;
        
        data_num[i].name = data_num_item;
        data_num[i].import(data_num_folder);
    }
    
    ///Computation of the cost function
    vec vexp = simcoon::calcV(data_exp, data_exp, nfiles, sizev);
    vec vnum = simcoon::calcV(data_num, data_exp, nfiles, sizev);
    vec W = simcoon::calcW(sizev, nfiles, weight_types, weight_files, weight_cols, data_weight, data_exp);

    return simcoon::calcC(vexp, vnum, W);
}
    
    
} //namepsace simpy
