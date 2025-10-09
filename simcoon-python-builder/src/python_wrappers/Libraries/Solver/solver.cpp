#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>

#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/solver.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path    
void solver(const std::string &umat_name_py, const py::array_t<double> &props_py, const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const std::string &path_data_py, const std::string &path_results_py, const std::string &pathfile_py, const std::string &outputfile_py) {
    
    vec props = simpy::arr_to_col(props_py);
    
    double div_tnew_dt_solver = 0.5;
    double mul_tnew_dt_solver = 2.;
    int miniter_solver = 10;
    int maxiter_solver = 100;
    int inforce_solver = 1;
    double precision_solver = 1.E-6;
    double lambda_solver = 10000.;
    
    simcoon::solver(umat_name_py, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, corate_type, div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data_py, path_results_py, pathfile_py, outputfile_py);
}

} //namepsace simpy
