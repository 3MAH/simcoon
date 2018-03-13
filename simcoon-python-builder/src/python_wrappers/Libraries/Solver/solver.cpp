
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/solver.hpp>
#include <simcoon/python_wrappers/Libraries/Solver/solver.hpp>

namespace bn = boost::python::numpy;
namespace bp = boost::python;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function computes the response of materials for an homogeneous mixed thermomechanical loading path    
void solver(const bp::str &umat_name_py, const bn::ndarray &props_py, const int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const bp::str &path_data_py, const bp::str &path_results_py, const bp::str &pathfile_py, const bp::str &outputfile_py) {
    
    vec props = array2vec(props_py);
    
    std::string umat_name = bp::extract<std::string>(umat_name_py);
    std::string path_data = bp::extract<std::string>(path_data_py);
    std::string path_results = bp::extract<std::string>(path_results_py);
    std::string pathfile = bp::extract<std::string>(pathfile_py);
    std::string outputfile = bp::extract<std::string>(outputfile_py);
    
    double div_tnew_dt_solver = 0.5;
    double mul_tnew_dt_solver = 2.;
    int miniter_solver = 10;
    int maxiter_solver = 100;
    int inforce_solver = 1;
    double precision_solver = 1.E-6;
    double lambda_solver = 10000.;
    
    simcoon::solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, path_results, pathfile, outputfile);
}

} //namepsace simpy