
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/objective_rates.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function computes the logarithmic strain velocity and the logarithmic spin, along with the correct rotation increment
py::tuple logarithmic(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat DR = zeros(3,3);
    mat D = zeros(3,3);
    mat Omega = zeros(3,3);
    simcoon::logarithmic(DR, D, Omega, DTime, F0_cpp, F1_cpp);
    return py::make_tuple(carma::mat_to_arr(D, copy), carma::mat_to_arr(DR, copy), carma::mat_to_arr(Omega, copy));
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
py::array_t<double> Delta_log_strain(const py::array_t<double> &D, const py::array_t<double> &Omega, const double &DTime, const bool &copy) {
    mat D_cpp = carma::arr_to_mat(D);
    mat Omega_cpp = carma::arr_to_mat(Omega);
    mat Delta_log_strain = simcoon::Delta_log_strain(D_cpp, Omega_cpp, DTime);
    return carma::mat_to_arr(Delta_log_strain, copy);
}
    
} //namepsace simpy
