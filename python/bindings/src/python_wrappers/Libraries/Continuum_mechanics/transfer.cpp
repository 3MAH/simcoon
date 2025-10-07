
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/transfer.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function transforms the strain Voigt vector into a 3*3 strain matrix
py::array_t<double> v2t_strain(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    mat m = simcoon::v2t_strain(v);
    return carma::mat_to_arr(m, copy);
}

//This function transforms a 3*3 strain matrix into a strain Voigt vector
py::array_t<double> t2v_strain (const py::array_t<double> &input, const bool &copy) {
    mat m = carma::arr_to_mat(input);
    vec v = simcoon::t2v_strain(m);
    return carma::col_to_arr(v, copy);
}

//This function transforms the stress Voigt vector into a 3*3 stress matrix
py::array_t<double> v2t_stress(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    mat m = simcoon::v2t_stress(v);
    return carma::mat_to_arr(m, copy);    
}

//This function transforms a 3*3 stress matrix into a stress Voigt vector
py::array_t<double> t2v_stress (const py::array_t<double> &input, const bool &copy) {
    mat m = carma::arr_to_mat(input);
    vec v = simcoon::t2v_stress(m);
    return carma::mat_to_arr(m, copy);    
}

} //namepsace simpy