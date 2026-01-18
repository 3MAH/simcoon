
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
    // Accept either (6,) or (6,1) shaped Voigt vectors from Python.
    py::array_t<double> in = input;
    if ((int)input.ndim() == 2 && input.shape(1) == 1) {
        in = input.attr("ravel")().cast<py::array_t<double>>();
    }
    vec v = carma::arr_to_col(in);
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
    // Accept either (6,) or (6,1) shaped Voigt vectors from Python.
    py::array_t<double> in = input;
    if ((int)input.ndim() == 2 && input.shape(1) == 1) {
        in = input.attr("ravel")().cast<py::array_t<double>>();
    }
    vec v = carma::arr_to_col(in);
    mat m = simcoon::v2t_stress(v);
    return carma::mat_to_arr(m, copy);    
}

//This function transforms a 3*3 stress matrix into a stress Voigt vector
py::array_t<double> t2v_stress (const py::array_t<double> &input, const bool &copy) {
    // Accept either a 3x3 matrix or a (6,) / (6,1) Voigt-style input
    // If the input is a 2D array with shape (6,1) it may be a column;
    // ensure we convert matrices to arma::mat and then return a 1D Voigt array.
    if ((int)input.ndim() == 2 && input.shape(0) == 6 && input.shape(1) == 1) {
        // Input is a column Voigt vector accidentally passed here; convert to 1D and return
        py::array_t<double> col = input.attr("ravel")().cast<py::array_t<double>>();
        // Return squeezed column as 1D
        return col;
    }

    mat m = carma::arr_to_mat(input);
    vec v = simcoon::t2v_stress(m);
    // Ensure the returned numpy object is a 1D array (6,)
    return carma::col_to_arr(v, copy);
}

} //namepsace simpy