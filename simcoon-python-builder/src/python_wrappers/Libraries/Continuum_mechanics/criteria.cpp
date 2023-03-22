
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/criteria.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;
namespace simpy {

//This function returns the Prager equivalent stress.
double Prager_stress(const py::array_t<double> &nd, const double &b, const double &n) {
    vec v = carma::arr_to_col(nd);
    return simcoon::Prager_stress(v,b,n);
}

//This function returns the derivative of the Prager equivalent stress.
py::array_t<double> dPrager_stress(const py::array_t<double> &nd, const double &b, const double &n, const bool &copy) {
    vec v = carma::arr_to_col(nd);
    vec t = simcoon::dPrager_stress(v,b,n);
    return carma::col_to_arr(t, copy);
}

//This function returns the Tresca equivalent stress.
double Tresca_stress(const py::array_t<double> &nd) {
    vec v = carma::arr_to_col(nd);
    return simcoon::Tresca_stress(v);
}

//This function returns the derivative of the Tresca equivalent stress.
py::array_t<double> dTresca_stress(const py::array_t<double> &nd, const bool &copy) {
    vec v = carma::arr_to_col(nd);
    vec t = simcoon::dTresca_stress(v);
    return carma::col_to_arr(t, copy);
}

//Provides an anisotropic configurational tensor P in the Voigt format (6x6 numpy array), given its vector representation
py::array_t<double> P_ani(const py::array_t<double> &nparams, const bool &copy) {
    mat m = carma::arr_to_mat(nparams);
    mat t = simcoon::P_ani(m);
    return carma::mat_to_arr(t, copy);
}

//Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion in the Voigt format (6x6 numpy array), given its vector representation
py::array_t<double> P_hill(const py::array_t<double> &nparams, const bool &copy) {
    mat m = carma::arr_to_mat(nparams);
    mat t = simcoon::P_hill(m);
    return carma::mat_to_arr(t, copy);
}

//This function computes the selected equivalent stress function
double Eq_stress(const py::array_t<double> &nd, const string &eq_type, const py::array_t<double> &ndparam) {
    vec v = carma::arr_to_col(nd);
    vec param = carma::arr_to_col(ndparam);
    return simcoon::Eq_stress(v,eq_type,param);
}

//This function computes the deriavtive of the selected equivalent stress function
py::array_t<double> dEq_stress(const py::array_t<double> &nd, const string &eq_type, const py::array_t<double> &ndparam, const bool &copy) {
    vec v = carma::arr_to_col(nd);
    vec param = carma::arr_to_col(ndparam);
    vec t = simcoon::dEq_stress(v,eq_type,param);
    return carma::col_to_arr(t, copy);
}

} //namepsace simpy