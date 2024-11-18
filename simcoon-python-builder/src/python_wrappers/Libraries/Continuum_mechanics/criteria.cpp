
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
double Prager_stress(const py::array_t<double> &input, const py::array_t<double> &props) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);    
    return simcoon::Prager_stress(v,props_cpp(0),props_cpp(1));
}

//This function returns the derivative of the Prager equivalent stress.
py::array_t<double> dPrager_stress(const py::array_t<double> &input, const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);
    vec t = simcoon::dPrager_stress(v,props_cpp(0),props_cpp(1));
    return carma::col_to_arr(t, copy);
}

//This function returns the Tresca equivalent stress.
double Tresca_stress(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::Tresca_stress(v);
}

//This function returns the derivative of the Tresca equivalent stress.
py::array_t<double> dTresca_stress(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec t = simcoon::dTresca_stress(v);
    return carma::col_to_arr(t, copy);
}

//Provides an anisotropic configurational tensor P in the Voigt format (6x6 numpy array), given its vector representation
py::array_t<double> P_Ani(const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(props);
    mat t = simcoon::P_Ani(v);
    return carma::mat_to_arr(t, copy);
}

//Provides an anisotropic configurational tensor considering the quadratic Hill yield criterion in the Voigt format (6x6 numpy array), given its vector representation
py::array_t<double> P_Hill(const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(props);
    mat t = simcoon::P_Hill(v);
    return carma::mat_to_arr(t, copy);
}

//This function returns the Hill equivalent stress.
double Hill_stress(const py::array_t<double> &input, const py::array_t<double> &props) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);
    return simcoon::Hill_stress(v, props_cpp);
}

//This function returns the derivative of the Hill equivalent stress.
py::array_t<double> dHill_stress(const py::array_t<double> &input, const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);    
    vec t = simcoon::dHill_stress(v, props_cpp);
    return carma::col_to_arr(t, copy);
}

//This function returns the anisotropic equivalent stress.
double Ani_stress(const py::array_t<double> &input, const py::array_t<double> &props) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);    
    return simcoon::Ani_stress(v, props_cpp);
}

//This function returns the derivative of the anisotropic equivalent stress.
py::array_t<double> dAni_stress(const py::array_t<double> &input, const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);      
    vec t = simcoon::dAni_stress(v, props_cpp);
    return carma::col_to_arr(t, copy);
}

//This function returns the DFA equivalent stress.
double DFA_stress(const py::array_t<double> &input, const py::array_t<double> &props) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);
    return simcoon::DFA_stress(v, props_cpp);
}

//This function returns the derivative of the DFA equivalent stress.
py::array_t<double> dDFA_stress(const py::array_t<double> &input, const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec props_cpp = carma::arr_to_col(props);
    vec t = simcoon::dDFA_stress(v, props_cpp);
    return carma::col_to_arr(t, copy);
}

//This function computes the selected equivalent stress function
double Eq_stress(const py::array_t<double> &input, const string &criteria, const py::array_t<double> &props) {
    vec v = carma::arr_to_col(input);
    vec param = carma::arr_to_col(props);
    return simcoon::Eq_stress(v,criteria,param);
}

//This function computes the deriavtive of the selected equivalent stress function
py::array_t<double> dEq_stress(const py::array_t<double> &input, const string &criteria, const py::array_t<double> &props, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec param = carma::arr_to_col(props);
    vec t = simcoon::dEq_stress(v,criteria,param);
    return carma::col_to_arr(t, copy);
}

} //namepsace simpy
