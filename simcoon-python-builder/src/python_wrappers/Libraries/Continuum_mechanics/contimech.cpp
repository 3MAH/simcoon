
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/contimech.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//This function returns the trace of the tensor v
double tr(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::tr(v);
}

//This function returns the deviatoric part of v
py::array_t<double> dev(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec t = simcoon::dev(v);
    return carma::col_to_arr(t, copy);
}

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress
double Mises_stress(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::Mises_stress(v);
}

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
py::array_t<double> eta_stress(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec t = simcoon::eta_stress(v);
    return carma::col_to_arr(t, copy);
}

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains
double Mises_strain(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::Mises_strain(v);
}

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
py::array_t<double> eta_strain(const py::array_t<double> &input, const bool &copy) {
    vec v = carma::arr_to_col(input);
    vec t = simcoon::eta_strain(v);
    return carma::col_to_arr(t, copy);
}

//Returns the secoinput invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::J2_stress(v);
}

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::J2_strain(v);
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::J3_stress(v);
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const py::array_t<double> &input) {
    vec v = carma::arr_to_col(input);
    return simcoon::J3_strain(v);
}

//This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)
double Macaulay_p(const double &d) {
    return simcoon::Macaulay_p(d);
}

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double Macaulay_n(const double &d) {
    return simcoon::Macaulay_n(d);
}

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double sign(const double &d) {
    return simcoon::sign(d);
}

//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
py::array_t<double> normal_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3, const bool &copy) {
    vec t = simcoon::normal_ellipsoid(u,v,a1,a2,a3);
    return carma::col_to_arr(t, copy);
}

//Provides the curvature of an ellipsoid with semi-principal axes of length a1, a2, a3 at the angle u,v.
double  curvature_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3) {
    return simcoon::curvature_ellipsoid(u,v,a1,a2,a3);
}

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
py::array_t<double> sigma_int(const py::array_t<double> &input, const double &u, const double &v, const double &a1, const double &a2, const double &a3, const bool &copy) {
    vec sigma_in = carma::arr_to_col(input);
    vec t = simcoon::sigma_int(sigma_in,u,v,a1,a2,a3);
    return carma::col_to_arr(t, copy);
}

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
py::array_t<double> p_ikjl(const py::array_t<double> &normal, const bool &copy) {
    vec a = carma::arr_to_col(normal);
    mat t = simcoon::p_ikjl(a);
    return carma::mat_to_arr(t, copy);
}

py::array_t<double> auto_sym_dyadic(const py::array_t<double> &input, const bool &copy) {
    mat a = carma::arr_to_mat(input);
    mat c = simcoon::auto_sym_dyadic(a);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> sym_dyadic(const py::array_t<double> &a, const py::array_t<double> &b, const bool &copy) {
    mat a_cpp = carma::arr_to_mat(a);
    mat b_cpp = carma::arr_to_mat(b);    
    mat c = simcoon::sym_dyadic(a_cpp, b_cpp);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> auto_dyadic(const py::array_t<double> &input, const bool &copy) {
    mat a = carma::arr_to_mat(input);
    mat c = simcoon::auto_dyadic(a);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> dyadic_4vectors_sym(const py::array_t<double> &a, const py::array_t<double> &b, const std::string &conv, const bool &copy=true) {
    vec a_cpp = carma::arr_to_col(a);
    vec b_cpp = carma::arr_to_col(b);    
    mat c = simcoon::dyadic_4vectors_sym(a_cpp, b_cpp, conv);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> dyadic(const py::array_t<double> &a, const py::array_t<double> &b, const bool &copy) {
    mat a_cpp = carma::arr_to_mat(a);
    mat b_cpp = carma::arr_to_mat(b);    
    mat c = simcoon::dyadic(a_cpp, b_cpp);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> auto_dyadic_operator(const py::array_t<double> &input, const bool &copy) {
    mat a = carma::arr_to_mat(input);
    mat c = simcoon::auto_dyadic_operator(a);
    return carma::mat_to_arr(c, copy);
}

py::array_t<double> sym_dyadic_operator(const py::array_t<double> &a, const py::array_t<double> &b, const bool &copy) {
    mat a_cpp = carma::arr_to_mat(a);
    mat b_cpp = carma::arr_to_mat(b);    
    mat c = simcoon::dyasym_dyadic_operatordic(a_cpp, b_cpp);
    return carma::mat_to_arr(c, copy);
}

} //namepsace simpy