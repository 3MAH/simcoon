
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>
#include <simcoon/python_wrappers/Libraries/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//Eshelby tensor for a sphere
py::array_t<double> Eshelby_sphere(const double &nu, const bool &copy) {
    mat m = simcoon::Eshelby_sphere(nu);
    return carma::mat_to_arr(m, copy);
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
py::array_t<double> Eshelby_cylinder(const double &nu, const bool &copy) {
    return carma::mat_to_arr(simcoon::Eshelby_cylinder(nu));
}

//	Eshelby tensor determination. The prolate shape is oriented in such a way that the axis direction is the 1 direction. a1>a2=a3 here
py::array_t<double> Eshelby_prolate(const double &nu, const double &ar, const bool &copy) {
    mat m = simcoon::Eshelby_prolate(nu,ar);
    return carma::mat_to_arr(m, copy);    
}

//	Eshelby tensor determination. The oblate shape is oriented in such a way that the axis direction is the 1 direction. a1<a2=a3 here
py::array_t<double> Eshelby_oblate(const double &nu, const double &ar, const bool &copy) {
    mat m = simcoon::Eshelby_oblate(nu,ar);
    return carma::mat_to_arr(m, copy);    
}

//Numerical Eshelby tensor determination
py::array_t<double> Eshelby(const py::array_t<double> &L, const double &a1, const double &a2, const double &a3, const int &mp, const int &np, const bool &copy) {
    mat L_cpp = carma::arr_to_mat(L);
    mat m = simcoon::Eshelby(L_cpp, a1, a2, a3, mp, np);
    return carma::mat_to_arr(m, copy);    
}
    
//Numerical Hill Interaction tensor determination
py::array_t<double> T_II(const py::array_t<double> &L, const double &a1, const double &a2, const double &a3, const int &mp, const int &np, const bool &copy) {
    mat L_cpp = carma::arr_to_mat(L);
    mat m = simcoon::T_II(L_cpp, a1, a2, a3, mp, np);
    return carma::mat_to_arr(m, copy);
}

} //namepsace simpy
