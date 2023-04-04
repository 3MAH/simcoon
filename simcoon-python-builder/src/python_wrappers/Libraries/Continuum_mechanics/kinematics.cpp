
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//Provides the transformation gradient, from the Green-Lagrange strain and the rotation. 
py::array_t<double> ER_to_F(const py::array_t<double> &E, const py::array_t<double> &R, const bool &copy) {
    mat E_cpp = carma::arr_to_mat(E);
    mat R_cpp = carma::arr_to_mat(R);
    mat F = simcoon::ER_to_F(E_cpp, R_cpp);    
    return carma::mat_to_arr(F, copy);
}

//Provides the transformation gradient, from the Green-Lagrange strain and the rotation. 
py::array_t<double> eR_to_F(const py::array_t<double> &e, const py::array_t<double> &R, const bool &copy) {
    mat e_cpp = carma::arr_to_mat(e);
    mat R_cpp = carma::arr_to_mat(R);
    mat F = simcoon::eR_to_F(e_cpp, R_cpp);    
    return carma::mat_to_arr(F, copy);
}

//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
py::array_t<double> G_UdX(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat G_UdX = simcoon::G_UdX(F_cpp);
    return carma::mat_to_arr(G_UdX, copy);
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
py::array_t<double> G_Udx(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat G_Udx = simcoon::G_Udx(F_cpp);
    return carma::mat_to_arr(G_Udx, copy);
}

//This function computes the Right Cauchy-Green C
py::array_t<double> R_Cauchy_Green(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat R_Cauchy_Green = simcoon::R_Cauchy_Green(F_cpp);
    return carma::mat_to_arr(R_Cauchy_Green, copy);
}

//This function computes the Left Cauchy-Green B
py::array_t<double> L_Cauchy_Green(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat L_Cauchy_Green = simcoon::L_Cauchy_Green(F_cpp);
    return carma::mat_to_arr(L_Cauchy_Green, copy);
}

//Provides the RU decomposition of the transformation gradient F
py::tuple RU_decomposition(const pybind11::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat R_cpp = zeros(3,3);
    mat U_cpp = zeros(3,3);
    simcoon::RU_decomposition(R_cpp,U_cpp,F_cpp);
    return py::make_tuple(carma::mat_to_arr(R_cpp, copy), carma::mat_to_arr(U_cpp, copy));
}

//Provides the VR decomposition of the transformation gradient F
py::tuple VR_decomposition(const pybind11::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat V_cpp = zeros(3,3);
    mat R_cpp = zeros(3,3);    
    simcoon::VR_decomposition(V_cpp,R_cpp,F_cpp);
    return py::make_tuple(carma::mat_to_arr(V_cpp, copy), carma::mat_to_arr(R_cpp, copy));
}

//This function computes the common Right (or Left) Cauchy-Green invariants
py::array_t<double> Inv_X(const py::array_t<double> &input, const bool &copy) {
    mat X = carma::arr_to_mat(input);
    vec Inv_X = simcoon::Inv_X(X);
    return carma::col_to_arr(Inv_X, copy);
}

//This function computes the Cauchy deformation tensor c
py::array_t<double> Cauchy(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat Cauchy = simcoon::Cauchy(F_cpp);
    return carma::mat_to_arr(Cauchy, copy);
}

//This function computes the Green-Lagrange finite strain tensor E
py::array_t<double> Green_Lagrange(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat Green_Lagrange = simcoon::Green_Lagrange(F_cpp);
    return carma::mat_to_arr(Green_Lagrange, copy);
}

//This function computes the Euler-Almansi finite strain tensor A
py::array_t<double> Euler_Almansi(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat Euler_Almansi = simcoon::Euler_Almansi(F_cpp);
    return carma::mat_to_arr(Euler_Almansi, copy);
}

//This function computes the logarithmic strain ln[V] = 1/2 ln[b] (b is the left Cauchy-Green Tensor)
py::array_t<double> Log_strain(const py::array_t<double> &F, const bool &copy) {
    mat F_cpp = carma::arr_to_mat(F);
    mat Log_strain = simcoon::Log_strain(F_cpp);
    return carma::mat_to_arr(Log_strain, copy);
}

//This function computes the velocity difference
py::array_t<double> finite_L(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat finite_L = simcoon::finite_L(F0_cpp,F1_cpp,DTime);
    return carma::mat_to_arr(finite_L, copy);
}

//This function computes the deformation rate D
py::array_t<double> finite_D(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat finite_D = simcoon::finite_D(F0_cpp,F1_cpp,DTime);
    return carma::mat_to_arr(finite_D, copy);
}

//This function computes the spin tensor W (correspond to Jaumann rate)
py::array_t<double> finite_W(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat finite_W = simcoon::finite_W(F0_cpp,F1_cpp,DTime);
    return carma::mat_to_arr(finite_W, copy);
}

//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the polar decomposition of the deformation gradient F
py::array_t<double> finite_Omega(const py::array_t<double> &F0, const py::array_t<double> &F1, const double &DTime, const bool &copy) {
    mat F0_cpp = carma::arr_to_mat(F0);
    mat F1_cpp = carma::arr_to_mat(F1);
    mat finite_Omega = simcoon::finite_Omega(F0_cpp,F1_cpp,DTime);
    return carma::mat_to_arr(finite_Omega, copy);
}

//This function computes the increment of finite rotation
py::array_t<double> finite_DQ(const py::array_t<double> &Omega0, const py::array_t<double> &Omega1, const double &DTime, const bool &copy) {
    mat Omega0_cpp = carma::arr_to_mat(Omega0);
    mat Omega1_cpp = carma::arr_to_mat(Omega1);
    mat finite_DQ = simcoon::finite_DQ(Omega0_cpp,Omega1_cpp,DTime);
    return carma::mat_to_arr(finite_DQ, copy);
}

} //namepsace simpy
