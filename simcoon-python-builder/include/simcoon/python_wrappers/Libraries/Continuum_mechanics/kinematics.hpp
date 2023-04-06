#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Provides the transformation gradient, from the Green-Lagrange strain and the rotation.   
pybind11::array_t<double> ER_to_F(const pybind11::array_t<double> &E, const pybind11::array_t<double> &R, const bool &copy=true);

//Provides the transformation gradient, from the logarithmic strain and the rotation.
pybind11::array_t<double> eR_to_F(const pybind11::array_t<double> &e, const pybind11::array_t<double> &R, const bool &copy=true);

//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
pybind11::array_t<double> G_UdX(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
pybind11::array_t<double> G_Udx(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the Right Cauchy-Green C
pybind11::array_t<double> R_Cauchy_Green(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the Left Cauchy-Green B
pybind11::array_t<double> L_Cauchy_Green(const pybind11::array_t<double> &F, const bool &copy=true);

//Provides the RU decomposition of the transformation gradient F
pybind11::tuple RU_decomposition(const pybind11::array_t<double> &F, const bool &copy=true);

//Provides the VR decomposition of the transformation gradient F
pybind11::tuple VR_decomposition(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the common Right (or Left) Cauchy-Green invariants
pybind11::array_t<double> Inv_X(const pybind11::array_t<double> &input, const bool &copy=true);

//This function computes the Cauchy deformation tensor c from the transformation gradient F
pybind11::array_t<double> Cauchy(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the Green-Lagrange finite strain tensor E
pybind11::array_t<double> Green_Lagrange(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the Euler-Almansi finite strain tensor A
pybind11::array_t<double> Euler_Almansi(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the logarithmic strain ln[V] = 1/2 ln[B] (B is the left Cauchy-Green Tensor)
pybind11::array_t<double> Log_strain(const pybind11::array_t<double> &F, const bool &copy=true);

//This function computes the velocity difference (F0,F1,DTime)
pybind11::array_t<double> finite_L(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);

//This function computes the deformation rate D (F0,F1,DTime)
pybind11::array_t<double> finite_D(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);

//This function computes the spin tensor W (correspond to Jaumann rate) (F0,F1,DTime)
pybind11::array_t<double> finite_W(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);

//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the RU or VR polar decomposition of the deformation gradient F (F0,F1,DTime)
pybind11::array_t<double> finite_Omega(const pybind11::array_t<double> &F0, const pybind11::array_t<double> &F1, const double &DTime, const bool &copy=true);

//This function computes the increment of finite rotation (Omega0, Omega1, DTime)
pybind11::array_t<double> finite_DQ(const pybind11::array_t<double> &Omega0, const pybind11::array_t<double> &Omega1, const double &DTime, const bool &copy=true);
    
} //namespace simpy
