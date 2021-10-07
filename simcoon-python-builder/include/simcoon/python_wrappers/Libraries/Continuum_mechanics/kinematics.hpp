#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
boost::python::numpy::ndarray G_UdX(const boost::python::numpy::ndarray &);

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
boost::python::numpy::ndarray G_Udx(const boost::python::numpy::ndarray &);

//This function computes the Right Cauchy-Green C
boost::python::numpy::ndarray R_Cauchy_Green(const boost::python::numpy::ndarray &);

//This function computes the Left Cauchy-Green B
boost::python::numpy::ndarray L_Cauchy_Green(const boost::python::numpy::ndarray &);

//This function computes the common Right (or Left) Cauchy-Green invariants
boost::python::numpy::ndarray Inv_X(const boost::python::numpy::ndarray &);

//This function computes the Cauchy deformation tensor c
boost::python::numpy::ndarray Cauchy(const boost::python::numpy::ndarray &);

//This function computes the Green-Lagrange finite strain tensor E
boost::python::numpy::ndarray Green_Lagrange(const boost::python::numpy::ndarray &);

//This function computes the Euler-Almansi finite strain tensor A
boost::python::numpy::ndarray Euler_Almansi(const boost::python::numpy::ndarray &);

//This function computes the logarithmic strain ln[V] = 1/2 ln[b] (b is the left Cauchy-Green Tensor)
boost::python::numpy::ndarray Log_strain(const boost::python::numpy::ndarray &);

//This function computes the velocity difference (F0,F1,DTime)
boost::python::numpy::ndarray finite_L(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);

//This function computes the spin tensor W (correspond to Jaumann rate) (F0,F1,DTime)
boost::python::numpy::ndarray finite_W(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);

//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the RU or VR polar decomposition of the deformation gradient F (F0,F1,DTime)
boost::python::numpy::ndarray finite_Omega(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);

//This function computes the deformation rate D (F0,F1,DTime)
boost::python::numpy::ndarray finite_D(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);

//This function computes the increment of finite rotation (Omega0, Omega1, DTime)
boost::python::numpy::ndarray finite_DQ(const boost::python::numpy::ndarray &, const boost::python::numpy::ndarray &, const double &);
    
} //namespace simpy
