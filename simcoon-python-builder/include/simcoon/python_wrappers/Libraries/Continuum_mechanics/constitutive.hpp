#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{

//Returns the fourth order identity tensor written in Voigt notation Ireal
pybind11::array_t<double> Ireal(const bool &copy=true);

//Returns the volumic of the identity tensor Ireal written in Voigt notation
pybind11::array_t<double> Ivol(const bool &copy=true);

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
pybind11::array_t<double> Idev(const bool &copy=true);

//Returns the fourth order identity tensor Iˆ written in Voigt notation
pybind11::array_t<double> Ireal2(const bool &copy=true);

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
pybind11::array_t<double> Idev2(const bool &copy=true);

//Returns the expansion vector
pybind11::array_t<double> Ith(const bool &copy=true);

//Returns the stress 2 strain operator
pybind11::array_t<double> Ir2(const bool &copy=true);

//Returns the strain 2 stress operator
pybind11::array_t<double> Ir05(const bool &copy=true);

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
pybind11::array_t<double> L_iso(const double &, const double &, const std::string &, const bool &copy=true);

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
pybind11::array_t<double> M_iso(const double &, const double &, const std::string &, const bool &copy=true);

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
pybind11::array_t<double> L_cubic(const double &, const double &, const double &, const std::string &, const bool &copy=true);

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
pybind11::array_t<double> M_cubic(const double &, const double &, const double &, const std::string &, const bool &copy=true);

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
pybind11::array_t<double> L_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string &, const bool &copy=true);

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
pybind11::array_t<double> M_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string &, const bool &copy=true);

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
pybind11::array_t<double> L_isotrans(const double &, const double &, const double &, const double &, const double &, const int &, const bool &copy=true);

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
pybind11::array_t<double> M_isotrans(const double &, const double &, const double &, const double &, const double &, const int &, const bool &copy=true);

//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
pybind11::array_t<double> H_iso(const double &, const double &, const bool &copy=true);
    
} //namespace simpy
