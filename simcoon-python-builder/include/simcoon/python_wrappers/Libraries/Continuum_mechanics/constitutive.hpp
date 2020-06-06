#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//Returns the fourth order identity tensor written in Voigt notation Ireal
boost::python::numpy::ndarray Ireal();

//Returns the volumic of the identity tensor Ireal written in Voigt notation
boost::python::numpy::ndarray Ivol();

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
boost::python::numpy::ndarray Idev();

//Returns the fourth order identity tensor Iˆ written in Voigt notation
boost::python::numpy::ndarray Ireal2();

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
boost::python::numpy::ndarray Idev2();

//Returns the expansion vector
boost::python::numpy::ndarray Ith();

//Returns the stress 2 strain operator
boost::python::numpy::ndarray Ir2();

//Returns the strain 2 stress operator
boost::python::numpy::ndarray Ir05();

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
boost::python::numpy::ndarray L_iso(const double &, const double &, const std::string &);

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
boost::python::numpy::ndarray M_iso(const double &, const double &, const std::string &);

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
boost::python::numpy::ndarray L_cubic(const double &, const double &, const double &, const std::string &);

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44.
boost::python::numpy::ndarray M_cubic(const double &, const double &, const double &, const std::string &);

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
boost::python::numpy::ndarray L_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string &);

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
boost::python::numpy::ndarray M_ortho(const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const double &, const std::string &);

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
boost::python::numpy::ndarray L_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
boost::python::numpy::ndarray M_isotrans(const double &, const double &, const double &, const double &, const double &, const int &);

//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
boost::python::numpy::ndarray H_iso(const double &, const double &);
    
} //namespace simpy
