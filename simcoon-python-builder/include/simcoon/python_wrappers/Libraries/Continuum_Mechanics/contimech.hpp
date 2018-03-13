#pragma once
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace simpy{
    
//This function returns the trace of the tensor v
double tr(const boost::python::numpy::ndarray &);

//This function returns the deviatoric part of v
boost::python::numpy::ndarray dev(const boost::python::numpy::ndarray &);

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress
double Mises_stress(const boost::python::numpy::ndarray &);

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
boost::python::numpy::ndarray eta_stress(const boost::python::numpy::ndarray &);

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains
double Mises_strain(const boost::python::numpy::ndarray &);

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
boost::python::numpy::ndarray eta_strain(const boost::python::numpy::ndarray &);

//This function transforms the strain Voigt vector into a 3*3 strain matrix
boost::python::numpy::ndarray v2t_strain(const boost::python::numpy::ndarray &);

//This function transforms a 3*3 strain matrix into a strain Voigt vector
boost::python::numpy::ndarray t2v_strain (const boost::python::numpy::ndarray &);

//This function transforms the stress Voigt vector into a 3*3 stress matrix
boost::python::numpy::ndarray v2t_stress(const boost::python::numpy::ndarray &);

//This function transforms a 3*3 stress matrix into a stress Voigt vector
boost::python::numpy::ndarray t2v_stress (const boost::python::numpy::ndarray &);

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const boost::python::numpy::ndarray &);

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const boost::python::numpy::ndarray &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const boost::python::numpy::ndarray &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const boost::python::numpy::ndarray &);

//This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)
double Macaulay_p(const double &);

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double Macaulay_n(const double &);

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double sign(const double &);

//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
boost::python::numpy::ndarray normal_ellipsoid(const double &, const double &, const double &, const double &, const double &);

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
boost::python::numpy::ndarray sigma_int(const boost::python::numpy::ndarray &, const double &, const double &, const double &, const double &, const double &);

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
boost::python::numpy::ndarray p_ikjl(const boost::python::numpy::ndarray &);
    
} //namespace simpy