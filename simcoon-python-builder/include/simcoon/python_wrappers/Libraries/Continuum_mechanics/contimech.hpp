#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace simpy{
    
//This function returns the trace of the tensor v
double tr(const pybind11::array_t<double> &);

//This function returns the deviatoric part of v
pybind11::array_t<double> dev(const pybind11::array_t<double> &, const bool &copy=true);

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress
double Mises_stress(const pybind11::array_t<double> &);

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
pybind11::array_t<double> eta_stress(const pybind11::array_t<double> &, const bool &copy=true);

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains
double Mises_strain(const pybind11::array_t<double> &);

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
pybind11::array_t<double> eta_strain(const pybind11::array_t<double> &, const bool &copy=true);

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const pybind11::array_t<double> &);

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const pybind11::array_t<double> &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const pybind11::array_t<double> &);

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const pybind11::array_t<double> &);

//This function returns the value if it's positive, zero if it's negative (Macaulay brackets <>+)
double Macaulay_p(const double &);

//This function returns the value if it's negative, zero if it's positive (Macaulay brackets <>-)
double Macaulay_n(const double &);

//This function returns the sign of a double
double sign(const double &);

//Returns the normalized vector normal to an ellipsoid with semi-principal axes of length a1, a2, a3. The direction of the normalized vector is set by angles u
pybind11::array_t<double> normal_ellipsoid(const double &, const double &, const double &, const double &, const double &, const bool &copy=true);

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
pybind11::array_t<double> sigma_int(const pybind11::array_t<double> &, const double &, const double &, const double &, const double &, const double &, const bool &copy=true);

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
pybind11::array_t<double> p_ikjl(const pybind11::array_t<double> &, const bool &copy=true);
    
} //namespace simpy