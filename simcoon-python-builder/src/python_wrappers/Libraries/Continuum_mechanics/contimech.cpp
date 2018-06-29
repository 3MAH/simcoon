
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/contimech.hpp>

namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function returns the trace of the tensor v
double tr(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::tr(v);
}

//This function returns the deviatoric part of v
bn::ndarray dev(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return vec2array(simcoon::dev(v));
}

//This function determines the Mises equivalent of a stress tensor, according to the Voigt convention for stress
double Mises_stress(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::Mises_stress(v);
}

//This function determines the strain flow (direction) from a stress tensor, according to the Voigt convention for strains
bn::ndarray eta_stress(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return vec2array(simcoon::eta_stress(v));
}

//This function determines the Mises equivalent of a strain tensor, according to the Voigt convention for strains
double Mises_strain(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::Mises_strain(v);
}

//This function determines the strain flow (direction) from a strain tensor, according to the Voigt convention for strains
bn::ndarray eta_strain(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return vec2array(simcoon::eta_strain(v));
}

//Returns the second invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J2_stress(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::J2_stress(v);
}

//Returns the second invariant of the deviatoric part of a second order strain tensor written as a Voigt vector
double J2_strain(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::J2_strain(v);
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_stress(const bn::ndarray &nd) {
    vec v = array2vec(nd);
    return simcoon::J3_stress(v);
}

//Returns the third invariant of the deviatoric part of a second order stress tensor written as a Voigt vector
double J3_strain(const bn::ndarray &nd) {
    vec v = array2vec(nd);
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
bn::ndarray normal_ellipsoid(const double &u, const double &v, const double &a1, const double &a2, const double &a3) {
    return vec2array(simcoon::normal_ellipsoid(u,v,a1,a2,a3));
}

//Returns the normal and tangent components of the stress vector in the normal direction n to an ellipsoid with axes a1, a2, a3. The direction of the normalized vector is set by angles u
bn::ndarray sigma_int(const bn::ndarray &nd, const double &a1, const double &a2, const double &a3, const double &u, const double &v) {
    vec sigma_in = array2vec(nd);
    return vec2array(simcoon::sigma_int(sigma_in,a1,a2,a3,u,v));
}

///This computes the Hill interfacial operator according to a normal a (see papers of Siredey and Entemeyer phD dissertation)
bn::ndarray p_ikjl(const bn::ndarray &nd) {
    vec a = array2vec(nd);
    return mat2array(simcoon::p_ikjl(a));
}

} //namepsace simpy