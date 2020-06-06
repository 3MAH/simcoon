
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/constitutive.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//Returns the fourth order identity tensor written in Voigt notation Ireal
bn::ndarray Ireal() {
    mat m = simcoon::Ireal();
    return mat2array(m);
}

//Returns the volumic of the identity tensor Ireal written in Voigt notation
bn::ndarray Ivol() {
    mat m = simcoon::Ivol();
    return mat2array(m);
}

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
bn::ndarray Idev() {
    mat m = simcoon::Idev();
    return mat2array(m);
}

//Returns the fourth order identity tensor Iˆ written in Voigt notation
bn::ndarray Ireal2() {
    mat m = simcoon::Ireal2();
    return mat2array(m);
}

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
bn::ndarray Idev2() {
    mat m = simcoon::Idev2();
    return mat2array(m);
}

//Returns the expansion vector
bn::ndarray Ith() {
    vec v = simcoon::Ith();
    return vec2array(v);
}

//Returns the stress 2 strain operator
bn::ndarray Ir2() {
    vec v = simcoon::Ir2();
    return vec2array(v);
}

//Returns the strain 2 stress operator
bn::ndarray Ir05() {
    vec v = simcoon::Ir05();
    return vec2array(v);
}

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
bn::ndarray L_iso(const double &C1, const double &C2, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_iso(C1,C2,conv_py);
    return mat2array(m);
}

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
bn::ndarray M_iso(const double &C1, const double &C2, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_iso(C1,C2,conv_py);
    return mat2array(m);
}

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44 ('Cii'), or E, nu and G ('EnuG')
bn::ndarray L_cubic(const double &C11, const double &C12, const double &C44, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_cubic(C11,C12,C44,conv_py);
    return mat2array(m);
}

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44 ('Cii'), or E, nu and G ('EnuG')
bn::ndarray M_cubic(const double &C11, const double &C12, const double &C44, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_cubic(C11,C12,C44,conv_py);
    return mat2array(m);
}

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
//C11,C12,C13,C22,C23,C13,C44,C55,C66 ('Cii')
//E1,E2,E3,nu12,nu13,nu23,G12,G13,G23 ('EnuG')
bn::ndarray L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_ortho(C11,C12,C13,C22,C23,C33,C44,C55,C66,conv_py);
    return mat2array(m);
}

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
//C11,C12,C13,C22,C23,C13,C44,C55,C66 ('Cii')
//E1,E2,E3,nu12,nu13,nu23,G12,G13,G23 ('EnuG')
bn::ndarray M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string &conv_py){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_ortho(C11,C12,C13,C22,C23,C33,C44,C55,C66,conv_py);
    return mat2array(m);
}

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
bn::ndarray L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
    mat m = simcoon::L_isotrans(EL,ET,nuTL,nuTT,GLT,axis);
    return mat2array(m);
}

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
bn::ndarray M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis){
    mat m = simcoon::M_isotrans(EL,ET,nuTL,nuTT,GLT,axis);
    return mat2array(m);
}
//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
bn::ndarray H_iso(const double &etaB, const double &etaS){
    mat m = simcoon::H_iso(etaB,etaS);
    return mat2array(m);
}

} //namepsace simpy
