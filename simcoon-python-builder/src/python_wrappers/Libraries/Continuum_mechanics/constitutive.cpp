
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <carma>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/constitutive.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

//Returns the fourth order identity tensor written in Voigt notation Ireal
py::array_t<double> Ireal(const bool &copy) {
    mat m = simcoon::Ireal();
    return carma::mat_to_arr(m, copy);
}

//Returns the volumic of the identity tensor Ireal written in Voigt notation
py::array_t<double> Ivol(const bool &copy) {
    mat m = simcoon::Ivol();
    return carma::mat_to_arr(m, copy);
}

//Returns the deviatoric of the identity tensor Ireal written in Voigt notation
py::array_t<double> Idev(const bool &copy) {
    mat m = simcoon::Idev();
    return carma::mat_to_arr(m, copy);
}

//Returns the fourth order identity tensor Iˆ written in Voigt notation
py::array_t<double> Ireal2(const bool &copy) {
    mat m = simcoon::Ireal2();
    return carma::mat_to_arr(m, copy);
}

//Returns the deviatoric of the identity tensor Iˆ written in Voigt notation
py::array_t<double> Idev2(const bool &copy) {
    mat m = simcoon::Idev2();
    return carma::mat_to_arr(m, copy);
}

//Returns the expansion vector
py::array_t<double> Ith(const bool &copy) {
    vec v = simcoon::Ith();
    return carma::col_to_arr(v, copy);
}

//Returns the stress 2 strain operator
py::array_t<double> Ir2(const bool &copy) {
    vec v = simcoon::Ir2();
    return carma::col_to_arr(v, copy);
}

//Returns the strain 2 stress operator
py::array_t<double> Ir05(const bool &copy) {
    vec v = simcoon::Ir05();
    return carma::col_to_arr(v, copy);
}

//Provides the elastic stiffness tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
// ‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
py::array_t<double> L_iso(const double &C1, const double &C2, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_iso(C1,C2,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Provides the elastic compliance tensor for an isotropic material.
//The two first arguments are a couple of Lamé coefficients. The third argument specify which couple has been provided and the order of coefficients.
//Exhaustive list of possible third argument :
//‘Enu’,’nuE,’Kmu’,’muK’, ‘KG’, ‘GK’, ‘lambdamu’, ‘mulambda’, ‘lambdaG’, ‘Glambda’.
py::array_t<double> M_iso(const double &C1, const double &C2, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_iso(C1,C2,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic stiffness tensor for a cubic material.
//Arguments are the stiffness coefficients C11, C12 and C44 ('Cii'), or E, nu and G ('EnuG')
py::array_t<double> L_cubic(const double &C11, const double &C12, const double &C44, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_cubic(C11,C12,C44,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic compliance tensor for an isotropic material.
//Arguments are the stiffness coefficients C11, C12 and C44 ('Cii'), or E, nu and G ('EnuG')
py::array_t<double> M_cubic(const double &C11, const double &C12, const double &C44, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_cubic(C11,C12,C44,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic stiffness tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
//C11,C12,C13,C22,C23,C13,C44,C55,C66 ('Cii')
//E1,E2,E3,nu12,nu13,nu23,G12,G13,G23 ('EnuG')
py::array_t<double> L_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::L_ortho(C11,C12,C13,C22,C23,C33,C44,C55,C66,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic compliance tensor for an orthotropic material.
//Arguments are the stiffness coefficients Cii or E and nu's
//C11,C12,C13,C22,C23,C13,C44,C55,C66 ('Cii')
//E1,E2,E3,nu12,nu13,nu23,G12,G13,G23 ('EnuG')
py::array_t<double> M_ortho(const double &C11, const double &C12, const double &C13, const double &C22, const double &C23, const double &C33, const double &C44, const double &C55, const double &C66, const std::string &conv_py, const bool &copy){
//    std::string conv = bp::extract<std::string>(conv_py);
    mat m = simcoon::M_ortho(C11,C12,C13,C22,C23,C33,C44,C55,C66,conv_py);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic stiffness tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
py::array_t<double> L_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis, const bool &copy){
    mat m = simcoon::L_isotrans(EL,ET,nuTL,nuTT,GLT,axis);
    return carma::mat_to_arr(m, copy);
}

//Returns the elastic compliance tensor for an isotropic transverse material.
//Arguments are longitudinal Young modulus EL, transverse young modulus, Poisson’s ratio for loading along the longitudinal axis nuTL, Poisson’s ratio for loading along the transverse axis nuTT, shear modulus GLT and the axis of symmetry.
py::array_t<double> M_isotrans(const double &EL, const double &ET, const double &nuTL, const double &nuTT, const double &GLT, const int &axis, const bool &copy){
    mat m = simcoon::M_isotrans(EL,ET,nuTL,nuTT,GLT,axis);
    return carma::mat_to_arr(m, copy);
}
//Provides the viscous tensor H an isotropic material.
//The two first arguments are a couple of viscous coefficients (the first is bulk, the second is shear).
py::array_t<double> H_iso(const double &etaB, const double &etaS, const bool &copy){
    mat m = simcoon::H_iso(etaB,etaS);
    return carma::mat_to_arr(m, copy);
}

} //namepsace simpy
