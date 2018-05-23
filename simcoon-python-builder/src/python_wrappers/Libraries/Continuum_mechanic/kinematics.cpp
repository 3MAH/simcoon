
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/python_wrappers/Libraries/Continuum_mechanics/kinematics.hpp>

namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

//This function computes the gradient of displacement (Lagrangian) from the deformation gradient tensor
bn::ndarray G_UdX(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::G_UdX(F));
}

//This function computes the gradient of displacement (Eulerian) from the deformation gradient tensor
bn::ndarray G_Udx(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::G_Udx(F));
}

//This function computes the Right Cauchy-Green C
bn::ndarray R_Cauchy_Green(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::R_Cauchy_Green(F));
}

//This function computes the Left Cauchy-Green B
bn::ndarray L_Cauchy_Green(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::L_Cauchy_Green(F));
}

//This function computes the common Right (or Left) Cauchy-Green invariants
bn::ndarray Inv_X(const bn::ndarray &nd) {
    mat X = array2mat(nd);
    return vec2array(simcoon::Inv_X(X));
}

//This function computes the Cauchy deformation tensor c
bn::ndarray Cauchy(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::Cauchy(F));
}

//This function computes the Green-Lagrange finite strain tensor E
bn::ndarray Green_Lagrange(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::Green_Lagrange(F));
}

//This function computes the Euler-Almansi finite strain tensor A
bn::ndarray Euler_Almansi(const bn::ndarray &nd) {
    mat F = array2mat(nd);
    return mat2array(simcoon::Euler_Almansi(F));
}

//This function computes the velocity difference
bn::ndarray finite_L(const bn::ndarray &ndF, const bn::ndarray &ndDF, const double &DTime) {
    mat F = array2mat(ndF);
    mat DF = array2mat(ndDF);
    return mat2array(simcoon::finite_L(F,DF,DTime));
}

//This function computes the spin tensor W (correspond to Jaumann rate)
bn::ndarray finite_W(const bn::ndarray &ndF, const bn::ndarray &ndDF, const double &DTime) {
    mat F = array2mat(ndF);
    mat DF = array2mat(ndDF);
    return mat2array(simcoon::finite_L(F,DF,DTime));
}

//This function computes the spin tensor Omega (corrspond to Green-Naghdi rate)
// Note : here R is the is the rigid body rotation in the polar decomposition of the deformation gradient F
bn::ndarray finite_Omega(const bn::ndarray &ndR, const bn::ndarray &ndDR, const double &DTime) {
    mat R = array2mat(ndR);
    mat DR = array2mat(ndDR);
    return mat2array(simcoon::finite_Omega(R,DR,DTime));
}

//This function computes the deformation rate D
bn::ndarray finite_D(const bn::ndarray &ndF, const bn::ndarray &ndDF, const double &DTime) {
    mat F = array2mat(ndF);
    mat DF = array2mat(ndDF);
    return mat2array(simcoon::finite_D(F,DF,DTime));
}

//This function computes the increment of finite rotation
bn::ndarray finite_DQ(const bn::ndarray &ndOmega0, const bn::ndarray &ndOmega1, const double &DTime) {
    mat Omega0 = array2mat(ndOmega0);
    mat Omega1 = array2mat(ndOmega1);
    return mat2array(simcoon::finite_DQ(Omega0,Omega1,DTime));
}
    
} //namepsace simpy
