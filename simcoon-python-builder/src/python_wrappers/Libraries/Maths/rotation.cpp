
#include <armadillo>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <simcoon/arma2numpy/numpy_arma.hpp>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>

namespace bn = boost::python::numpy;
using namespace std;
using namespace arma;
using namespace arma2numpy;

namespace simpy {

bn::ndarray rotate_vec_R(const bn::ndarray &ndv, const bn::ndarray &ndR) {
    vec v = array2vec(ndv);
    mat R = array2mat(ndR);
    return vec2array(simcoon::rotate_vec(v,R));
}

bn::ndarray rotate_vec_angle(const bn::ndarray &ndv, const double &alpha, const int &axis) {
    vec v = array2vec(ndv);
    return vec2array(simcoon::rotate_vec(v,alpha,axis));
}

bn::ndarray rotate_mat_R(const bn::ndarray &ndm, const bn::ndarray &ndR) {
    mat m = array2mat(ndm);
    mat R = array2mat(ndR);
    return mat2array(simcoon::rotate_mat(m,R));
}

bn::ndarray rotate_mat_angle(const bn::ndarray &ndm, const double &alpha, const int &axis) {
    
    mat m = array2mat(ndm);
    return mat2array(simcoon::rotate_mat(m,alpha,axis));
}

bn::ndarray fillR_angle(const double &alpha, const int &axis, const bool &active) {
    
    return mat2array(simcoon::fillR(alpha,axis,active));
}
    
//This function returns the 3*3 rotation matrix
bn::ndarray fillR_euler(const double &psi, const double &theta, const double &phi, const bool &active, const string &conv) {
    return mat2array(simcoon::fillR(psi,theta,phi,active,conv));
}

//This function returns the 6*6 rotation matrix of a vector of type 'stress'
bn::ndarray fillQS_angle(const double &alpha, const int &axis, const bool &active) {
    return mat2array(simcoon::fillQS(alpha,axis,active));
}

//This function returns the 6*6 rotation matrix of a vector of type 'stress'
bn::ndarray fillQS_R(const bn::ndarray &ndR, const bool &active) {
    mat R = array2mat(ndR);
    return mat2array(simcoon::fillQS(R,active));
}
    
//This function returns the 6*6 rotation matrix of a vector of type 'strain'
bn::ndarray fillQE_angle(const double &alpha, const int &axis, const bool &active) {
    return mat2array(simcoon::fillQE(alpha,axis,active));
}

//This function returns the 6*6 rotation matrix of a vector of type 'strain'
bn::ndarray fillQE_R(const bn::ndarray &ndR, const bool &active) {
    mat R = array2mat(ndR);
    return mat2array(simcoon::fillQE(R,active));
}

//This function rotates a 6*6 stiffness matrix (L)
bn::ndarray rotateL_angle(const bn::ndarray &ndL, const double &alpha, const int & axis, const bool &active) {
    mat L = array2mat(ndL);
    return mat2array(simcoon::rotateL(L,alpha,axis,active));
}

//This function rotates a 6*6 stiffness matrix (L)
bn::ndarray rotateL_R(const bn::ndarray &ndL, const bn::ndarray &ndR, const bool &active) {
    mat L = array2mat(ndL);
    mat R = array2mat(ndR);
    return mat2array(simcoon::rotateL(L,R,active));
}

//This function rotates a 6*6 compliance matrix (M)
bn::ndarray rotateM_angle(const bn::ndarray &ndM, const double &alpha, const int &axis, const bool &active) {
    mat M = array2mat(ndM);
    return mat2array(simcoon::rotateM(M,alpha,axis,active));
}

//This function rotates a 6*6 compliance matrix (M)
bn::ndarray rotateM_R(const bn::ndarray &ndM, const bn::ndarray &ndR, const bool &active) {
    mat M = array2mat(ndM);
    mat R = array2mat(ndR);
    return mat2array(simcoon::rotateM(M,R,active));
}

//This function rotates a 6*6 strain concentration (A)
bn::ndarray rotateA_angle(const bn::ndarray &ndA, const double &alpha, const int &axis, const bool &active) {
    mat A = array2mat(ndA);
    return mat2array(simcoon::rotateA(A,alpha,axis,active));
}

//This function rotates a 6*6 strain concentration (A)
bn::ndarray rotateA_R(const bn::ndarray &ndA, const bn::ndarray &ndR, const bool &active) {
    mat A = array2mat(ndA);
    mat R = array2mat(ndR);
    return mat2array(simcoon::rotateA(A,R,active));
}

//This function rotates a 6*6 stress concentration (B)
bn::ndarray rotateB_angle(const bn::ndarray &ndB, const double &alpha, const int &axis, const bool &active) {
    mat B = array2mat(ndB);
    return mat2array(simcoon::rotateB(B,alpha,axis,active));
}

//This function rotates a 6*6 stress concentration (B)
bn::ndarray rotateB_R(const bn::ndarray &ndB, const bn::ndarray &ndR, const bool &active) {
    mat B = array2mat(ndB);
    mat R = array2mat(ndR);
    return mat2array(simcoon::rotateB(B,R,active));
}

//This function rotates stress vectors
bn::ndarray rotate_stress_angle(const bn::ndarray &nd, const double &alpha, const int &axis, const bool &active) {
    vec v = array2vec(nd);
    return vec2array(simcoon::rotate_stress(v,alpha,axis,active));
}
                 
//This function rotates stress vectors
bn::ndarray rotate_stress_R(const bn::ndarray &nd, const bn::ndarray &ndR, const bool &active) {
    vec v = array2vec(nd);
    mat R = array2mat(ndR);
    return vec2array(simcoon::rotate_stress(v,R,active));
}

//This function rotates strain vectors
bn::ndarray rotate_strain_angle(const bn::ndarray &nd, const double &alpha, const int &axis, const bool &active) {
    vec v = array2vec(nd);
    return vec2array(simcoon::rotate_strain(v,alpha,axis,active));
}

//This function rotates stress vectors
bn::ndarray rotate_strain_R(const bn::ndarray &nd, const bn::ndarray &ndR, const bool &active) {
    vec v = array2vec(nd);
    mat R = array2mat(ndR);
    return vec2array(simcoon::rotate_strain(v,R,active));
}

//This function rotates strain vectors from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_strain(const bn::ndarray &ndE, const double &psi, const double &theta, const double &phi) {
    vec E = array2vec(ndE);
    return vec2array(simcoon::rotate_l2g_strain(E,psi,theta,phi));
}

//This function rotates strain vectors from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_strain(const bn::ndarray &ndE, const double &psi, const double &theta, const double &phi) {
    vec E = array2vec(ndE);
    return vec2array(simcoon::rotate_g2l_strain(E,psi,theta,phi));
}

//This function rotates stress vectors from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_stress(const bn::ndarray &ndS, const double &psi, const double &theta, const double &phi) {
 vec S = array2vec(ndS);
 return vec2array(simcoon::rotate_l2g_strain(S,psi,theta,phi));
}

//This function rotates stress vectors from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_stress(const bn::ndarray &ndS, const double &psi, const double &theta, const double &phi) {
 vec S = array2vec(ndS);
 return vec2array(simcoon::rotate_g2l_stress(S,psi,theta,phi));
}

//This function rotates stiffness matrices (L) from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_L(const bn::ndarray &ndL, const double &psi, const double &theta, const double &phi) {
    mat L = array2mat(ndL);
    return mat2array(simcoon::rotate_l2g_L(L,psi,theta,phi));
}

//This function rotates stiffness matrices (L) from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_L(const bn::ndarray &ndL, const double &psi, const double &theta, const double &phi) {
    mat L = array2mat(ndL);
    return mat2array(simcoon::rotate_g2l_L(L,psi,theta,phi));
}

//This function rotates compliance matrices (M) from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_M(const bn::ndarray &ndM, const double &psi, const double &theta, const double &phi) {
    mat M = array2mat(ndM);
    return mat2array(simcoon::rotate_l2g_M(M,psi,theta,phi));
}

//This function rotates compliance matrices (M) from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_M(const bn::ndarray &ndM, const double &psi, const double &theta, const double &phi) {
    mat M = array2mat(ndM);
    return mat2array(simcoon::rotate_g2l_M(M,psi,theta,phi));
}

//This function rotates strain concentration matrices (A) from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_A(const bn::ndarray &ndA, const double &psi, const double &theta, const double &phi) {
    mat A = array2mat(ndA);
    return mat2array(simcoon::rotate_l2g_A(A,psi,theta,phi));
}

//This function rotates strain concentration matrices (A) from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_A(const bn::ndarray &ndA, const double &psi, const double &theta, const double &phi) {
    mat A = array2mat(ndA);
    return mat2array(simcoon::rotate_g2l_A(A,psi,theta,phi));
}

//This function rotates stress concentration matrices (B) from a local to global set of coordinates (using Euler angles)
bn::ndarray rotate_l2g_B(const bn::ndarray &ndB, const double &psi, const double &theta, const double &phi) {
 mat B = array2mat(ndB);
 return mat2array(simcoon::rotate_l2g_B(B,psi,theta,phi));
}

//This function rotates stress concentration matrices (B) from a global to local set of coordinates (using Euler angles)
bn::ndarray rotate_g2l_B(const bn::ndarray &ndB, const double &psi, const double &theta, const double &phi) {
 mat B = array2mat(ndB);
 return mat2array(simcoon::rotate_g2l_B(B,psi,theta,phi));
}
    
} //namepsace simpy
