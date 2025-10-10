#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <armadillo>
#include <simcoon/python_wrappers/conversion_helpers.hpp>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

py::array_t<double> rotate_vec_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &copy) {
    vec v_cpp = simpy::arr_to_col(input);
    mat R_cpp = simpy::arr_to_mat(R);
    vec rotated_v = simcoon::rotate_vec(v_cpp,R_cpp);   
    return simpy::col_to_arr(rotated_v, copy);
}

py::array_t<double> rotate_vec_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &copy) {
    vec v_cpp = simpy::arr_to_col(input);
    vec rotated_v = simcoon::rotate_vec(v_cpp,angle,axis);     
    return simpy::col_to_arr(rotated_v, copy);
}

py::array_t<double> rotate_mat_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &copy) {
    mat m_cpp = simpy::arr_to_mat(input);
    mat R_cpp = simpy::arr_to_mat(R);
    mat rotated_m = simcoon::rotate_mat(m_cpp,R_cpp);
    return simpy::mat_to_arr(rotated_m, copy);
}

py::array_t<double> rotate_mat_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &copy) {    
    mat m_cpp = simpy::arr_to_mat(input);
    mat rotated_m = simcoon::rotate_mat(m_cpp,angle,axis);    
    return simpy::mat_to_arr(rotated_m, copy);
}

py::array_t<double> fillR_angle(const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat R = simcoon::fillR(angle,axis,active);
    return simpy::mat_to_arr(R, copy);
}
    
//This function returns the 3*3 rotation matrix
py::array_t<double> fillR_euler(const double &psi, const double &theta, const double &phi, const bool &active, const string &conv, const bool &copy) {
    mat R = simcoon::fillR(psi,theta,phi,active,conv);
    return simpy::mat_to_arr(R, copy);
}

//This function returns the 6*6 rotation matrix of a vector of type 'stress'
py::array_t<double> fillQS_angle(const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat QS = simcoon::fillQS(angle,axis,active);
    return simpy::mat_to_arr(QS, copy);
}

//This function returns the 6*6 rotation matrix of a vector of type 'stress'
py::array_t<double> fillQS_R(const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat R_cpp = simpy::arr_to_mat(R);
    mat QS = simcoon::fillQS(R_cpp,active);
    return simpy::mat_to_arr(QS, copy);
}
    
//This function returns the 6*6 rotation matrix of a vector of type 'strain'
py::array_t<double> fillQE_angle(const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat QE = simcoon::fillQE(angle,axis,active);
    return simpy::mat_to_arr(QE, copy);
}

//This function returns the 6*6 rotation matrix of a vector of type 'strain'
py::array_t<double> fillQE_R(const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat R_cpp = simpy::arr_to_mat(R);
    mat QE = simcoon::fillQE(R_cpp,active);
    return simpy::mat_to_arr(QE, copy);
}

//This function rotates a 6*6 stiffness matrix (L)
py::array_t<double> rotateL_angle(const py::array_t<double> &input, const double &angle, const int & axis, const bool &active, const bool &copy) {
    mat L_cpp = simpy::arr_to_mat(input);
    mat rotated_L = simcoon::rotateL(L_cpp,angle,axis,active); 
    return simpy::mat_to_arr(rotated_L, copy);
}

//This function rotates a 6*6 stiffness matrix (L)
py::array_t<double> rotateL_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat L_cpp = simpy::arr_to_mat(input);
    mat R_cpp = simpy::arr_to_mat(R);
    mat rotated_L = simcoon::rotateL(L_cpp,R_cpp,active);
    return simpy::mat_to_arr(rotated_L, copy);
}

//This function rotates a 6*6 compliance matrix (M)
py::array_t<double> rotateM_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat M_cpp = simpy::arr_to_mat(input);
    mat rotated_M = simcoon::rotateM(M_cpp,angle,axis,active);
    return simpy::mat_to_arr(rotated_M, copy);
}

//This function rotates a 6*6 compliance matrix (M)
py::array_t<double> rotateM_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat M_cpp = simpy::arr_to_mat(input);
    mat R_cpp = simpy::arr_to_mat(R);
    mat rotated_M = simcoon::rotateM(M_cpp,R_cpp,active);
    return simpy::mat_to_arr(rotated_M);
}

//This function rotates a 6*6 strain concentration (A)
py::array_t<double> rotateA_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat A_cpp = simpy::arr_to_mat(input);
    mat rotated_A = simcoon::rotateA(A_cpp,angle,axis,active);
    return simpy::mat_to_arr(rotated_A, copy);
}

//This function rotates a 6*6 strain concentration (A)
py::array_t<double> rotateA_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat A_cpp = simpy::arr_to_mat(input);
    mat R_cpp = simpy::arr_to_mat(R);
    mat rotated_A = simcoon::rotateA(A_cpp,R_cpp,active);  
    return simpy::mat_to_arr(rotated_A, copy);
}

//This function rotates a 6*6 stress concentration (B)
py::array_t<double> rotateB_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat B_cpp = simpy::arr_to_mat(input);
    mat rotated_B = simcoon::rotateB(B_cpp,angle,axis,active);  
    return simpy::mat_to_arr(rotated_B, copy);
}

//This function rotates a 6*6 stress concentration (B)
py::array_t<double> rotateB_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat B_cpp = simpy::arr_to_mat(input);
    mat R_cpp = simpy::arr_to_mat(R);
    mat rotated_B = simcoon::rotateB(B_cpp,R_cpp,active);
    return simpy::mat_to_arr(rotated_B, copy);
}

//This function rotates stress vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_stress_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec v = simpy::arr_to_col(input);
        vec rotated_stress = simcoon::rotate_stress(v,angle,axis,active);
        return simpy::col_to_arr(rotated_stress, copy);
    }
    else if (input.ndim() == 2){
        mat m = simpy::arr_to_mat_view(input);
        int nb_points = m.n_cols;
        mat rotated_stress(6,nb_points); 
        for (int pt = 0; pt < nb_points; pt++) {
            vec v = m.unsafe_col(pt);
            rotated_stress.col(pt) = simcoon::rotate_stress(v,angle,axis,active);
        }
        return simpy::mat_to_arr(rotated_stress, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}
                 
//This function rotates stress vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_stress_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec v = simpy::arr_to_col(input);
        mat R_cpp = simpy::arr_to_mat(R);
        vec rotated_stress = simcoon::rotate_stress(v,R_cpp,active);
        return simpy::col_to_arr(rotated_stress, copy);
    }
    else if (input.ndim() == 2){
        mat m = simpy::arr_to_mat_view(input);
        cube R_cpp = simpy::arr_to_cube_view(R);
        int nb_points = m.n_cols;
        mat rotated_stress(6,nb_points); 
        for (int pt = 0; pt < nb_points; pt++) {
            vec v = m.unsafe_col(pt);
            rotated_stress.col(pt) = simcoon::rotate_stress(v,R_cpp.slice(pt),active);
        }
        return simpy::mat_to_arr(rotated_stress, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates strain vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_strain_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec v = simpy::arr_to_col(input);
        vec rotated_strain = simcoon::rotate_strain(v,angle,axis,active);
        return simpy::col_to_arr(rotated_strain, copy);
    }
    else if (input.ndim() == 2){
        mat m = simpy::arr_to_mat_view(input);
        int nb_points = m.n_cols;
        mat rotated_strain(6,nb_points); 
        for (int pt = 0; pt < nb_points; pt++) {
            vec v = m.unsafe_col(pt);
            rotated_strain.col(pt) = simcoon::rotate_strain(v,angle,axis,active);
        }
        return simpy::mat_to_arr(rotated_strain, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates strain vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_strain_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec v = simpy::arr_to_col(input);
        mat R_cpp = simpy::arr_to_mat(R);
        vec rotated_strain = simcoon::rotate_strain(v,R_cpp,active);
        return simpy::col_to_arr(rotated_strain, copy);
    }
    else if (input.ndim() == 2){
        mat m = simpy::arr_to_mat_view(input);
        cube R_cpp = simpy::arr_to_cube_view(R);
        int nb_points = m.n_cols;
        mat rotated_strain(6,nb_points); 
        for (int pt = 0; pt < nb_points; pt++) {
            vec v = m.unsafe_col(pt);
            rotated_strain.col(pt) = simcoon::rotate_strain(v,R_cpp.slice(pt),active);
        }
        return simpy::mat_to_arr(rotated_strain, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates strain vectors from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_strain(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    vec E_cpp = simpy::arr_to_col(input);
    vec rotated_E = simcoon::rotate_l2g_strain(E_cpp,psi,theta,phi);
    return simpy::col_to_arr(rotated_E, copy);
}

//This function rotates strain vectors from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_strain(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    vec E_cpp = simpy::arr_to_col(input);
    vec rotated_E = simcoon::rotate_g2l_strain(E_cpp,psi,theta,phi);
    return simpy::col_to_arr(rotated_E, copy);
}

//This function rotates stress vectors from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_stress(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
 vec S_cpp = simpy::arr_to_col(input);
 vec rotated_S = simcoon::rotate_l2g_strain(S_cpp,psi,theta,phi);
 return simpy::col_to_arr(rotated_S, copy);
}

//This function rotates stress vectors from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_stress(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
 vec S_cpp = simpy::arr_to_col(input);
 vec rotated_S = simcoon::rotate_g2l_stress(S_cpp,psi,theta,phi);
 return simpy::col_to_arr(rotated_S, copy);
}

//This function rotates stiffness matrices (L) from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_L(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat L_cpp = simpy::arr_to_mat(input);
    mat rotated_L = simcoon::rotate_l2g_L(L_cpp,psi,theta,phi);
    return simpy::mat_to_arr(rotated_L, copy);
}

//This function rotates stiffness matrices (L) from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_L(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat L_cpp = simpy::arr_to_mat(input);
    mat rotated_L = simcoon::rotate_g2l_L(L_cpp,psi,theta,phi);
    return simpy::mat_to_arr(rotated_L, copy);
}

//This function rotates compliance matrices (M) from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_M(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat M_cpp = simpy::arr_to_mat(input);
    mat rotated_M = simcoon::rotate_l2g_M(M_cpp,psi,theta,phi);
    return simpy::mat_to_arr(rotated_M, copy);
}

//This function rotates compliance matrices (M) from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_M(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat M_cpp = simpy::arr_to_mat(input);
    mat rotated_M = simcoon::rotate_g2l_M(M_cpp,psi,theta,phi);   
    return simpy::mat_to_arr(rotated_M, copy);
}

//This function rotates strain concentration matrices (A) from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_A(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat A_cpp = simpy::arr_to_mat(input);
    mat rotated_A = simcoon::rotate_l2g_A(A_cpp,psi,theta,phi);   
    return simpy::mat_to_arr(rotated_A, copy);
}

//This function rotates strain concentration matrices (A) from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_A(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat A_cpp = simpy::arr_to_mat(input);
    mat rotated_A = simcoon::rotate_g2l_A(A_cpp,psi,theta,phi);    
    return simpy::mat_to_arr(rotated_A, copy);
}

//This function rotates stress concentration matrices (B) from a local to global set of coordinates (using Euler angles)
py::array_t<double> rotate_l2g_B(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat B_cpp = simpy::arr_to_mat(input);
    mat rotated_B = simcoon::rotate_l2g_B(B_cpp,psi,theta,phi); 
    return simpy::mat_to_arr(rotated_B, copy);
}

//This function rotates stress concentration matrices (B) from a global to local set of coordinates (using Euler angles)
py::array_t<double> rotate_g2l_B(const py::array_t<double> &input, const double &psi, const double &theta, const double &phi, const bool &copy) {
    mat B_cpp = simpy::arr_to_mat(input);
    mat rotated_B = simcoon::rotate_g2l_B(B_cpp,psi,theta,phi); 
    return simpy::mat_to_arr(rotated_B, copy);
}
    
} //namepsace simpy
