#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>
#include <carma>
#include <armadillo>

#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/python_wrappers/Libraries/Maths/rotation.hpp>

using namespace std;
using namespace arma;
namespace py=pybind11;

namespace simpy {

py::array_t<double> rotate_vec_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &copy) {
    vec v_cpp = carma::arr_to_col(input);
    mat R_cpp = carma::arr_to_mat(R);
    vec rotated_vec = simcoon::rotate_vec(v_cpp, R_cpp);
    return carma::col_to_arr(rotated_vec, copy);
}

py::array_t<double> rotate_vec_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &copy) {
    vec v_cpp = carma::arr_to_col(input);
    vec rotated_vec = simcoon::rotate_vec(v_cpp, angle, axis);
    return carma::col_to_arr(rotated_vec, copy);
}

py::array_t<double> rotate_mat_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &copy) {
    mat m_cpp = carma::arr_to_mat(input);
    mat R_cpp = carma::arr_to_mat(R);
    mat rotated_mat = simcoon::rotate_mat(m_cpp, R_cpp);
    return carma::mat_to_arr(rotated_mat, copy);
}

py::array_t<double> rotate_mat_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &copy) {
    mat m_cpp = carma::arr_to_mat(input);
    mat rotated_mat = simcoon::rotate_mat(m_cpp, angle, axis);
    return carma::mat_to_arr(rotated_mat, copy);
}

//This function rotates a 6x6 stiffness matrix
py::array_t<double> rotate_stiffness_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat stiffness = carma::arr_to_mat(input);
    mat rotated_stiffness = simcoon::rotate_stiffness(stiffness, angle, axis, active);
    return carma::mat_to_arr(rotated_stiffness, copy);
}

//This function rotates a 6x6 stiffness matrix
py::array_t<double> rotate_stiffness_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat stiffness = carma::arr_to_mat(input);
    mat R_cpp = carma::arr_to_mat(R);
    mat rotated_stiffness = simcoon::rotate_stiffness(stiffness, R_cpp, active);
    return carma::mat_to_arr(rotated_stiffness, copy);
}

//This function rotates a 6x6 compliance matrix
py::array_t<double> rotate_compliance_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat compliance = carma::arr_to_mat(input);
    mat rotated_compliance = simcoon::rotate_compliance(compliance, angle, axis, active);
    return carma::mat_to_arr(rotated_compliance, copy);
}

//This function rotates a 6x6 compliance matrix
py::array_t<double> rotate_compliance_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat compliance = carma::arr_to_mat(input);
    mat R_cpp = carma::arr_to_mat(R);
    mat rotated_compliance = simcoon::rotate_compliance(compliance, R_cpp, active);
    return carma::mat_to_arr(rotated_compliance, copy);
}

//This function rotates a 6x6 strain concentration tensor
py::array_t<double> rotate_strain_concentration_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat strain_concentration = carma::arr_to_mat(input);
    mat rotated_strain_concentration = simcoon::rotate_strain_concentration(strain_concentration, angle, axis, active);
    return carma::mat_to_arr(rotated_strain_concentration, copy);
}

//This function rotates a 6x6 strain concentration tensor
py::array_t<double> rotate_strain_concentration_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat strain_concentration = carma::arr_to_mat(input);
    mat R_cpp = carma::arr_to_mat(R);
    mat rotated_strain_concentration = simcoon::rotate_strain_concentration(strain_concentration, R_cpp, active);
    return carma::mat_to_arr(rotated_strain_concentration, copy);
}

//This function rotates a 6x6 stress concentration tensor
py::array_t<double> rotate_stress_concentration_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    mat stress_concentration = carma::arr_to_mat(input);
    mat rotated_stress_concentration = simcoon::rotate_stress_concentration(stress_concentration, angle, axis, active);
    return carma::mat_to_arr(rotated_stress_concentration, copy);
}

//This function rotates a 6x6 stress concentration tensor
py::array_t<double> rotate_stress_concentration_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    mat stress_concentration = carma::arr_to_mat(input);
    mat R_cpp = carma::arr_to_mat(R);
    mat rotated_stress_concentration = simcoon::rotate_stress_concentration(stress_concentration, R_cpp, active);
    return carma::mat_to_arr(rotated_stress_concentration, copy);
}

//This function rotates stress vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_stress_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec stress = carma::arr_to_col(input);
        vec rotated_stress = simcoon::rotate_stress(stress, angle, axis, active);
        return carma::col_to_arr(rotated_stress, copy);
    }
    else if (input.ndim() == 2){
        mat stress_batch = carma::arr_to_mat_view(input);
        int nb_points = stress_batch.n_cols;
        mat rotated_stress(6, nb_points);
        for (int pt = 0; pt < nb_points; pt++) {
            vec stress = stress_batch.unsafe_col(pt);
            rotated_stress.col(pt) = simcoon::rotate_stress(stress, angle, axis, active);
        }
        return carma::mat_to_arr(rotated_stress, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates stress vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_stress_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec stress = carma::arr_to_col(input);
        mat R_cpp = carma::arr_to_mat(R);
        vec rotated_stress = simcoon::rotate_stress(stress, R_cpp, active);
        return carma::col_to_arr(rotated_stress, copy);
    }
    else if (input.ndim() == 2){
        mat stress_batch = carma::arr_to_mat_view(input);
        cube R_cpp = carma::arr_to_cube_view(R);
        int nb_points = stress_batch.n_cols;
        mat rotated_stress(6, nb_points);
        for (int pt = 0; pt < nb_points; pt++) {
            vec stress = stress_batch.unsafe_col(pt);
            rotated_stress.col(pt) = simcoon::rotate_stress(stress, R_cpp.slice(pt), active);
        }
        return carma::mat_to_arr(rotated_stress, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates strain vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_strain_angle(const py::array_t<double> &input, const double &angle, const int &axis, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec strain = carma::arr_to_col(input);
        vec rotated_strain = simcoon::rotate_strain(strain, angle, axis, active);
        return carma::col_to_arr(rotated_strain, copy);
    }
    else if (input.ndim() == 2){
        mat strain_batch = carma::arr_to_mat_view(input);
        int nb_points = strain_batch.n_cols;
        mat rotated_strain(6, nb_points);
        for (int pt = 0; pt < nb_points; pt++) {
            vec strain = strain_batch.unsafe_col(pt);
            rotated_strain.col(pt) = simcoon::rotate_strain(strain, angle, axis, active);
        }
        return carma::mat_to_arr(rotated_strain, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

//This function rotates strain vectors - Can be used with stack of arrays (vectorized)
py::array_t<double> rotate_strain_R(const py::array_t<double> &input, const py::array_t<double> &R, const bool &active, const bool &copy) {
    if (input.ndim()==1){
        vec strain = carma::arr_to_col(input);
        mat R_cpp = carma::arr_to_mat(R);
        vec rotated_strain = simcoon::rotate_strain(strain, R_cpp, active);
        return carma::col_to_arr(rotated_strain, copy);
    }
    else if (input.ndim() == 2){
        mat strain_batch = carma::arr_to_mat_view(input);
        cube R_cpp = carma::arr_to_cube_view(R);
        int nb_points = strain_batch.n_cols;
        mat rotated_strain(6, nb_points);
        for (int pt = 0; pt < nb_points; pt++) {
            vec strain = strain_batch.unsafe_col(pt);
            rotated_strain.col(pt) = simcoon::rotate_strain(strain, R_cpp.slice(pt), active);
        }
        return carma::mat_to_arr(rotated_strain, copy);
    }
    else{
        throw std::invalid_argument("input.ndim should be 1 or 2");
    }
}

} //namespace simpy
