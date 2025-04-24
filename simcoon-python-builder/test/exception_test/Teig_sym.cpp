/* This file is part of simcoon.
 
 simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 simcoon is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with simcoon.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file Teig_sym.cpp
///@brief Test for the exceptions raised with simcoon
///@version 1.0

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <carma>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>

using namespace std;
namespace py = pybind11;
namespace testexception {

py::array_t<double> test_eig_sym_val_affect(const py::array_t<double> &matrix_fail_eigen, const bool &copy) {

    //Try eigen decomposition : fail since decomposition is not unique
    arma::mat arma_fail_eigen = carma::arr_to_mat(matrix_fail_eigen);

    arma::vec eigval;

    try {
        eigval = arma::eig_sym(arma_fail_eigen);
    } catch (const std::runtime_error &e) {
        cerr << "Error in eig_sym: " << e.what() << endl;
        throw simcoon::exception_eig_sym("Error in eig_sym function inside test_eig_sym_val_affect.");
    }
}    

py::array_t<double> test_eig_sym_val_modify(const py::array_t<double> &matrix_fail_eigen, const bool &copy) {

    //Try eigen decomposition : fail since decomposition is not unique
    arma::mat arma_fail_eigen = carma::arr_to_mat(matrix_fail_eigen);
    arma::vec eigval;

    bool success = arma::eig_sym(eigval, arma_fail_eigen);

    if (success) {
        eigval.print("Eigenvalues:");
        return carma::col_to_arr(eigval, copy);       
    } else {
        cout << "Eigen decomposition failed." << endl;
        throw simcoon::exception_eig_sym("Error in eig_sym function inside test_eig_sym_val_modify.");
    }
}

py::tuple test_eig_sym_val_vec(const py::array_t<double> &matrix_fail_eigen, const bool &copy) {

    //Try eigen decomposition : fail since decomposition is not unique
    arma::mat arma_fail_eigen = carma::arr_to_mat(matrix_fail_eigen);
    arma::vec eigval;
    arma::mat eigvec;

    bool success = arma::eig_sym(eigval, eigvec, arma_fail_eigen);

    if (success) {
        eigval.print("Eigenvalues:");
        eigvec.print("Eigenvectors:");
        return py::make_tuple(carma::col_to_arr(eigval, copy), carma::mat_to_arr(eigvec, copy));        
    } else {
        cout << "Eigen decomposition failed." << endl;
        throw simcoon::exception_eig_sym("Error in eig_sym function inside test_eig_sym_val_vec.");
    }
}





}  // namespace testexception
