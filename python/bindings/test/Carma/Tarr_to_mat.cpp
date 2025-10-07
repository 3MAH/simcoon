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

///@file Troundtrip.cpp
///@brief Test for the carma library
///@version 1.0

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <carma>
#include <armadillo>

namespace py = pybind11;

namespace testcarma {

int test_arr_to_mat_double(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    auto arr_p = arr.unchecked<2>();

    // get buffer for raw pointer
    py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t ci = 0; ci < arr_S1; ci++) {
        for (size_t ri = 0; ri < arr_S0; ri++) {
            arr_sum += arr_p(ri, ci);
        }
    }

    // call function to be tested
    arma::Mat<double> M = carma::arr_to_mat<double>(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_mat_double */

int test_arr_to_mat_long(py::array_t<int64_t>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    auto arr_p = arr.unchecked<2>();

    // get buffer for raw pointer
    py::buffer_info info = arr.request();

    // compute sum of array
    int64_t arr_sum = 0;
    for (size_t ci = 0; ci < arr_S1; ci++) {
        for (size_t ri = 0; ri < arr_S0; ri++) {
            arr_sum += arr_p(ri, ci);
        }
    }

    // call function to be tested
    arma::Mat<int64_t> M = carma::arr_to_mat<int64_t>(arr, copy);

    // ---------------------------------------------------------------
    int64_t mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (arr_sum != mat_sum)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_mat_long */

int test_arr_to_mat_double_copy(const py::array_t<double>& arr) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    auto arr_p = arr.unchecked<2>();

    // get buffer for raw pointer
    const py::buffer_info pre_info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t ci = 0; ci < arr_S1; ci++) {
        for (size_t ri = 0; ri < arr_S0; ri++) {
            arr_sum += arr_p(ri, ci);
        }
    }

    // call function to be tested
    arma::Mat<double> M = carma::arr_to_mat<double>(arr);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (pre_info.ptr == M.memptr())
        return 5;
    return 0;
} /* arr_to_mat_double_copy */

int test_arr_to_mat_1d(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    auto arr_p = arr.unchecked<1>();

    // get buffer for raw pointer
    py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t i = 0; i < arr_N; i++)
        arr_sum += arr_p[i];

    // call function to be tested
    arma::Mat<double> M = carma::arr_to_mat<double>(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_mat_1d */

int test_arr_to_col(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    auto arr_p = arr.unchecked();

    // get buffer for raw pointer
    const py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t i = 0; i < arr_N; i++)
        arr_sum += arr_p[i];

    // call function to be tested
    arma::Col<double> M = carma::arr_to_col<double>(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_col */

int test_arr_to_row(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    auto arr_p = arr.unchecked();

    // get buffer for raw pointer
    const py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    if (arr.ndim() == 1) {
        for (size_t i = 0; i < arr_N; i++)
            arr_sum += arr_p[i];
    } else {
        for (size_t i = 0; i < arr_N; i++)
            arr_sum += arr_p(0, i);
    }

    // call function to be tested
    arma::Row<double> M = carma::arr_to_row<double>(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_col_double */

int test_arr_to_cube(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    size_t arr_S2 = arr.shape(2);
    auto arr_p = arr.unchecked<3>();

    // get buffer for raw pointer
    py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;

    for (size_t si = 0; si < arr_S2; si++) {
        for (size_t ci = 0; ci < arr_S1; ci++) {
            for (size_t ri = 0; ri < arr_S0; ri++) {
                arr_sum += arr_p(ri, ci, si);
            }
        }
    }

    // call function to be tested
    arma::Cube<double> M = carma::arr_to_cube<double>(arr, copy);

    // ---------------------------------------------------------------
    double cube_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (arr_S2 != M.n_slices)
        return 3;
    if (std::abs(arr_sum - cube_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_cube_double */

int test_to_arma_mat(py::array_t<double>& arr, bool copy) {
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    auto arr_p = arr.unchecked<2>();

    // get buffer for raw pointer
    const py::buffer_info pre_info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t ci = 0; ci < arr_S1; ci++) {
        for (size_t ri = 0; ri < arr_S0; ri++) {
            arr_sum += arr_p(ri, ci);
        }
    }

    // call function to be tested
    arma::Mat<double> M = carma::to_arma<arma::Mat<double>>::from(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (pre_info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_arma_mat */

int test_to_arma_cube(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    size_t arr_S0 = arr.shape(0);
    size_t arr_S1 = arr.shape(1);
    size_t arr_S2 = arr.shape(2);
    auto arr_p = arr.unchecked<3>();

    // get buffer for raw pointer
    py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;

    for (size_t si = 0; si < arr_S2; si++) {
        for (size_t ci = 0; ci < arr_S1; ci++) {
            for (size_t ri = 0; ri < arr_S0; ri++) {
                arr_sum += arr_p(ri, ci, si);
            }
        }
    }

    // call function to be tested
    arma::Cube<double> M = carma::to_arma<arma::Cube<double>>::from(arr, copy);

    // ---------------------------------------------------------------
    double cube_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (arr_S0 != M.n_rows)
        return 2;
    if (arr_S1 != M.n_cols)
        return 3;
    if (arr_S2 != M.n_slices)
        return 3;
    if (std::abs(arr_sum - cube_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_arma_cube */

int test_to_arma_col(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    auto arr_p = arr.unchecked();

    // get buffer for raw pointer
    const py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t i = 0; i < arr_N; i++)
        arr_sum += arr_p[i];

    // call function to be tested
    arma::Col<double> M = carma::to_arma<arma::Col<double>>::from(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_arma_col */

int test_to_arma_row(py::array_t<double>& arr, bool copy) {
    // attributes of the numpy array
    size_t arr_N = arr.size();
    auto arr_p = arr.unchecked();

    // get buffer for raw pointer
    const py::buffer_info info = arr.request();

    // compute sum of array
    double arr_sum = 0;
    for (size_t i = 0; i < arr_N; i++)
        arr_sum += arr_p[i];

    // call function to be tested
    arma::Row<double> M = carma::to_arma<arma::Row<double>>::from(arr, copy);

    // ---------------------------------------------------------------
    double mat_sum = arma::accu(M);

    // variable for test status
    if (arr_N != M.n_elem)
        return 1;
    if (std::abs(arr_sum - mat_sum) > 1e-8)
        return 4;
    if (info.ptr != M.memptr())
        return 5;
    return 0;
} /* arr_to_arma_row */

py::array_t<double> debug_arr_to_mat(py::array_t<double>& arr, int copy) {
    if (copy < 0) {
        return carma::mat_to_arr(carma::arr_to_mat<double>(std::move(arr)));
    }
    return carma::mat_to_arr(carma::arr_to_mat<double>(arr, copy));
}

}  // namespace testcarma
