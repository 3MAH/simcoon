#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <armadillo>
#include <algorithm>

namespace py = pybind11;

namespace simpy {

/**
 * Conversion helpers to replace carma functionality.
 * These functions provide explicit copying between numpy arrays and armadillo containers,
 * avoiding the allocator mismatch issues that occur with carma's zero-copy approach on Windows.
 */

// ============================================================================
// numpy array → armadillo vector
// ============================================================================

/**
 * Convert 1D numpy array to armadillo column vector.
 * Replaces: carma::arr_to_col(arr)
 */
inline arma::vec arr_to_col(const py::array_t<double>& arr) {
    auto r = arr.request();
    const double* ptr = static_cast<const double*>(r.ptr);
    arma::vec v(r.shape[0]);
    std::copy(ptr, ptr + r.shape[0], v.begin());
    return v;
}

/**
 * Convert 1D numpy array to armadillo column vector (view version).
 * Note: This creates a COPY, not a view, to avoid allocator mismatch.
 * Replaces: carma::arr_to_col_view(arr)
 */
inline arma::vec arr_to_col_view(const py::array_t<double>& arr) {
    return arr_to_col(arr);  // Same as copy for safety
}

// ============================================================================
// numpy array → armadillo matrix
// ============================================================================

/**
 * Convert 2D numpy array to armadillo matrix.
 * Replaces: carma::arr_to_mat(arr)
 * Note: Handles row-major (numpy) to column-major (armadillo) conversion
 */
inline arma::mat arr_to_mat(const py::array_t<double>& arr) {
    auto r = arr.request();
    const double* ptr = static_cast<const double*>(r.ptr);
    arma::mat m(r.shape[0], r.shape[1]);

    // Copy element by element (numpy is row-major, armadillo is column-major)
    for (py::ssize_t i = 0; i < r.shape[0]; i++) {
        for (py::ssize_t j = 0; j < r.shape[1]; j++) {
            m(i, j) = ptr[i * r.shape[1] + j];
        }
    }
    return m;
}

/**
 * Convert 2D numpy array to armadillo matrix (view version).
 * Note: This creates a COPY, not a view, to avoid allocator mismatch.
 * Replaces: carma::arr_to_mat_view(arr)
 */
inline arma::mat arr_to_mat_view(const py::array_t<double>& arr) {
    return arr_to_mat(arr);  // Same as copy for safety
}

// ============================================================================
// numpy array → armadillo cube
// ============================================================================

/**
 * Convert 3D numpy array to armadillo cube.
 * Replaces: carma::arr_to_cube(arr)
 */
inline arma::cube arr_to_cube(const py::array_t<double>& arr) {
    auto r = arr.request();
    const double* ptr = static_cast<const double*>(r.ptr);
    arma::cube c(r.shape[0], r.shape[1], r.shape[2]);

    // Copy element by element
    for (py::ssize_t i = 0; i < r.shape[0]; i++) {
        for (py::ssize_t j = 0; j < r.shape[1]; j++) {
            for (py::ssize_t k = 0; k < r.shape[2]; k++) {
                c(i, j, k) = ptr[i * r.shape[1] * r.shape[2] + j * r.shape[2] + k];
            }
        }
    }
    return c;
}

/**
 * Convert 3D numpy array to armadillo cube (view version).
 * Note: This creates a COPY, not a view, to avoid allocator mismatch.
 * Replaces: carma::arr_to_cube_view(arr)
 */
inline arma::cube arr_to_cube_view(const py::array_t<double>& arr) {
    return arr_to_cube(arr);  // Same as copy for safety
}

// ============================================================================
// armadillo vector → numpy array
// ============================================================================

/**
 * Convert armadillo column vector to 1D numpy array.
 * Replaces: carma::col_to_arr(v, copy)
 */
inline py::array_t<double> col_to_arr(const arma::vec& v, bool copy = true) {
    py::array_t<double> result(v.n_elem);
    auto r = result.request();
    double* ptr = static_cast<double*>(r.ptr);
    std::copy(v.begin(), v.end(), ptr);
    return result;
}

// ============================================================================
// armadillo matrix → numpy array
// ============================================================================

/**
 * Convert armadillo matrix to 2D numpy array.
 * Replaces: carma::mat_to_arr(m, copy)
 * Note: Armadillo uses column-major, numpy default is row-major, so we copy element-by-element
 */
inline py::array_t<double> mat_to_arr(const arma::mat& m, bool copy = true) {
    py::array_t<double> result({m.n_rows, m.n_cols});
    auto r = result.request();
    double* ptr = static_cast<double*>(r.ptr);

    // Copy element by element (armadillo is column-major, numpy is row-major)
    for (size_t i = 0; i < m.n_rows; i++) {
        for (size_t j = 0; j < m.n_cols; j++) {
            ptr[i * m.n_cols + j] = m(i, j);
        }
    }
    return result;
}

// ============================================================================
// armadillo cube → numpy array
// ============================================================================

/**
 * Convert armadillo cube to 3D numpy array.
 * Replaces: carma::cube_to_arr(c, copy)
 */
inline py::array_t<double> cube_to_arr(const arma::cube& c, bool copy = true) {
    py::array_t<double> result({c.n_rows, c.n_cols, c.n_slices});
    auto r = result.request();
    double* ptr = static_cast<double*>(r.ptr);

    // Copy element by element
    for (size_t i = 0; i < c.n_rows; i++) {
        for (size_t j = 0; j < c.n_cols; j++) {
            for (size_t k = 0; k < c.n_slices; k++) {
                ptr[i * c.n_cols * c.n_slices + j * c.n_slices + k] = c(i, j, k);
            }
        }
    }
    return result;
}

} // namespace simpy
