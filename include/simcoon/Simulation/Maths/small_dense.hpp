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

///@file small_dense.hpp
///@brief Fixed-order 6x6 products for the per-point kernels
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup maths
 *  @{
 */

/**
 * @brief \f$ \mathbf{y} = \mathbf{A}\,\mathbf{x} \f$ for a column-major 6x6 @p A and a 6-vector @p x.
 *
 * Armadillo hands every product above 4x4 to BLAS; at this size the call costs more than the
 * 36 multiply-adds, and the summation order depends on the BLAS backend. This loop has a fixed
 * order (column by column), so the result is the same on every platform.
 * @p y may alias @p x.
 */
inline void mul6(const double *A, const double *x, double *y) {
    double r[6];
    for (int i = 0; i < 6; ++i) r[i] = A[i]*x[0];
    for (int k = 1; k < 6; ++k)
        for (int i = 0; i < 6; ++i) r[i] += A[i + 6*k]*x[k];
    for (int i = 0; i < 6; ++i) y[i] = r[i];
}

/// \f$ \mathbf{y} = \mathbf{A}\,\mathbf{x} \f$ for a 6x6 @p A and a 6-vector @p x (see the pointer form).
inline void mul6(const arma::mat &A, const arma::vec &x, arma::vec &y) {
    if (&y != &x) y.set_size(6);
    mul6(A.memptr(), x.memptr(), y.memptr());
}

/// \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$ for 6x6 operands, fixed summation order (see mul6).
inline arma::mat::fixed<6,6> mul66(const arma::mat &A, const arma::mat &B) {
    const double *a = A.memptr();
    const double *b = B.memptr();
    arma::mat::fixed<6,6> C;
    double *c = C.memptr();
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += a[i + 6*k]*b[k + 6*j];
            c[i + 6*j] = s;
        }
    return C;
}

/// \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B}^T \f$ for 6x6 operands, without forming the transpose.
inline arma::mat::fixed<6,6> mul66_ABt(const arma::mat &A, const arma::mat &B) {
    const double *a = A.memptr();
    const double *b = B.memptr();
    arma::mat::fixed<6,6> C;
    double *c = C.memptr();
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += a[i + 6*k]*b[j + 6*k];
            c[i + 6*j] = s;
        }
    return C;
}

/** @} */ // end of maths group

} //namespace simcoon
