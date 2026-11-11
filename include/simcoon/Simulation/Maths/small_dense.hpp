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

/// \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B} \f$ for 3x3 operands, each entry the three products summed in order.
inline arma::mat::fixed<3,3> mul3(const arma::mat::fixed<3,3> &A, const arma::mat::fixed<3,3> &B) {
    arma::mat::fixed<3,3> C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(0,j), x1 = B(1,j), x2 = B(2,j);
        for (int i = 0; i < 3; ++i) C(i,j) = A(i,0)*x0 + A(i,1)*x1 + A(i,2)*x2;
    }
    return C;
}

/// \f$ \mathbf{C} = \mathbf{A}^T\mathbf{B} \f$ for 3x3 operands, without forming the transpose.
inline arma::mat::fixed<3,3> mul3_AtB(const arma::mat::fixed<3,3> &A, const arma::mat::fixed<3,3> &B) {
    arma::mat::fixed<3,3> C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(0,j), x1 = B(1,j), x2 = B(2,j);
        for (int i = 0; i < 3; ++i) C(i,j) = A(0,i)*x0 + A(1,i)*x1 + A(2,i)*x2;
    }
    return C;
}

/// \f$ \mathbf{C} = \mathbf{A}\,\mathbf{B}^T \f$ for 3x3 operands, without forming the transpose.
inline arma::mat::fixed<3,3> mul3_ABt(const arma::mat::fixed<3,3> &A, const arma::mat::fixed<3,3> &B) {
    arma::mat::fixed<3,3> C;
    for (int j = 0; j < 3; ++j) {
        const double x0 = B(j,0), x1 = B(j,1), x2 = B(j,2);
        for (int i = 0; i < 3; ++i) C(i,j) = A(i,0)*x0 + A(i,1)*x1 + A(i,2)*x2;
    }
    return C;
}

/**
 * @brief Contraction \f$ \mathbf{A}:\mathbf{B} \f$ of two symmetric tensors stored as engineering Voigt
 * 6-vectors: \f$ \sum_{k<3} a_k b_k + w \sum_{k \ge 3} a_k b_k \f$. The shear weight @p w is 2 for two
 * stresses, 1/2 for two strains and 1 for a stress with a strain (see voigt_shear_weight in tensor.hpp).
 */
inline double contract_voigt(const double *a, const double *b, const double w) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + w*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]);
}

/// Deviator of a Voigt 6-vector, \f$ d = v - \tfrac{1}{3}\mathrm{tr}(v)\,\mathbf{I} \f$.
inline void dev_voigt(const double *v, double *d) {
    const double m = (v[0] + v[1] + v[2])/3.;
    d[0] = v[0] - m; d[1] = v[1] - m; d[2] = v[2] - m;
    d[3] = v[3]; d[4] = v[4]; d[5] = v[5];
}

/**
 * @brief Determinant of the symmetric 3x3 tensor stored as a Voigt 6-vector @p v whose shear entries map to
 * the tensor with the factor @p shear (voigt_shear_in: 1 for a stress, 1/2 for a strain); closed form, no
 * 3x3 built and no LAPACK dispatch.
 */
inline double det_sym_voigt(const double *v, const double shear) {
    const double s3 = shear*v[3], s4 = shear*v[4], s5 = shear*v[5];
    return v[0]*(v[1]*v[2] - s5*s5) - s3*(s3*v[2] - s4*s5) + s4*(s3*s5 - v[1]*s4);
}

/// \f$ \mathbf{C} = \mathbf{a}\,\mathbf{b}^T \f$ for two 6-vectors: one product per entry, no BLAS dispatch.
inline arma::mat::fixed<6,6> outer6(const double *a, const double *b) {
    arma::mat::fixed<6,6> C;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            C(i,j) = a[i]*b[j];
    return C;
}

/** @} */ // end of maths group

} //namespace simcoon
