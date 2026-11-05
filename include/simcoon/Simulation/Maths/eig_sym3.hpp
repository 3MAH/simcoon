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

///@file eig_sym3.hpp
///@brief Fixed-size spectral decomposition of symmetric 3x3 matrices and the spectral functions built on it
///@version 1.0

#pragma once
#include <armadillo>

namespace simcoon{

/** @addtogroup maths
 *  @{
 */

/**
 * @brief Eigen-decomposition of a symmetric 3x3 matrix, \f$ \mathbf{A} = \mathbf{Q}\,\mathrm{diag}(\boldsymbol{\lambda})\,\mathbf{Q}^T \f$.
 *
 * Cyclic Jacobi rotations (Kopp 2008, dsyevj3) on fixed-size storage: no heap allocation and no
 * LAPACK call, about three times cheaper than arma::eig_sym (dsyevd) on a 3x3, and accurate to
 * the last digits for every eigenvalue, the small ones of a stiff SPD matrix included. Only the
 * upper triangle of @p A is read. The eigenvalues come in ascending order, as with
 * arma::eig_sym; the eigenvectors are the orthonormal columns of @p Q (sign arbitrary). A
 * diagonal @p A is returned exactly (@p Q a permutation).
 *
 * This is the solver behind the per-point spectral functions of the finite-strain kinematics
 * (polar decompositions, logarithmic strains, Daleckii-Krein derivatives).
 *
 * @param[out] lam eigenvalues, ascending
 * @param[out] Q eigenvectors (columns)
 * @param[in] A symmetric matrix (upper triangle read)
 * @return false when @p A is not finite or the iteration does not converge (as arma::eig_sym)
 */
bool eig_sym_3x3(arma::vec::fixed<3> &lam, arma::mat::fixed<3,3> &Q, const arma::mat::fixed<3,3> &A);

/**
 * @brief Symmetric matrix of given eigenvalues on a given eigenbasis, \f$ \mathbf{Q}\,\mathrm{diag}(\mathbf{f})\,\mathbf{Q}^T \f$.
 *
 * The upper triangle is computed and mirrored: the result is exactly symmetric.
 *
 * @param Q eigenvectors (columns)
 * @param f eigenvalues of the result
 * @return the 3x3 symmetric matrix
 */
arma::mat::fixed<3,3> spectral_compose_3x3(const arma::mat::fixed<3,3> &Q, const arma::vec::fixed<3> &f);

/**
 * @brief Spectral congruence of a 3x3 matrix on an eigenbasis,
 * \f$ \sum_{a,b} K_{ab}\,\mathbf{P}_a \mathbf{X} \mathbf{P}_b
 *     = \mathbf{N}\left(\mathbf{K} \circ \mathbf{N}^T \mathbf{X} \mathbf{N}\right)\mathbf{N}^T \f$.
 *
 * \f$ \mathbf{P}_a = \mathbf{n}_a \otimes \mathbf{n}_a \f$ are the eigenprojectors of the
 * orthonormal columns of @p N and \f$ \circ \f$ the entrywise product: the form of the
 * Daleckii-Krein derivatives of isotropic matrix functions, with \f$ K_{ab} \f$ their divided
 * differences.
 *
 * @param N eigenvectors (orthonormal columns), e.g. from eig_sym_3x3
 * @param K coefficients \f$ K_{ab} \f$
 * @param X the 3x3 matrix
 * @return \f$ \mathbf{N}(\mathbf{K} \circ \mathbf{N}^T \mathbf{X} \mathbf{N})\mathbf{N}^T \f$
 */
inline arma::mat::fixed<3,3> spectral_congruence(const arma::mat::fixed<3,3> &N, const arma::mat::fixed<3,3> &K,
                                                 const arma::mat::fixed<3,3> &X) {
    const arma::mat::fixed<3,3> Xp = N.t()*X*N;
    return N*(K % Xp)*N.t();
}

/**
 * @brief Square root of a symmetric positive semi-definite 3x3 matrix (as arma::sqrtmat_sympd).
 * @param[out] out the square root
 * @param[in] A symmetric matrix
 * @return false when the decomposition fails or an eigenvalue is negative
 */
bool sqrtmat_sympd_3x3(arma::mat::fixed<3,3> &out, const arma::mat::fixed<3,3> &A);

/**
 * @brief Logarithm of a symmetric positive definite 3x3 matrix (as arma::logmat_sympd).
 * @param[out] out the logarithm
 * @param[in] A symmetric matrix
 * @return false when the decomposition fails or an eigenvalue is not positive
 */
bool logmat_sympd_3x3(arma::mat::fixed<3,3> &out, const arma::mat::fixed<3,3> &A);

/**
 * @brief Exponential of a symmetric 3x3 matrix (as arma::expmat_sym).
 * @param[out] out the exponential
 * @param[in] A symmetric matrix
 * @return false when the decomposition fails
 */
bool expmat_sym_3x3(arma::mat::fixed<3,3> &out, const arma::mat::fixed<3,3> &A);

/** @} */ // end of maths group

} //namespace simcoon
