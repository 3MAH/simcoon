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

///@file derivatives.hpp
///@brief A set of functions that performs the derivative of specific tensors attributes (Invariants, determinant, inverse...)
///@version 1.0

#pragma once
#include <armadillo>
#include <simcoon/FTensor.hpp>

namespace simcoon{

/**
 * @file derivatives.hpp
 * @brief Tensor derivative functions for continuum mechanics
 * @author Yves Chemisky
 * @version 1.0
 *
 * This file provides functions to compute derivatives of tensor invariants,
 * trace, determinant, and inverse with respect to the tensor components.
 * These derivatives are essential for constitutive model development,
 * particularly for hyperelastic materials and return mapping algorithms.
 */

/** @addtogroup functions
 *  @{
 */

/**
 * @brief Computes the derivative of the first invariant (trace) with respect to the tensor
 *
 * @details The first invariant of a tensor \f$ \mathbf{S} \f$ is the trace:
 * \f[
 * I_1 = \text{tr}(\mathbf{S}) = S_{ii}
 * \f]
 *
 * The derivative with respect to the tensor components is:
 * \f[
 * \frac{\partial I_1}{\partial \mathbf{S}} = \mathbf{I}
 * \f]
 * where \f$ \mathbf{I} \f$ is the second-order identity tensor.
 *
 * @param S Input tensor (3×3 matrix)
 * @return 3×3 matrix representing the derivative (identity tensor)
 *
 * @note For a symmetric tensor, this returns the identity matrix
 *
 * @code
 * arma::mat S = {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}};
 * arma::mat dI1 = dI1DS(S);  // Returns identity matrix
 * @endcode
 */
arma::mat dI1DS(const arma::mat &S);

/**
 * @brief Computes the derivative of the second invariant with respect to the tensor
 *
 * @details The second invariant is defined as:
 * \f[
 * I_2 = \frac{1}{2} S_{ij} S_{ij} = \frac{1}{2} \mathbf{S} : \mathbf{S}
 * \f]
 *
 * The derivative with respect to the tensor components is:
 * \f[
 * \frac{\partial I_2}{\partial \mathbf{S}} = \mathbf{S}
 * \f]
 *
 * @param S Input tensor (3×3 matrix)
 * @return 3×3 matrix representing the derivative (the tensor itself)
 *
 * @note This invariant definition differs from the classical second principal invariant.
 *       Here \f$ I_2 = \frac{1}{2} \|\mathbf{S}\|^2 \f$ (squared Frobenius norm divided by 2)
 *
 * @code
 * arma::mat S = {{1, 2, 0}, {2, 3, 0}, {0, 0, 1}};
 * arma::mat dI2 = dI2DS(S);  // Returns S
 * @endcode
 */
arma::mat dI2DS(const arma::mat &S);

/**
 * @brief Computes the derivative of the third invariant with respect to the tensor
 *
 * @details The third invariant is defined as:
 * \f[
 * I_3 = \frac{1}{3} S_{ij} S_{jk} S_{ki} = \frac{1}{3} \text{tr}(\mathbf{S}^3)
 * \f]
 *
 * The derivative with respect to the tensor components is:
 * \f[
 * \frac{\partial I_3}{\partial \mathbf{S}} = \mathbf{S}^2
 * \f]
 *
 * @param S Input tensor (3×3 matrix)
 * @return 3×3 matrix representing the derivative (\f$ \mathbf{S}^2 \f$)
 *
 * @note This invariant definition differs from the determinant, which is often selected as the third principal invariant (determinant).
 *       Here \f$ I_3 = \frac{1}{3} \text{tr}(\mathbf{S}^3) \f$
 *
 * @code
 * arma::mat S = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
 * arma::mat dI3 = dI3DS(S);  // Returns S^2 = {{4,0,0},{0,4,0},{0,0,4}}
 * @endcode
 */
arma::mat dI3DS(const arma::mat &S);

/**
 * @brief Computes the derivative of the second invariant of the deviatoric tensor
 *
 * @details The second invariant of the deviatoric part of the tensor \f$\mathbf{S}\f$ is
 * \f[
 * J_2 = \frac{1}{2} S_{ij} S_{ij}
 * \f]
 * Its derivative with respect to \f$\mathbf{S}\f$ in Voigt notation is simply the deviatoric part:
 * \f[
 * \frac{\partial J_2}{\partial \mathbf{S}} = \mathbf{S}^{\text{dev}}
 * \f]
 *
 * @param S Input 3×3 stress tensor
 * @return 3×3 matrix representing the derivative of J2
 *
 * @see dev() for computing the deviatoric part of a tensor
 *
 * @code
 * arma::mat S = arma::randu(3,3);
 * arma::mat dJ2 = dJ2DS(S);  // Returns the deviatoric part of S
 * @endcode
 */
arma::mat dJ2DS(const arma::mat &S);

/**
 * @brief Computes the derivative of the third invariant of the deviatoric tensor
 *
 * @details The third invariant of the deviatoric part of the tensor \f$\mathbf{S}\f$ is
 * \f[
 * J_3 = \frac{1}{3} S_{ij} S_{jk} S_{ki}
 * \f]
 * Its derivative with respect to \f$\mathbf{S}\f$ in Voigt notation is:
 * \f[
 * \frac{\partial J_3}{\partial \mathbf{S}} = \mathbf{S}^2 - \frac{1}{3} \text{tr}(\mathbf{S}^2) \mathbf{I}
 * \f]
 * where \f$\mathbf{S}\f$ is the deviatoric part of the input tensor.
 *
 * @param S Input 3×3 stress tensor
 * @return 3×3 matrix representing the derivative of J3
 *
 * @see dev() for computing the deviatoric part of a tensor
 *
 * @code
 * arma::mat S = arma::randu(3,3);
 * arma::mat dJ3 = dJ3DS(S);  // Returns the derivative of the third invariant
 * @endcode
 */
arma::mat dJ3DS(const arma::mat &S);
/**
 * @brief Computes the derivative of the trace with respect to the tensor
 *
 * @details This is equivalent to dI1DS(). The trace derivative is:
 * \f[
 * \frac{\partial \text{tr}(\mathbf{S})}{\partial \mathbf{S}} = \mathbf{I}
 * \f]
 *
 * @param S Input tensor (3×3 matrix)
 * @return 3×3 identity matrix
 *
 * @see dI1DS() for equivalent functionality
 *
 * @code
 * arma::mat S = arma::randu(3,3);
 * arma::mat dtr = dtrSdS(S);  // Returns identity matrix
 * @endcode
 */
arma::mat dtrSdS(const arma::mat &S);

/**
 * @brief Computes the derivative of the determinant with respect to the tensor
 *
 * @details For a tensor \f$ \mathbf{S} \f$, the derivative of the determinant is:
 * \f[
 * \frac{\partial \det(\mathbf{S})}{\partial \mathbf{S}} = \det(\mathbf{S}) \mathbf{S}^{-T}
 * \f]
 *
 * Using the cofactor matrix \f$ \text{cof}(\mathbf{S}) = \det(\mathbf{S}) \mathbf{S}^{-T} \f$:
 * \f[
 * \frac{\partial \det(\mathbf{S})}{\partial \mathbf{S}} = \text{cof}(\mathbf{S})
 * \f]
 *
 * @param S Input tensor (3×3 matrix), must be invertible
 * @return 3×3 matrix representing the cofactor matrix
 *
 * @warning The input tensor must be invertible (non-zero determinant)
 *
 * @note This derivative is essential for hyperelastic models involving \f$ J = \det(\mathbf{F}) \f$
 *
 * @code
 * arma::mat F = {{1.1, 0.1, 0}, {0, 1.0, 0}, {0, 0, 1.0}};
 * arma::mat ddet = ddetSdS(F);  // Returns cofactor matrix
 * @endcode
 */
arma::mat ddetSdS(const arma::mat &S);

/**
 * @brief Computes the derivative of the inverse of a symmetric tensor
 *
 * @details For a symmetric tensor \f$ \mathbf{S} \f$, the derivative of the inverse is:
 * \f[
 * \frac{\partial \mathbf{S}^{-1}}{\partial S_{kl}} = -S^{-1}_{ik} S^{-1}_{lj}
 * \f]
 *
 * In index notation, for symmetric tensors:
 * \f[
 * \frac{\partial S^{-1}_{ij}}{\partial S_{kl}} = -\frac{1}{2}\left( S^{-1}_{ik} S^{-1}_{jl} + S^{-1}_{il} S^{-1}_{jk} \right)
 * \f]
 *
 * @param S Input symmetric tensor (3×3 matrix), must be invertible
 * @return 3×3 matrix containing the derivative (returned in a specific contracted form)
 *
 * @warning The input tensor must be symmetric and invertible
 *
 * @note This function is optimized for symmetric tensors. For general tensors,
 *       the derivative is a fourth-order tensor.
 *
 * @code
 * arma::mat S = {{2, 1, 0}, {1, 3, 0}, {0, 0, 1}};  // Symmetric
 * arma::mat dinv = dinvSdSsym(S);
 * @endcode
 *
 * @see Holzapfel, G.A. (2000). Nonlinear Solid Mechanics, Wiley. Section 1.6.3
 */
arma::mat dinvSdSsym(const arma::mat &S);

/** @} */ // end of functions group

} //namespace simcoon
