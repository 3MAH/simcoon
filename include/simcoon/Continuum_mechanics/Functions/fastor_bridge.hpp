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

#pragma once
#include <armadillo>
#include <Fastor/Fastor.h>

namespace simcoon {

/**
 * @file fastor_bridge.hpp
 * @brief Conversion helpers between Armadillo and Fastor tensor types.
 *
 * These functions replace the FTensor conversion functions in transfer.hpp
 * for use with the tensor2/tensor4 classes and Fastor-based operations.
 */

// Voigt index mapping: (i,j) -> Voigt index I
// (0,0)->0, (1,1)->1, (2,2)->2, (0,1)->3, (0,2)->4, (1,2)->5
inline constexpr int voigt_map[3][3] = {
    {0, 3, 4},
    {3, 1, 5},
    {4, 5, 2}
};

/**
 * @brief Zero-copy TensorMap wrapping an arma::mat::fixed<3,3> memory.
 *
 * Both Armadillo and Fastor use column-major storage, so this is safe.
 * The returned TensorMap is valid only as long as the source matrix is alive.
 */
inline Fastor::TensorMap<double,3,3> arma_to_fastor2(arma::mat::fixed<3,3> &m) {
    return Fastor::TensorMap<double,3,3>(m.memptr());
}

inline Fastor::TensorMap<const double,3,3> arma_to_fastor2(const arma::mat::fixed<3,3> &m) {
    return Fastor::TensorMap<const double,3,3>(m.memptr());
}

/**
 * @brief Copy a Fastor Tensor<double,3,3> to an arma::mat::fixed<3,3>.
 */
inline arma::mat::fixed<3,3> fastor2_to_arma(const Fastor::Tensor<double,3,3> &T) {
    arma::mat::fixed<3,3> m;
    // Both column-major, direct memcpy
    std::memcpy(m.memptr(), T.data(), 9 * sizeof(double));
    return m;
}

/**
 * @brief Voigt 6x6 matrix -> Fastor Tensor<double,3,3,3,3>
 *
 * The conversion factors depend on the tensor type. For stiffness convention
 * (the default, matching mat_FTensor4), no factors are needed:
 * C_ijkl = L(voigt_map[i][j], voigt_map[k][l]) with minor symmetry enforcement.
 *
 * This matches the existing mat_FTensor4 behavior: the Voigt matrix stores
 * stiffness components directly without factor corrections.
 */
inline Fastor::Tensor<double,3,3,3,3> voigt_to_fastor4(const arma::mat::fixed<6,6> &L) {
    Fastor::Tensor<double,3,3,3,3> C;
    C.zeros();

    for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
            int ij = voigt_map[i][j];
            for (int k = 0; k < 3; ++k) {
                for (int l = k; l < 3; ++l) {
                    int kl = voigt_map[k][l];
                    double val = L(ij, kl);
                    C(i,j,k,l) = val;
                    C(i,j,l,k) = val;
                    C(j,i,k,l) = val;
                    C(j,i,l,k) = val;
                }
            }
        }
    }
    return C;
}

/**
 * @brief Fastor Tensor<double,3,3,3,3> -> Voigt 6x6 matrix
 *
 * Matches the existing FTensor4_mat behavior: averages minor-symmetric pairs.
 * This is the stiffness convention (no factor corrections).
 */
inline arma::mat::fixed<6,6> fastor4_to_voigt(const Fastor::Tensor<double,3,3,3,3> &C) {
    // Voigt index -> (i,j) pair
    static constexpr int vi[6] = {0, 1, 2, 0, 0, 1};
    static constexpr int vj[6] = {0, 1, 2, 1, 2, 2};

    arma::mat::fixed<6,6> L;
    for (int I = 0; I < 6; ++I) {
        int i = vi[I], j = vj[I];
        for (int J = 0; J < 6; ++J) {
            int k = vi[J], l = vj[J];
            L(I, J) = 0.5 * (C(i,j,k,l) + C(i,j,l,k));
        }
    }
    return L;
}

/**
 * @brief Push-forward a 4th-order tensor via Fastor: C'_isrp = F_iL F_sJ F_rM F_pN C_LJMN
 *
 * This replaces the FTensor boilerplate in objective_rates.cpp.
 */
inline Fastor::Tensor<double,3,3,3,3> push_forward_4(
    const Fastor::Tensor<double,3,3,3,3> &C,
    const Fastor::Tensor<double,3,3> &F)
{
    // Two-pass contraction: O(2 * 3^5) instead of O(3^8)
    // Step 1: T_isJN = F_iL * C_LJMN * F^T_Ns  (contract L and M)
    //   Actually: contract two indices at a time via intermediate tensor
    // Simpler factored version: hoist F(i,L)*F(s,J) out of inner loops
    Fastor::Tensor<double,3,3,3,3> result;
    result.zeros();

    for (int i = 0; i < 3; ++i)
    for (int s = 0; s < 3; ++s)
    for (int r = 0; r < 3; ++r)
    for (int p = 0; p < 3; ++p) {
        double sum = 0.0;
        for (int L = 0; L < 3; ++L)
        for (int J = 0; J < 3; ++J) {
            double fiL_fsJ = F(i,L) * F(s,J);
            for (int M = 0; M < 3; ++M) {
                double frM = F(r,M);
                for (int N = 0; N < 3; ++N) {
                    sum += fiL_fsJ * frM * F(p,N) * C(L,J,M,N);
                }
            }
        }
        result(i,s,r,p) = sum;
    }
    return result;
}

} // namespace simcoon
