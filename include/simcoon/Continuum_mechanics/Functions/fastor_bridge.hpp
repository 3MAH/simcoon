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
    arma::mat::fixed<6,6> L;
    L.zeros();

    // Row 0: (0,0,*,*)
    L(0,0) = C(0,0,0,0);
    L(0,1) = C(0,0,1,1);
    L(0,2) = C(0,0,2,2);
    L(0,3) = 0.5*(C(0,0,0,1) + C(0,0,1,0));
    L(0,4) = 0.5*(C(0,0,0,2) + C(0,0,2,0));
    L(0,5) = 0.5*(C(0,0,1,2) + C(0,0,2,1));

    // Row 1: (1,1,*,*)
    L(1,0) = C(1,1,0,0);
    L(1,1) = C(1,1,1,1);
    L(1,2) = C(1,1,2,2);
    L(1,3) = 0.5*(C(1,1,0,1) + C(1,1,1,0));
    L(1,4) = 0.5*(C(1,1,0,2) + C(1,1,2,0));
    L(1,5) = 0.5*(C(1,1,1,2) + C(1,1,2,1));

    // Row 2: (2,2,*,*)
    L(2,0) = C(2,2,0,0);
    L(2,1) = C(2,2,1,1);
    L(2,2) = C(2,2,2,2);
    L(2,3) = 0.5*(C(2,2,0,1) + C(2,2,1,0));
    L(2,4) = 0.5*(C(2,2,0,2) + C(2,2,2,0));
    L(2,5) = 0.5*(C(2,2,1,2) + C(2,2,2,1));

    // Row 3: (0,1,*,*)
    L(3,0) = C(0,1,0,0);
    L(3,1) = C(0,1,1,1);
    L(3,2) = C(0,1,2,2);
    L(3,3) = 0.5*(C(0,1,0,1) + C(0,1,1,0));
    L(3,4) = 0.5*(C(0,1,0,2) + C(0,1,2,0));
    L(3,5) = 0.5*(C(0,1,1,2) + C(0,1,2,1));

    // Row 4: (0,2,*,*)
    L(4,0) = C(0,2,0,0);
    L(4,1) = C(0,2,1,1);
    L(4,2) = C(0,2,2,2);
    L(4,3) = 0.5*(C(0,2,0,1) + C(0,2,1,0));
    L(4,4) = 0.5*(C(0,2,0,2) + C(0,2,2,0));
    L(4,5) = 0.5*(C(0,2,1,2) + C(0,2,2,1));

    // Row 5: (1,2,*,*)
    L(5,0) = C(1,2,0,0);
    L(5,1) = C(1,2,1,1);
    L(5,2) = C(1,2,2,2);
    L(5,3) = 0.5*(C(1,2,0,1) + C(1,2,1,0));
    L(5,4) = 0.5*(C(1,2,0,2) + C(1,2,2,0));
    L(5,5) = 0.5*(C(1,2,1,2) + C(1,2,2,1));

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
    Fastor::Tensor<double,3,3,3,3> result;
    result.zeros();

    for (int i = 0; i < 3; ++i)
    for (int s = 0; s < 3; ++s)
    for (int r = 0; r < 3; ++r)
    for (int p = 0; p < 3; ++p) {
        double sum = 0.0;
        for (int L = 0; L < 3; ++L)
        for (int J = 0; J < 3; ++J)
        for (int M = 0; M < 3; ++M)
        for (int N = 0; N < 3; ++N) {
            sum += F(i,L) * F(s,J) * F(r,M) * F(p,N) * C(L,J,M,N);
        }
        result(i,s,r,p) = sum;
    }
    return result;
}

/**
 * @brief Pull-back a 4th-order tensor via Fastor: C'_LJMN = F^-1_lN F^-1_kM F^-1_jJ F^-1_iL C_ijkl
 */
inline Fastor::Tensor<double,3,3,3,3> pull_back_4(
    const Fastor::Tensor<double,3,3,3,3> &C,
    const Fastor::Tensor<double,3,3> &invF)
{
    // Pull-back is just push-forward with F^{-1}
    return push_forward_4(C, invF);
}

} // namespace simcoon
