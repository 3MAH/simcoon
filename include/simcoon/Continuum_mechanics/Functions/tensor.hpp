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
#include <string>
#include <Fastor/Fastor.h>

namespace simcoon {

// Forward declarations
class Rotation;
class tensor4;

/**
 * @brief Type tag for 2nd-order tensors, determines Voigt conversion factors and rotation rules.
 *
 * - stress:  Voigt = [s11, s22, s33, s12, s13, s23] — shear as-is
 * - strain:  Voigt = [e11, e22, e33, 2*e12, 2*e13, 2*e23] — shear doubled
 * - generic: same storage as stress (symmetric tensor, no physical convention)
 * - none:    non-symmetric tensor (e.g. F), .voigt() throws
 */
enum class VoigtType {
    stress,
    strain,
    generic,
    none
};

/**
 * @brief Type tag for 4th-order tensors, determines rotation rules and Voigt factor conventions.
 *
 * Rotation dispatch:
 * - stiffness:             QS * L * QS^T
 * - compliance:            QE * M * QE^T
 * - strain_concentration:  QE * A * QS^T
 * - stress_concentration:  QS * B * QE^T
 * - generic:               QS * C * QS^T (default)
 */
enum class Tensor4Type {
    stiffness,
    compliance,
    strain_concentration,
    stress_concentration,
    generic
};

/**
 * @brief A 2nd-order tensor with type tag for Voigt convention and rotation dispatch.
 *
 * Storage: arma::mat::fixed<3,3> (column-major, 72 bytes).
 * The Voigt vector is computed on the fly (no cache).
 */
class tensor2 {
private:
    arma::mat::fixed<3,3> _mat;
    VoigtType _vtype;
    mutable bool _symmetric;
    mutable bool _symmetry_checked;

public:
    // Constructors
    tensor2();
    explicit tensor2(VoigtType vtype);
    tensor2(const arma::mat::fixed<3,3> &m, VoigtType vtype);
    tensor2(const arma::mat &m, VoigtType vtype);

    // Construct from Voigt vector
    static tensor2 from_voigt(const arma::vec::fixed<6> &v, VoigtType vtype);
    static tensor2 from_voigt(const arma::vec &v, VoigtType vtype);

    // Static factories
    static tensor2 zeros(VoigtType vtype = VoigtType::stress);
    static tensor2 identity(VoigtType vtype = VoigtType::stress);

    // Copy/move
    tensor2(const tensor2 &other) = default;
    tensor2(tensor2 &&other) noexcept = default;
    tensor2& operator=(const tensor2 &other) = default;
    tensor2& operator=(tensor2 &&other) noexcept = default;
    ~tensor2() = default;

    // Accessors
    const arma::mat::fixed<3,3>& mat() const { return _mat; }
    arma::mat::fixed<3,3>& mat_mut();
    void set_mat(const arma::mat::fixed<3,3> &m);
    void set_mat(const arma::mat &m);
    void set_voigt(const arma::vec::fixed<6> &v);
    void set_voigt(const arma::vec &v);

    VoigtType vtype() const { return _vtype; }

    /**
     * @brief Compute Voigt vector on the fly (6 lookups + factor-2 on shear for strain).
     * @return vec::fixed<6> by value
     * @throws std::runtime_error if VoigtType::none
     */
    arma::vec::fixed<6> voigt() const;

    /**
     * @brief Zero-copy Fastor TensorMap wrapping the 3x3 matrix memory.
     */
    Fastor::TensorMap<double,3,3> fastor();
    Fastor::TensorMap<const double,3,3> fastor() const;

    // Symmetry
    bool is_symmetric(double tol = 1e-12) const;

    // Conversion helpers
    arma::mat to_arma_mat() const { return arma::mat(_mat); }
    arma::vec to_arma_voigt() const { return arma::vec(voigt()); }

    // Rotation
    tensor2 rotate(const Rotation &R, bool active = true) const;

    // Push-forward / pull-back (dispatch on VoigtType)
    tensor2 push_forward(const arma::mat::fixed<3,3> &F) const;
    tensor2 pull_back(const arma::mat::fixed<3,3> &F) const;

    // Arithmetic
    tensor2 operator+(const tensor2 &other) const;
    tensor2 operator-(const tensor2 &other) const;
    tensor2 operator*(double scalar) const;
    tensor2& operator+=(const tensor2 &other);
    tensor2& operator-=(const tensor2 &other);
    tensor2& operator*=(double scalar);
    friend tensor2 operator*(double scalar, const tensor2 &t);

    // Comparison
    bool operator==(const tensor2 &other) const;
};

// Free functions for tensor2
tensor2 stress(const arma::mat::fixed<3,3> &m);
tensor2 stress(const arma::vec::fixed<6> &v);
tensor2 strain(const arma::mat::fixed<3,3> &m);
tensor2 strain(const arma::vec::fixed<6> &v);
arma::vec::fixed<6> dev(const tensor2 &t);
double Mises(const tensor2 &t);
double trace(const tensor2 &t);

/**
 * @brief A 4th-order tensor with type tag for rotation dispatch and lazy Fastor cache.
 *
 * Storage: arma::mat::fixed<6,6> Voigt matrix (primary).
 * Fastor Tensor<double,3,3,3,3> is lazily computed and cached.
 */
class tensor4 {
private:
    arma::mat::fixed<6,6> _voigt;
    Tensor4Type _type;
    mutable Fastor::Tensor<double,3,3,3,3> _fastor;
    mutable bool _fastor_valid;

    void _invalidate_fastor() { _fastor_valid = false; }
    void _ensure_fastor() const;

public:
    // Constructors
    tensor4();
    explicit tensor4(Tensor4Type type);
    tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type);
    tensor4(const arma::mat &m, Tensor4Type type);

    // Static factories wrapping constitutive.hpp functions
    static tensor4 identity(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 volumetric(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 deviatoric(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 identity2(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 deviatoric2(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 zeros(Tensor4Type type = Tensor4Type::stiffness);

    // Copy/move
    tensor4(const tensor4 &other);
    tensor4(tensor4 &&other) noexcept;
    tensor4& operator=(const tensor4 &other);
    tensor4& operator=(tensor4 &&other) noexcept;
    ~tensor4() = default;

    // Accessors
    const arma::mat::fixed<6,6>& mat() const { return _voigt; }
    arma::mat::fixed<6,6>& mat_mut();
    void set_mat(const arma::mat::fixed<6,6> &m);
    void set_mat(const arma::mat &m);

    Tensor4Type type() const { return _type; }

    /**
     * @brief Get the cached Fastor 3x3x3x3 tensor (lazy recomputation if dirty).
     */
    const Fastor::Tensor<double,3,3,3,3>& fastor() const;

    // Conversion
    arma::mat to_arma_mat() const { return arma::mat(_voigt); }

    /**
     * @brief Contract with a tensor2: result_voigt = mat * t.voigt()
     *
     * Output VoigtType is inferred from Tensor4Type:
     * - stiffness (sigma = L : epsilon) -> stress
     * - compliance (epsilon = M : sigma) -> strain
     * - strain_concentration (epsilon = A : epsilon) -> strain
     * - stress_concentration (sigma = B : sigma) -> stress
     * - generic -> stress (default)
     */
    tensor2 contract(const tensor2 &t) const;

    /**
     * @brief Push-forward: C'_isrp = F_iL F_sJ F_rM F_pN C_LJMN
     */
    tensor4 push_forward(const arma::mat::fixed<3,3> &F) const;

    /**
     * @brief Pull-back: C'_LJMN = invF_lN invF_kM invF_jJ invF_iL C_ijkl
     */
    tensor4 pull_back(const arma::mat::fixed<3,3> &F) const;

    // Rotation
    tensor4 rotate(const Rotation &R, bool active = true) const;

    // Arithmetic
    tensor4 operator+(const tensor4 &other) const;
    tensor4 operator-(const tensor4 &other) const;
    tensor4 operator*(double scalar) const;
    tensor4& operator+=(const tensor4 &other);
    tensor4& operator-=(const tensor4 &other);
    tensor4& operator*=(double scalar);
    friend tensor4 operator*(double scalar, const tensor4 &t);

    // Comparison
    bool operator==(const tensor4 &other) const;
};

// Free functions for tensor4
tensor4 dyadic(const tensor2 &a, const tensor2 &b);
tensor4 auto_dyadic(const tensor2 &a);

} // namespace simcoon
