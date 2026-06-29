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
#include <optional>
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
 * @brief Type tag for 4th-order tensors: selects the engineering<->Mandel congruence,
 *        the contraction output type, and the inverse type swap.
 *
 * Storage is Kelvin-Mandel internally (see tensor4). The tag fixes how the engineering
 * Voigt matrix passed to the constructor maps to the stored Mandel matrix (and back via
 * mat()), with \f$ N = \mathrm{diag}(1,1,1,\sqrt2,\sqrt2,\sqrt2) \f$:
 * - stiffness / generic:    \f$ X_m = N\,X_{eng}\,N \f$            (\f$\sigma = L:\varepsilon\f$)
 * - compliance:             \f$ X_m = N^{-1} X_{eng}\,N^{-1} \f$   (\f$\varepsilon = M:\sigma\f$)
 * - strain_concentration:   \f$ X_m = N^{-1} X_{eng}\,N \f$        (\f$\varepsilon = A:\varepsilon\f$)
 * - stress_concentration:   \f$ X_m = N\,X_{eng}\,N^{-1} \f$       (\f$\sigma = B:\sigma\f$)
 *
 * In Mandel every type shares one identity (eye(6)), one orthogonal rotation
 * (R6 * X * R6^T) and a plain matrix inverse; the engineering asymmetry lives only in
 * the boundary congruence above.
 */
enum class Tensor4Type {
    stiffness,
    compliance,
    strain_concentration,
    stress_concentration,
    generic
};

/**
 * @brief Corotational rate type for tangent modulus push-forward.
 *
 * Controls the B-correction applied after the F-push-forward of a stiffness
 * tangent to make the spatial tangent consistent with a given objective rate
 * of Kirchhoff stress:
 * - lie:            Pure push-forward (Truesdell/Oldroyd rate), no correction
 * - jaumann:        Jaumann (co-rotational) W-spin correction
 * - green_naghdi:   Green-Naghdi spin correction (R from polar decomposition)
 * - logarithmic:    BXM logarithmic rate correction (default in simcoon solver)
 * - logarithmic_R:  Logarithmic, R-transport framework
 * - logarithmic_F:  Logarithmic, F-transport framework
 *
 * @note `logarithmic`, `logarithmic_R`, and `logarithmic_F` share the same
 *       algorithmic tangent here: they all reduce to the F-push-forward of the
 *       material tangent followed by the same B^(4) correction (cf. Chemisky
 *       et al. preprint, Sec. "Two integrable logarithmic frameworks": both
 *       integrable frameworks share B^(4); the difference lies in the *stress
 *       integrator*'s transport operator). Their genuinely distinct stress
 *       updates live in Continuum_mechanics/Functions/objective_rates.cpp
 *       (`logarithmic`, `logarithmic_R`, `logarithmic_F`), not here.
 */
enum class CoRate {
    lie,
    jaumann,
    green_naghdi,
    logarithmic,
    logarithmic_R,
    logarithmic_F
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

public:
    // Constructors
    tensor2();
    explicit tensor2(VoigtType vtype);
    tensor2(const arma::mat::fixed<3,3> &m, VoigtType vtype);
    tensor2(const arma::mat &m, VoigtType vtype);

    // Construct from Voigt vector
    static tensor2 from_voigt(const arma::vec::fixed<6> &v, VoigtType vtype);
    static tensor2 from_voigt(const arma::vec &v, VoigtType vtype);
    static tensor2 from_voigt(const arma::vec &v, const std::string &type_str);

    // Construct from a Kelvin-Mandel vector (shear scaled by sqrt(2), identical for stress/strain).
    static tensor2 from_mandel(const arma::vec::fixed<6> &v, VoigtType vtype);

    // Static factories
    static tensor2 zeros(VoigtType vtype = VoigtType::stress);
    static tensor2 zeros(const std::string &type_str);
    static tensor2 identity(VoigtType vtype = VoigtType::stress);
    static tensor2 identity(const std::string &type_str);

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
     * @brief Kelvin-Mandel 6-vector: \f$[t_{11},t_{22},t_{33},\sqrt2\,t_{12},\sqrt2\,t_{13},\sqrt2\,t_{23}]\f$.
     *
     * Convention-symmetric (identical for stress and strain); the natural input to a
     * tensor4 Mandel contraction.
     */
    arma::vec::fixed<6> mandel() const;

    /**
     * @brief Convert to Fastor::Tensor<double,3,3>, handling row/col-major layout.
     */
    Fastor::Tensor<double,3,3> fastor() const;

    // Symmetry
    bool is_symmetric(double tol = 1e-12) const;

    // Conversion helpers
    arma::mat to_arma_mat() const { return arma::mat(_mat); }
    arma::vec to_arma_voigt() const { return arma::vec(voigt()); }

    // Rotation
    tensor2 rotate(const Rotation &R, bool active = true) const;

    // Push-forward / pull-back (dispatch on VoigtType)
    // metric=true (default): includes J=det(F) factor (Piola transformation)
    // metric=false: pure geometric transport
    tensor2 push_forward(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    tensor2 push_forward(const arma::mat &F, bool metric = true) const;
    tensor2 pull_back(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    tensor2 pull_back(const arma::mat &F, bool metric = true) const;

    // Arithmetic
    tensor2 operator+(const tensor2 &other) const;
    tensor2 operator-(const tensor2 &other) const;
    tensor2 operator-() const;   // unary minus
    tensor2 operator*(double scalar) const;
    tensor2 operator/(double scalar) const;
    tensor2& operator+=(const tensor2 &other);
    tensor2& operator-=(const tensor2 &other);
    tensor2& operator*=(double scalar);
    tensor2& operator/=(double scalar);
    friend tensor2 operator*(double scalar, const tensor2 &t);

    /// Element-wise (Schur) product on the 3x3 matrix, like Armadillo %.
    /// Use sum(a % b) for the double contraction scalar.
    tensor2 operator%(const tensor2 &other) const;

    /// Element-wise division on the 3x3 matrix.
    tensor2 operator/(const tensor2 &other) const;

    // Comparison
    bool operator==(const tensor2 &other) const;
    bool operator!=(const tensor2 &other) const;
};

// Free functions for tensor2
tensor2 stress(const arma::mat::fixed<3,3> &m);
tensor2 stress(const arma::vec::fixed<6> &v);
tensor2 stress(const arma::mat &m);
tensor2 stress(const arma::vec &v);
tensor2 strain(const arma::mat::fixed<3,3> &m);
tensor2 strain(const arma::vec::fixed<6> &v);
tensor2 strain(const arma::mat &m);
tensor2 strain(const arma::vec &v);
tensor2 dev(const tensor2 &t);
double Mises(const tensor2 &t);
double trace(const tensor2 &t);

// Reduction free functions (Armadillo-style)
double sum(const tensor2 &t);
double accu(const tensor2 &t);
double norm(const tensor2 &t);       // Frobenius norm
double det(const tensor2 &t);
tensor2 abs(const tensor2 &t);
tensor2 trans(const tensor2 &t);     // transpose

/**
 * @brief A 4th-order tensor stored in the Kelvin-Mandel convention, with a type tag and lazy Fastor cache.
 *
 * Storage: arma::mat::fixed<6,6> in Mandel (sqrt2 on shear rows/cols), so matrix algebra
 * equals tensor algebra (identity = eye(6), inverse = arma::inv, rotation = orthogonal
 * R6 * X * R6^T). The engineering Voigt form is exposed at the boundary via mat()/set_mat();
 * the per-type engineering<->Mandel congruence is documented on Tensor4Type.
 * Fastor Tensor<double,3,3,3,3> (the convention-free full-index tensor) is lazily computed
 * and cached for push_forward/pull_back.
 *
 * @warning NOT thread-safe. The mutable Fastor cache is lazily populated by const
 *          methods (fastor(), push_forward(), pull_back()) without synchronization.
 *          Do not share a single tensor4 instance across threads without external
 *          locking. Creating separate tensor4 objects per thread is safe.
 */
class tensor4 {
private:
    arma::mat::fixed<6,6> _mandel;   ///< Kelvin-Mandel 6x6 (sqrt2 on shear); engineering form via mat()
    Tensor4Type _type;
    mutable std::optional<Fastor::Tensor<double,3,3,3,3>> _fastor;

    void _invalidate_fastor() { _fastor.reset(); }
    void _ensure_fastor() const;

    /// Internal: wrap an already-Mandel 6x6 directly (skips the engineering->Mandel congruence).
    static tensor4 from_mandel(const arma::mat::fixed<6,6> &m_mandel, Tensor4Type type);

public:
    // Constructors
    tensor4();
    explicit tensor4(Tensor4Type type);
    tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type);
    tensor4(const arma::mat &m, Tensor4Type type);
    tensor4(const arma::mat &m, const std::string &type_str);

    // Static projector factories. Stored in Mandel, every type shares identity = eye(6) and
    // deviatoric = eye(6) - volumetric. The engineering form via mat() recovers the classic
    // convention: identity() reads back as Ireal (stiffness), Ireal2 (compliance), eye(6)
    // (strain/stress_concentration).
    static tensor4 identity(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 volumetric(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 deviatoric(Tensor4Type type = Tensor4Type::stiffness);
    static tensor4 zeros(Tensor4Type type = Tensor4Type::stiffness);

    // Copy/move
    tensor4(const tensor4 &other);
    tensor4(tensor4 &&other) noexcept;
    tensor4& operator=(const tensor4 &other);
    tensor4& operator=(tensor4 &&other) noexcept;
    ~tensor4() = default;

    // Accessors. mat() returns the engineering Voigt 6x6 by value (converted from internal
    // Mandel); set_mat() takes the engineering form. mat_mut() is intentionally absent: there
    // is no coherent mutable engineering view over Mandel storage -- use set_mat().
    arma::mat::fixed<6,6> mat() const;
    void set_mat(const arma::mat::fixed<6,6> &m);
    void set_mat(const arma::mat &m);

    Tensor4Type type() const { return _type; }

    /**
     * @brief Get the cached Fastor 3x3x3x3 tensor (lazy recomputation if dirty).
     */
    const Fastor::Tensor<double,3,3,3,3>& fastor() const;

    // Conversion
    arma::mat to_arma_mat() const { return arma::mat(mat()); }

    /**
     * @brief Contract with a tensor2 (double contraction), e.g. sigma = L : eps.
     *        Internally a plain Mandel matrix-vector product.
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
     * @brief Push-forward to the spatial configuration. Kernel branches on Tensor4Type:
     *   - stiffness / generic (contravariant):
     *       C'_isrp = F_iL F_sJ F_rM F_pN C_LJMN, scaled by 1/J if metric=true
     *   - compliance (covariant):
     *       M'_isrp = (F^{-T})_iL (F^{-T})_sJ (F^{-T})_rM (F^{-T})_pN M_LJMN, scaled by J
     *   - strain_concentration / stress_concentration: not implemented (mixed indices)
     *
     * The compliance kernel is what makes inverse(K).push_forward(F) ≡ inverse(K.push_forward(F))
     * for non-orthogonal F.
     *
     * @param metric If true (default), includes the J=det(F) Piola scaling
     * @throws std::runtime_error for concentration tensor types
     */
    tensor4 push_forward(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    tensor4 push_forward(const arma::mat &F, bool metric = true) const;

    /**
     * @brief Push-forward with corotational rate correction.
     *
     * Computes the spatial tangent consistent with the selected objective rate.
     * For CoRate::lie this is equivalent to the plain push_forward.
     * Other rates apply a correction that requires the Kirchhoff stress tau.
     *
     * @param F    Deformation gradient (3x3)
     * @param rate Corotational rate type
     * @param tau  Kirchhoff stress (needed for all rates except lie)
     * @param metric If true, includes J factor
     */
    tensor4 push_forward(const arma::mat::fixed<3,3> &F, CoRate rate,
                         const tensor2 &tau, bool metric = true) const;
    tensor4 push_forward(const arma::mat &F, CoRate rate,
                         const tensor2 &tau, bool metric = true) const;

    /**
     * @brief Pull-back to the reference configuration. Inverts push_forward(F):
     *   - stiffness / generic: kernel is F^{-1}, scaled by J if metric=true
     *   - compliance:         kernel is F^T,    scaled by 1/J
     *   - concentration types: not implemented (mixed indices)
     *
     * @param metric If true (default), includes the J=det(F) Piola scaling
     * @throws std::runtime_error for concentration tensor types
     */
    tensor4 pull_back(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    tensor4 pull_back(const arma::mat &F, bool metric = true) const;

    /**
     * @brief Invert the 6x6 Voigt matrix.
     *
     * Type inference: stiffness <-> compliance, others stay the same.
     */
    tensor4 inverse() const;

    // Rotation
    tensor4 rotate(const Rotation &R, bool active = true) const;

    // Arithmetic
    tensor4 operator+(const tensor4 &other) const;
    tensor4 operator-(const tensor4 &other) const;
    tensor4 operator-() const;   // unary minus
    tensor2 operator*(const tensor2 &t) const;   // contraction: L * eps → sigma
    tensor4 operator*(double scalar) const;
    tensor4 operator/(double scalar) const;
    tensor4& operator+=(const tensor4 &other);
    tensor4& operator-=(const tensor4 &other);
    tensor4& operator*=(double scalar);
    tensor4& operator/=(double scalar);
    friend tensor4 operator*(double scalar, const tensor4 &t);

    /// Element-wise (Schur) product on the 6x6 Voigt matrix.
    tensor4 operator%(const tensor4 &other) const;

    // Comparison
    bool operator==(const tensor4 &other) const;
    bool operator!=(const tensor4 &other) const;
};

// Free functions for tensor4
tensor4 dyadic(const tensor2 &a, const tensor2 &b);
tensor4 auto_dyadic(const tensor2 &a);
tensor4 sym_dyadic(const tensor2 &a, const tensor2 &b);
tensor4 auto_sym_dyadic(const tensor2 &a);

// ============================================================================
// Batch operations
// ============================================================================
// Convention: voigt (6, N) columns, matrices/tensors (3,3,N) or (6,6,N) cubes.
// rot/F with n_slices==1 are broadcast to all N points.

/// Batch rotate N tensor2 objects. voigt:(6,N), rot:(3,3,N_r).
arma::mat batch_rotate(const arma::mat &voigt, VoigtType vtype,
                       const arma::cube &rot_matrices, bool active = true);

/// Batch push-forward N tensor2. voigt:(6,N), F:(3,3,N_f).
arma::mat batch_push_forward(const arma::mat &voigt, VoigtType vtype,
                             const arma::cube &F, bool metric = true);

/// Batch pull-back N tensor2. voigt:(6,N), F:(3,3,N_f).
arma::mat batch_pull_back(const arma::mat &voigt, VoigtType vtype,
                          const arma::cube &F, bool metric = true);

/// Batch von Mises for N tensor2. voigt:(6,N) → (N).
arma::vec batch_mises(const arma::mat &voigt, VoigtType vtype);

/// Batch trace for N tensor2. voigt:(6,N) → (N).
arma::vec batch_trace(const arma::mat &voigt, VoigtType vtype);

/// Infers the output VoigtType of a tensor4 @ tensor2 contraction from the tensor4 type.
VoigtType infer_contraction_vtype(Tensor4Type t4type);
/// Batch contract tensor4 @ tensor2. t4:(6,6,N4), t2:(6,N2) → (6,N); pair with infer_contraction_vtype for the output type.
arma::mat batch_contract(const arma::cube &t4, Tensor4Type t4type,
                         const arma::mat &t2, VoigtType t2_vtype);

/// Batch rotate N tensor4. t4:(6,6,N), rot:(3,3,N_r).
arma::cube batch_rotate_t4(const arma::cube &t4, Tensor4Type t4type,
                           const arma::cube &rot_matrices, bool active = true);

/// Batch push-forward N tensor4. t4:(6,6,N), F:(3,3,N_f).
arma::cube batch_push_forward_t4(const arma::cube &t4, Tensor4Type t4type,
                                 const arma::cube &F, bool metric = true);

/// Batch pull-back N tensor4. t4:(6,6,N), F:(3,3,N_f).
arma::cube batch_pull_back_t4(const arma::cube &t4, Tensor4Type t4type,
                              const arma::cube &F, bool metric = true);

/// Infers the resulting Tensor4Type of a tensor4 inverse from the input type.
Tensor4Type infer_inverse_type(Tensor4Type t4type);
/// Batch inverse N tensor4. t4:(6,6,N) → (6,6,N); pair with infer_inverse_type for the output type.
arma::cube batch_inverse_t4(const arma::cube &t4, Tensor4Type t4type);

} // namespace simcoon
