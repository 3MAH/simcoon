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
#include <cmath>
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

/// Per-type engineering->Mandel shear factors: the congruence multiplies shear rows
/// (3..5) by @p row and shear cols by @p col (Mandel->engineering divides). Encodes the
/// Tensor4Type table above; this is the ONLY place the \f$\sqrt2\f$ patterns live.
inline void mandel_factors(Tensor4Type type, double &row, double &col) {
    const double s2 = std::sqrt(2.0);
    switch (type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::generic:               row = s2;     col = s2;     break;
        case Tensor4Type::compliance:            row = 1.0/s2; col = 1.0/s2; break;
        case Tensor4Type::strain_concentration:  row = 1.0/s2; col = s2;     break;
        case Tensor4Type::stress_concentration:  row = s2;     col = 1.0/s2; break;
    }
}

/// Engineering Voigt 6x6 -> Kelvin-Mandel 6x6, per-type congruence (see Tensor4Type).
inline arma::mat::fixed<6,6> eng_to_mandel(arma::mat::fixed<6,6> X, Tensor4Type type) {
    double r, c; mandel_factors(type, r, c);
    for (int I = 3; I < 6; ++I) X.row(I) *= r;
    for (int J = 3; J < 6; ++J) X.col(J) *= c;
    return X;
}

/// Kelvin-Mandel 6x6 -> engineering Voigt 6x6 (exact inverse of eng_to_mandel).
inline arma::mat::fixed<6,6> mandel_to_eng(arma::mat::fixed<6,6> X, Tensor4Type type) {
    double r, c; mandel_factors(type, r, c);
    for (int I = 3; I < 6; ++I) X.row(I) /= r;
    for (int J = 3; J < 6; ++J) X.col(J) /= c;
    return X;
}

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
    /// Zero tensor, VoigtType::stress.
    tensor2();
    /// Zero tensor with the given type tag.
    explicit tensor2(VoigtType vtype);
    /// Wrap a 3x3 matrix (true components, no Voigt factor) with a type tag.
    tensor2(const arma::mat::fixed<3,3> &m, VoigtType vtype);
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2(const arma::mat &m, VoigtType vtype);

    /// Build from an engineering Voigt 6-vector \f$[t_{11},t_{22},t_{33},t_{12},t_{13},t_{23}]\f$;
    /// strain halves the shear entries (engineering \f$\gamma_{ij}=2\varepsilon_{ij}\f$), stress uses them as-is.
    static tensor2 from_voigt(const arma::vec::fixed<6> &v, VoigtType vtype);
    /// Same, from a dynamic arma::vec (must have 6 elements).
    static tensor2 from_voigt(const arma::vec &v, VoigtType vtype);
    /// Same, with a string tag ("stress"/"strain"/"generic") — parses per call; prefer the enum overload in loops.
    static tensor2 from_voigt(const arma::vec &v, const std::string &type_str);

    /// Build from a Kelvin-Mandel vector (shear scaled by \f$\sqrt2\f$, identical for stress/strain).
    static tensor2 from_mandel(const arma::vec::fixed<6> &v, VoigtType vtype);

    /// Zero tensor factory.
    static tensor2 zeros(VoigtType vtype = VoigtType::stress);
    /// Zero tensor factory (string tag).
    static tensor2 zeros(const std::string &type_str);
    /// Identity tensor \f$\mathbf{I}\f$ factory.
    static tensor2 identity(VoigtType vtype = VoigtType::stress);
    /// Identity tensor factory (string tag).
    static tensor2 identity(const std::string &type_str);

    // Copy/move
    tensor2(const tensor2 &other) = default;
    tensor2(tensor2 &&other) noexcept = default;
    tensor2& operator=(const tensor2 &other) = default;
    tensor2& operator=(tensor2 &&other) noexcept = default;
    ~tensor2() = default;

    /// The stored 3x3 matrix (true components, no Voigt factor).
    const arma::mat::fixed<3,3>& mat() const { return _mat; }
    /// Mutable access to the stored 3x3.
    arma::mat::fixed<3,3>& mat_mut();
    /// Replace the stored 3x3.
    void set_mat(const arma::mat::fixed<3,3> &m);
    /// Replace the stored 3x3 (dynamic arma::mat, must be 3x3).
    void set_mat(const arma::mat &m);
    /// Replace components from an engineering Voigt 6-vector (respects the current type tag).
    void set_voigt(const arma::vec::fixed<6> &v);
    /// Same, from a dynamic arma::vec (must have 6 elements).
    void set_voigt(const arma::vec &v);

    /// The Voigt type tag.
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

    /// True if the stored 3x3 is symmetric within @p tol.
    bool is_symmetric(double tol = 1e-12) const;

    /// Copy out as a dynamic arma::mat.
    arma::mat to_arma_mat() const { return arma::mat(_mat); }
    /// Copy out the engineering Voigt vector as a dynamic arma::vec.
    arma::vec to_arma_voigt() const { return arma::vec(voigt()); }

    /// Rotate: \f$\mathbf{Q}\,\mathbf{X}\,\mathbf{Q}^T\f$ (passive uses \f$\mathbf{Q}^T\f$).
    tensor2 rotate(const Rotation &R, bool active = true) const;

    /// Push-forward to the spatial configuration; kernel dispatches on VoigtType
    /// (stress: \f$\mathbf{F}\mathbf{X}\mathbf{F}^T\f$, strain: \f$\mathbf{F}^{-T}\mathbf{X}\mathbf{F}^{-1}\f$).
    /// metric=true (default) includes the \f$J=\det\mathbf{F}\f$ Piola factor; metric=false = pure transport.
    tensor2 push_forward(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2 push_forward(const arma::mat &F, bool metric = true) const;
    /// Pull-back to the reference configuration (exact inverse of push_forward).
    tensor2 pull_back(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2 pull_back(const arma::mat &F, bool metric = true) const;

    /// Componentwise sum (same type tag required).
    tensor2 operator+(const tensor2 &other) const;
    /// Componentwise difference (same type tag required).
    tensor2 operator-(const tensor2 &other) const;
    /// Unary minus.
    tensor2 operator-() const;
    /// Scalar multiple.
    tensor2 operator*(double scalar) const;
    /// Scalar division.
    tensor2 operator/(double scalar) const;
    /// In-place sum.
    tensor2& operator+=(const tensor2 &other);
    /// In-place difference.
    tensor2& operator-=(const tensor2 &other);
    /// In-place scalar multiple.
    tensor2& operator*=(double scalar);
    /// In-place scalar division.
    tensor2& operator/=(double scalar);
    /// Scalar multiple (scalar on the left).
    friend tensor2 operator*(double scalar, const tensor2 &t);

    /// Element-wise (Schur) product on the 3x3 matrix, like Armadillo %.
    /// Use sum(a % b) for the double contraction scalar.
    tensor2 operator%(const tensor2 &other) const;

    /// Element-wise division on the 3x3 matrix.
    tensor2 operator/(const tensor2 &other) const;

    /// Exact componentwise equality (matrix and type tag).
    bool operator==(const tensor2 &other) const;
    /// Negation of operator==.
    bool operator!=(const tensor2 &other) const;
};

/// Stress-typed tensor2 from a 3x3 matrix.
tensor2 stress(const arma::mat::fixed<3,3> &m);
/// Stress-typed tensor2 from an engineering Voigt 6-vector.
tensor2 stress(const arma::vec::fixed<6> &v);
/// Stress-typed tensor2 from a dynamic 3x3 arma::mat.
tensor2 stress(const arma::mat &m);
/// Stress-typed tensor2 from a dynamic 6-element arma::vec.
tensor2 stress(const arma::vec &v);
/// Strain-typed tensor2 from a 3x3 matrix.
tensor2 strain(const arma::mat::fixed<3,3> &m);
/// Strain-typed tensor2 from an engineering Voigt 6-vector (shear = \f$\gamma_{ij}=2\varepsilon_{ij}\f$).
tensor2 strain(const arma::vec::fixed<6> &v);
/// Strain-typed tensor2 from a dynamic 3x3 arma::mat.
tensor2 strain(const arma::mat &m);
/// Strain-typed tensor2 from a dynamic 6-element arma::vec.
tensor2 strain(const arma::vec &v);
/// Deviatoric part \f$\mathbf{X} - \tfrac{1}{3}\mathrm{tr}(\mathbf{X})\,\mathbf{I}\f$ (keeps the type tag).
tensor2 dev(const tensor2 &t);
/// Von Mises equivalent: stress \f$\sqrt{\tfrac{3}{2}\,\mathbf{s}:\mathbf{s}}\f$, strain \f$\sqrt{\tfrac{2}{3}\,\mathbf{e}:\mathbf{e}}\f$ (dispatch on the type tag).
double Mises(const tensor2 &t);
/// Trace \f$\mathrm{tr}(\mathbf{X})\f$.
double trace(const tensor2 &t);

/**
 * @brief Normal to a von Mises / Drucker–Prager yield surface, as a plastic-flow direction.
 *
 * Strain-typed gradient of the yield surface
 * \f$ \Phi(\sigma) = \sigma_{eq} + \alpha\,\mathrm{tr}\,\sigma - k \f$:
 * \f[ \frac{\partial\Phi}{\partial\sigma}
 *      = \frac{3}{2}\,\frac{\mathrm{dev}\,\sigma}{\sigma_{eq}} + \alpha\,I . \f]
 * The default \f$\alpha=0\f$ gives the pure von Mises normal
 * \f$(3/2)\,\mathrm{dev}\,\sigma/\sigma_{eq}\f$ (identical to eta_stress()). When
 * \f$\sigma_{eq}\le\f$ @c iota the deviatoric normal is undefined (cone vertex) and
 * is dropped; only the volumetric term \f$\alpha I\f$ is kept.
 *
 * @param sigma stress tensor
 * @param alpha pressure-sensitivity (friction) coefficient; 0 → von Mises
 * @return \f$\partial\Phi/\partial\sigma\f$, VoigtType::strain
 * @see flow() for the plastic-potential normal (non-associative flow)
 */
tensor2 flow_normal(const tensor2 &sigma, double alpha = 0.0);

/**
 * @brief Normal to the plastic potential — the actual plastic-flow direction.
 *
 * Same form as flow_normal() but for the potential
 * \f$ H(\sigma) = \sigma_{eq} + \beta\,\mathrm{tr}\,\sigma \f$, so that
 * \f$ \dot\varepsilon^p = \dot\lambda\,\partial H/\partial\sigma \f$ with
 * \f[ \frac{\partial H}{\partial\sigma}
 *      = \frac{3}{2}\,\frac{\mathrm{dev}\,\sigma}{\sigma_{eq}} + \beta\,I . \f]
 * Flow is non-associative when \f$H\neq\Phi\f$ (\f$\beta\neq\alpha\f$). Non-negative
 * dissipation \f$\sigma:\partial H/\partial\sigma\ge 0\f$ (2nd law) is guaranteed when
 * \f$H\f$ is convex and dominates the yield surface, \f$H\ge\Phi\f$ — here
 * \f$0\le\beta\le\alpha\f$. \f$\beta=0\f$ recovers associative von Mises.
 *
 * @param sigma stress tensor
 * @param beta  dilatancy coefficient (\f$\le\alpha\f$ for the 2nd law)
 * @return \f$\partial H/\partial\sigma\f$, VoigtType::strain
 * @see flow_normal() for the yield-surface normal
 */
tensor2 flow(const tensor2 &sigma, double beta = 0.0);

// Reduction free functions (Armadillo-style)
/// Sum of all 9 components; sum(a % b) is the double contraction \f$\mathbf{a}:\mathbf{b}\f$.
double sum(const tensor2 &t);
/// Alias of sum() (Armadillo naming).
double accu(const tensor2 &t);
/// Frobenius norm \f$\sqrt{\mathbf{X}:\mathbf{X}}\f$.
double norm(const tensor2 &t);
/// Determinant of the 3x3.
double det(const tensor2 &t);
/// Componentwise absolute value (keeps the type tag).
tensor2 abs(const tensor2 &t);
/// Transpose (keeps the type tag).
tensor2 trans(const tensor2 &t);

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

public:
    /// Wrap an already-Mandel 6x6 directly (skips the engineering->Mandel congruence).
    static tensor4 from_mandel(const arma::mat::fixed<6,6> &m_mandel, Tensor4Type type);

    /// Wrap an ENGINEERING Voigt 6x6 (same as the matrix constructor, but the call site
    /// names the convention — explicit counterpart of from_mandel).
    static tensor4 from_voigt(const arma::mat::fixed<6,6> &m, Tensor4Type type);
    /// Same, from a dynamic arma::mat (must be 6x6).
    static tensor4 from_voigt(const arma::mat &m, Tensor4Type type);
    /// Same, with a string tag ("stiffness"/"compliance"/...).
    static tensor4 from_voigt(const arma::mat &m, const std::string &type_str);

    /// Zero tensor, Tensor4Type::stiffness.
    tensor4();
    /// Zero tensor with the given type tag.
    explicit tensor4(Tensor4Type type);
    /// Wrap an ENGINEERING Voigt 6x6 with a type tag (converted to Mandel internally).
    tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type);
    /// Same, from a dynamic arma::mat (must be 6x6).
    tensor4(const arma::mat &m, Tensor4Type type);
    /// Same, with a string tag ("stiffness"/"compliance"/...) — parses per call; prefer the enum overload in loops.
    tensor4(const arma::mat &m, const std::string &type_str);

    /// Identity projector. Stored in Mandel, every type shares identity = eye(6); the
    /// engineering form via mat() recovers the classic convention (Ireal for stiffness,
    /// Ireal2 for compliance, eye(6) for the concentration types).
    static tensor4 identity(Tensor4Type type = Tensor4Type::stiffness);
    /// Volumetric (spherical) projector \f$\tfrac{1}{3}\,\mathbf{I}\otimes\mathbf{I}\f$.
    static tensor4 volumetric(Tensor4Type type = Tensor4Type::stiffness);
    /// Deviatoric projector: identity() - volumetric().
    static tensor4 deviatoric(Tensor4Type type = Tensor4Type::stiffness);
    /// Zero tensor factory.
    static tensor4 zeros(Tensor4Type type = Tensor4Type::stiffness);

    // Copy/move
    tensor4(const tensor4 &other);
    tensor4(tensor4 &&other) noexcept;
    tensor4& operator=(const tensor4 &other);
    tensor4& operator=(tensor4 &&other) noexcept;
    ~tensor4() = default;

    /// The ENGINEERING Voigt 6x6, by value (converted from the internal Mandel storage).
    /// mat_mut() is intentionally absent: there is no coherent mutable engineering view
    /// over Mandel storage — use set_mat().
    arma::mat::fixed<6,6> mat() const;
    /// Replace components from an engineering Voigt 6x6.
    void set_mat(const arma::mat::fixed<6,6> &m);
    /// Same, from a dynamic arma::mat (must be 6x6).
    void set_mat(const arma::mat &m);

    /// The Kelvin-Mandel 6x6 — the internal storage, by reference (no conversion).
    /// Convention-symmetric across types: identity = eye(6), tensor inverse = arma::inv,
    /// rotation = orthogonal congruence, for every Tensor4Type.
    const arma::mat::fixed<6,6>& mandel() const { return _mandel; }

    /// The tensor4 type tag.
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

    /// Rotate all four indices: exact Mandel congruence \f$R_6\,\mathbf{X}\,R_6^T\f$,
    /// valid for every Tensor4Type (passive uses \f$\mathbf{Q}^T\f$).
    tensor4 rotate(const Rotation &R, bool active = true) const;

    /// Componentwise sum (same type tag required).
    tensor4 operator+(const tensor4 &other) const;
    /// Componentwise difference (same type tag required).
    tensor4 operator-(const tensor4 &other) const;
    /// Unary minus.
    tensor4 operator-() const;
    /// Double contraction, alias of contract(): \f$\boldsymbol{\sigma} = \mathbf{L}:\boldsymbol{\varepsilon}\f$.
    tensor2 operator*(const tensor2 &t) const;
    /// Scalar multiple.
    tensor4 operator*(double scalar) const;
    /// Scalar division.
    tensor4 operator/(double scalar) const;
    /// In-place sum.
    tensor4& operator+=(const tensor4 &other);
    /// In-place difference.
    tensor4& operator-=(const tensor4 &other);
    /// In-place scalar multiple.
    tensor4& operator*=(double scalar);
    /// In-place scalar division.
    tensor4& operator/=(double scalar);
    /// Scalar multiple (scalar on the left).
    friend tensor4 operator*(double scalar, const tensor4 &t);

    /// Element-wise (Schur) product on the 6x6 ENGINEERING Voigt matrix (convention-bound,
    /// so it is computed on mat(), not on the Mandel storage).
    tensor4 operator%(const tensor4 &other) const;

    /// Exact componentwise equality (Mandel matrix and type tag).
    bool operator==(const tensor4 &other) const;
    /// Negation of operator==.
    bool operator!=(const tensor4 &other) const;
};

/// Dyadic product \f$C_{ijkl} = a_{ij}\,b_{kl}\f$, returned as a stiffness-typed tensor4.
tensor4 dyadic(const tensor2 &a, const tensor2 &b);
/// Dyadic product of a tensor with itself: \f$C_{ijkl} = a_{ij}\,a_{kl}\f$.
tensor4 auto_dyadic(const tensor2 &a);
/// Dyadic product of the symmetric parts, as the outer product of symmetrized Voigt
/// vectors: \f$\mathbf{C} = v(\mathbf{a})\,v(\mathbf{b})^T\f$ (see contimech sym_dyadic).
tensor4 sym_dyadic(const tensor2 &a, const tensor2 &b);
/// sym_dyadic of a tensor with itself.
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
