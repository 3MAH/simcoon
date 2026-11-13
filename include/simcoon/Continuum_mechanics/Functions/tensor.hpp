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
#include <stdexcept>
#include <simcoon/Simulation/Maths/small_voigt.hpp>
#include <simcoon/Simulation/Maths/small_dense.hpp>

/**
 * @file tensor.hpp
 * @brief Typed 2nd- and 4th-order tensor API: type-tagged Voigt/Mandel conversions,
 * convention-aware contractions and rotations (Voigt order [11,22,33,12,13,23], units MPa).
 */

namespace simcoon {

// Forward declarations
class Rotation;
class tensor4;

/**
 * @brief Type tag for 2nd-order tensors, determines Voigt conversion factors and rotation rules.
 *
 * - stress:  Voigt = [s11, s22, s33, s12, s13, s23] — shear as-is
 * - strain:  Voigt = [e11, e22, e33, 2*e12, 2*e13, 2*e23] — shear doubled
 * - symmetric: same storage as stress (symmetric tensor, no physical convention)
 * - none:      no Voigt convention: any 3x3 (e.g. F, R, PK1), .voigt() throws
 *
 * `generic` is the pre-2.2 name of `symmetric`, kept as an alias for one release.
 */
enum class Tensor2Type {
    stress,
    strain,
    symmetric,
    none,
    generic = symmetric   ///< deprecated alias of symmetric (pre-2.2 name)
};

/// Deprecated pre-2.0 alias — the tag was named after the Voigt encoding
/// rather than the tensor rank (asymmetric with Tensor4Type). Kept for
/// external users (e.g. fedoo); prefer Tensor2Type.
using VoigtType = Tensor2Type;

/**
 * @brief Shear weights of the engineering Voigt notation per Tensor2Type, for the stack
 * conversions mat_from_voigt / voigt_from_mat (small_voigt.hpp): a strain carries
 * \f$ \gamma_{ij} = 2\varepsilon_{ij} \f$ (0.5 on the way to the 3x3, a sum on the way back),
 * a stress or a symmetric tensor its true shears (1, then the average). The only place
 * where the type decides the factor.
 */
inline double voigt_shear_in(Tensor2Type t)  { return (t == Tensor2Type::strain) ? 0.5 : 1.0; }
/// See voigt_shear_in.
inline double voigt_shear_out(Tensor2Type t) { return (t == Tensor2Type::strain) ? 1.0 : 0.5; }
/**
 * @brief Shear weight of the contraction \f$ \mathbf{A}:\mathbf{B} \f$ of two Voigt 6-vectors of types @p a and
 * @p b, \f$ w = 2\,\kappa_a \kappa_b \f$ with \f$ \kappa \f$ = voigt_shear_in: 2 for two stresses, 1/2 for two
 * strains, 1 for a stress with a strain (see contract_voigt). Tensor2Type::none reads as a stress, as in
 * voigt_shear_in.
 */
inline double voigt_shear_weight(Tensor2Type a, Tensor2Type b) {
    return 2.*voigt_shear_in(a)*voigt_shear_in(b);
}
/// Self-contraction weight, voigt_shear_weight(t, t).
inline double voigt_shear_weight(Tensor2Type t) { return voigt_shear_weight(t, t); }

/// Coefficient of the Mises equivalent, \f$ \sigma_{eq} = \sqrt{\tfrac{3}{2}\,\mathbf{s}':\mathbf{s}'} \f$ for a stress,
/// \f$ \varepsilon_{eq} = \sqrt{\tfrac{2}{3}\,\mathbf{e}':\mathbf{e}'} \f$ for a strain.
inline double mises_factor(Tensor2Type t) { return (t == Tensor2Type::strain) ? 2./3. : 3./2.; }

/**
 * @brief Mises equivalent of the Voigt 6-vector @p v of type @p t: \f$ \sqrt{c\,\mathbf{d}:\mathbf{d}} \f$ with
 * \f$ \mathbf{d} \f$ its deviator (written to @p d), \f$ c \f$ = mises_factor and the contraction weight
 * voigt_shear_weight. The one kernel behind Mises_stress, Mises_strain, eta_stress, eta_strain and batch_mises.
 */
inline double mises_voigt(const double *v, Tensor2Type t, double *d) {
    dev_voigt(v, d);
    return std::sqrt(mises_factor(t)*contract_voigt(d, d, voigt_shear_weight(t)));
}

/// Flow \f$ c\,\mathbf{x}/n \f$ of the Voigt 6-vector @p x with its shear entries scaled by @p shear (2 for a
/// strain-like output of a stress-like @p x), zero at \f$ n \le 0 \f$: the kernel of eta_stress and of the
/// fused von Mises flow of the modular yield criterion.
inline arma::vec::fixed<6> flow_voigt(const double *x, const double n, const double c, const double shear) {
    if (n <= 0.) return arma::vec::fixed<6>(arma::fill::zeros);
    const double k = c/n;
    arma::vec::fixed<6> eta;
    for (int i=0; i<3; i++) eta[i] = k*x[i];
    for (int i=3; i<6; i++) eta[i] = shear*k*x[i];
    return eta;
}

/// tensor2::from_voigt(v, t).mat() on the stack, without the object; @p v must hold 6 values.
inline arma::mat::fixed<3,3> voigt_to_mat(const arma::vec &v, Tensor2Type t) {
    if (v.n_elem != 6) {
        throw std::invalid_argument("voigt_to_mat: expected 6-element vector, got " + std::to_string(v.n_elem));
    }
    return mat_from_voigt(v.memptr(), voigt_shear_in(t));
}

/// tensor2(m, t).voigt() into @p v (resized to 6), without the object.
inline void mat_to_voigt(const arma::mat::fixed<3,3> &m, Tensor2Type t, arma::vec &v) {
    v.set_size(6);
    voigt_from_mat(v.memptr(), m, voigt_shear_out(t));
}

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

/// Per-type engineering->Mandel shear exponents: the congruence multiplies shear rows (3..5) by
/// \f$ 2^{e_r/2} \f$ and shear cols by \f$ 2^{e_c/2} \f$, \f$ e = \pm 1 \f$ (Mandel->engineering divides).
/// Encodes the Tensor4Type table above; this is the ONLY place the \f$\sqrt2\f$ patterns live.
inline void mandel_shear_exponents(Tensor4Type type, int &er, int &ec) {
    switch (type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::generic:               er = +1; ec = +1; break;
        case Tensor4Type::compliance:            er = -1; ec = -1; break;
        case Tensor4Type::strain_concentration:  er = -1; ec = +1; break;
        case Tensor4Type::stress_concentration:  er = +1; ec = -1; break;
    }
}

/// The shear factors \f$ 2^{e_r/2}, 2^{e_c/2} \f$ of mandel_shear_exponents.
inline void mandel_factors(Tensor4Type type, double &row, double &col) {
    int er, ec; mandel_shear_exponents(type, er, ec);
    const double s2 = std::sqrt(2.0);
    row = (er > 0) ? s2 : 1.0/s2;
    col = (ec > 0) ? s2 : 1.0/s2;
}

/// The exact product of the two mandel_factors, \f$ 2^{(e_r+e_c)/2} \f$ (2 stiffness, 1/2 compliance, 1
/// concentration): the engineering shear diagonal of an isotropic operator is its deviatoric coefficient
/// divided by it.
inline double mandel_shear_product(Tensor4Type type) {
    int er, ec; mandel_shear_exponents(type, er, ec);
    return std::ldexp(1.0, (er + ec)/2);
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
 * @note The correction is the solver's own per-corate dispatch
 *       (Dtau_LieDD_2_DtauDe_corate in objective_rates.cpp), so a tangent built
 *       here with `metric = false` is exactly the (Kirchhoff) box tangent the solver
 *       integrates the stress against; the default `metric = true` divides it by
 *       \f$ J \f$ (Cauchy level):
 *       - `logarithmic` (corate 2) and `logarithmic_R` (corate 3) share the exact
 *         spectral map (Dtau_LieDD_Dtau_logarithmicDD), so they give the same
 *         tangent. Their stress updates differ only in the transport operator.
 *       - `logarithmic_F` (corate 5) is the Jaumann box chained through
 *         \f$ (\mathbf{A}^F)^{-1} \f$ (De = A^F:D dt), distinct from both.
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
 * Storage: arma::mat::fixed<3,3> (column-major, 72 bytes of payload) plus the type tag.
 * The Voigt vector is computed on the fly (no cache).
 *
 * @details **Basis.** The components are those in the fixed orthonormal lab basis
 * \f$ \mathbf{e}_i \f$; the class records no basis. Consequences:
 * - `rotate(R, active = true)` returns another tensor,
 *   \f$ \mathbf{Q}\,\mathbf{X}\,\mathbf{Q}^T \f$, in the same basis;
 *   `active = false` returns the components \f$ \mathbf{Q}^T\,\mathbf{X}\,\mathbf{Q} \f$ of the
 *   same tensor in the frame \f$ \mathbf{e}'_i = \mathbf{Q}\,\mathbf{e}_i \f$, and the result
 *   does not remember that frame.
 * - In `push_forward` / `pull_back` the type tag is the variance: stress is contravariant
 *   (\f$ \mathbf{F}\,\mathbf{X}\,\mathbf{F}^T \f$), strain covariant
 *   (\f$ \mathbf{F}^{-T}\,\mathbf{X}\,\mathbf{F}^{-1} \f$). This reading holds for the conjugate
 *   pairs \f$ (\mathbf{S}, \mathbf{E}) \f$ and \f$ (\boldsymbol{\tau}, \mathbf{e}) \f$
 *   (Green-Lagrange / Almansi), not for a logarithmic strain.
 * - A push-forward is a transport (another tensor, in the current configuration) that carries
 *   the basis along: the components of \f$ \boldsymbol{\tau} \f$ in the convected basis
 *   \f$ \mathbf{g}_i = \mathbf{F}\,\mathbf{G}_i \f$ are those of \f$ \mathbf{S} \f$ in
 *   \f$ \mathbf{G}_i \f$, \f$ \tau^{ij} = S^{IJ} \f$; `push_forward` returns them in the lab.
 * - The `metric` flag of the transports is the Piola weight \f$ J = \det\mathbf{F} \f$
 *   (Kirchhoff to Cauchy), not a metric tensor.
 * - The variance is a convention of the representation, not a property of the tensor. In C++
 *   it is read from the type tag and there is no index raising or lowering (the Python classes
 *   carry it as a tag, `variance`, and change it with `to_variance`). Two-point tensors
 *   (\f$ \mathbf{F} \f$, \f$ \mathbf{R} \f$, \f$ \Delta\mathbf{R} \f$, \f$ \mathbf{P} \f$) have
 *   one index in each configuration and no single basis: they are held as lab-lab components
 *   (type `none`) or as plain arma::mat, never transported by these functions.
 *
 * The Python classes `simcoon.Tensor2` / `simcoon.Tensor4` can carry a basis
 * (`simcoon.Basis`: orthonormal, or natural with its metric); these C++ value types do not.
 */
class tensor2 {
private:
    arma::mat::fixed<3,3> _mat;
    Tensor2Type _vtype;

public:
    /// Zero tensor, Tensor2Type::stress.
    tensor2();
    /// Zero tensor with the given type tag.
    explicit tensor2(Tensor2Type vtype);
    /// Wrap a 3x3 matrix (true components, no Voigt factor) with a type tag.
    tensor2(const arma::mat::fixed<3,3> &m, Tensor2Type vtype);
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2(const arma::mat &m, Tensor2Type vtype);

    /// Build from an engineering Voigt 6-vector \f$[t_{11},t_{22},t_{33},t_{12},t_{13},t_{23}]\f$;
    /// strain halves the shear entries (engineering \f$\gamma_{ij}=2\varepsilon_{ij}\f$), stress uses them as-is.
    static tensor2 from_voigt(const arma::vec::fixed<6> &v, Tensor2Type vtype);
    /// Same, from a dynamic arma::vec (must have 6 elements).
    static tensor2 from_voigt(const arma::vec &v, Tensor2Type vtype);
    /// Same, with a string tag ("stress"/"strain"/"symmetric"/"none") — parses per call; prefer the enum overload in loops.
    static tensor2 from_voigt(const arma::vec &v, const std::string &type_str);

    /// Build from a Kelvin-Mandel vector (shear scaled by \f$\sqrt2\f$, identical for stress/strain).
    static tensor2 from_mandel(const arma::vec::fixed<6> &v, Tensor2Type vtype);

    /// Zero tensor factory.
    static tensor2 zeros(Tensor2Type vtype = Tensor2Type::stress);
    /// Zero tensor factory (string tag).
    static tensor2 zeros(const std::string &type_str);
    /// Identity tensor \f$\mathbf{I}\f$ factory.
    static tensor2 identity(Tensor2Type vtype = Tensor2Type::stress);
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
    Tensor2Type vtype() const { return _vtype; }

    /**
     * @brief Compute Voigt vector on the fly (6 lookups + factor-2 on shear for strain).
     * @return vec::fixed<6> by value
     * @throws std::runtime_error if Tensor2Type::none
     */
    arma::vec::fixed<6> voigt() const;

    /**
     * @brief Kelvin-Mandel 6-vector: \f$[t_{11},t_{22},t_{33},\sqrt2\,t_{12},\sqrt2\,t_{13},\sqrt2\,t_{23}]\f$.
     *
     * Convention-symmetric (identical for stress and strain); the natural input to a
     * tensor4 Mandel contraction.
     */
    arma::vec::fixed<6> mandel() const;

    /// True if the stored 3x3 is symmetric within @p tol.
    bool is_symmetric(double tol = 1e-12) const;

    /// Copy out as a dynamic arma::mat.
    arma::mat to_arma_mat() const { return arma::mat(_mat); }
    /// Copy out the engineering Voigt vector as a dynamic arma::vec.
    arma::vec to_arma_voigt() const { return arma::vec(voigt()); }

    /// Rotate: \f$\mathbf{Q}\,\mathbf{X}\,\mathbf{Q}^T\f$ (passive uses \f$\mathbf{Q}^T\f$).
    tensor2 rotate(const Rotation &R, bool active = true) const;

    /// Push-forward to the spatial configuration; kernel dispatches on Tensor2Type
    /// (stress: \f$\mathbf{F}\mathbf{X}\mathbf{F}^T\f$, strain: \f$\mathbf{F}^{-T}\mathbf{X}\mathbf{F}^{-1}\f$).
    /// metric=true (default) includes the \f$J=\det\mathbf{F}\f$ Piola factor; metric=false = pure transport.
    tensor2 push_forward(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2 push_forward(const arma::mat &F, bool metric = true) const;
    /// Pull-back to the reference configuration (exact inverse of push_forward).
    tensor2 pull_back(const arma::mat::fixed<3,3> &F, bool metric = true) const;
    /// Same, from a dynamic arma::mat (must be 3x3).
    tensor2 pull_back(const arma::mat &F, bool metric = true) const;

    /// Componentwise sum. The type tags are not checked: the result takes the left one.
    tensor2 operator+(const tensor2 &other) const;
    /// Componentwise difference. The type tags are not checked: the result takes the left one.
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
 * @return \f$\partial\Phi/\partial\sigma\f$, Tensor2Type::strain
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
 * @return \f$\partial H/\partial\sigma\f$, Tensor2Type::strain
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
 * @brief A 4th-order tensor stored in the Kelvin-Mandel convention, with a type tag.
 *
 * Storage: arma::mat::fixed<6,6> in Mandel (sqrt2 on shear rows/cols), so matrix algebra
 * equals tensor algebra (identity = eye(6), inverse = arma::inv, rotation = orthogonal
 * R6 * X * R6^T). The engineering Voigt form is exposed at the boundary via mat()/set_mat();
 * the per-type engineering<->Mandel congruence is documented on Tensor4Type.
 *
 * Every transport is a Mandel congruence \f$ \mathbb{Q}\,\mathbb{C}\,\mathbb{Q}^T \f$ with
 * \f$ \mathbb{Q} \f$ the 6x6 Mandel matrix of \f$ \mathbf{X} \mapsto \mathbf{Q}\mathbf{X}\mathbf{Q}^T \f$
 * (see push_forward): the full-index products of the theory are never formed.
 *
 * @details **Basis.** As for tensor2, the components are those in the fixed orthonormal lab
 * basis. The Kelvin-Mandel shortcuts above (identity = eye(6), projectors, rotation as an
 * orthogonal congruence, Frobenius norm) hold only in an orthonormal basis: in a natural basis
 * of metric \f$ g_{ij} \f$ the stiffness-type identity is
 * \f$ \tfrac12 (g^{ik} g^{jl} + g^{il} g^{jk}) \f$. The Tensor4Type is the variance used by
 * `push_forward` / `pull_back`: stiffness contravariant on its four indices, compliance
 * covariant, concentration tensors mixed.
 *
 * Const methods do not mutate the object: a tensor4 can be read from several threads.
 */
class tensor4 {
private:
    arma::mat::fixed<6,6> _mandel;   ///< Kelvin-Mandel 6x6 (sqrt2 on shear); engineering form via mat()
    Tensor4Type _type;

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

    tensor4(const tensor4 &other) = default;
    tensor4(tensor4 &&other) noexcept = default;
    tensor4& operator=(const tensor4 &other) = default;
    tensor4& operator=(tensor4 &&other) noexcept = default;
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

    // Conversion
    arma::mat to_arma_mat() const { return arma::mat(mat()); }

    /**
     * @brief Contract with a tensor2 (double contraction), e.g. sigma = L : eps.
     *        Internally a plain Mandel matrix-vector product.
     *
     * Output Tensor2Type is inferred from Tensor4Type:
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
     *       \f$ C'^{ijkl} = F_{iA} F_{jB} F_{kC} F_{lD}\, C^{ABCD} \f$, scaled by \f$ 1/J \f$ if metric=true
     *   - compliance (covariant):
     *       \f$ M'_{ijkl} = F^{-T}_{iA} F^{-T}_{jB} F^{-T}_{kC} F^{-T}_{lD}\, M_{ABCD} \f$, scaled by \f$ J \f$
     *   - strain_concentration / stress_concentration: not implemented (mixed indices)
     *
     * Computed as the Mandel congruence \f$ \mathbb{Q}_{\mathbf{K}}\,\mathbb{C}\,\mathbb{Q}_{\mathbf{K}}^T \f$,
     * \f$ \mathbf{K} \f$ the kernel above (\f$ \mathbf{F} \f$ or \f$ \mathbf{F}^{-T} \f$) and
     * \f[
     *   (\mathbb{Q}_{\mathbf{K}})_{IJ} = \frac{c_I}{c_J}\left( K_{ik} K_{jl} + (1-\delta_{kl})\, K_{il} K_{jk} \right),
     *   \quad I = (ij),\ J = (kl),\ c = (1,1,1,\sqrt2,\sqrt2,\sqrt2),
     * \f]
     * the Mandel matrix of the linear map \f$ \mathbf{X} \mapsto \mathbf{K}\mathbf{X}\mathbf{K}^T \f$ on the
     * symmetric second-order tensors. The congruence is the index product above, exactly, for ANY
     * \f$ \mathbf{F} \f$ (neither symmetric nor orthogonal is assumed) and a tensor with both minor
     * symmetries. \f$ \mathbb{Q}_{\mathbf{K}} \f$ is not orthogonal in general:
     * \f$ \mathbb{Q}_{\mathbf{K}}^{-1} = \mathbb{Q}_{\mathbf{K}^{-1}} \neq \mathbb{Q}_{\mathbf{K}}^T \f$, which is
     * why pull_back() applies the inverse kernel and never the transpose. It is orthogonal exactly
     * when \f$ \mathbf{K} \f$ is a rotation: rotate() shares the formula, not the property. Mixed
     * variance (one kernel per pair) is not a congruence: the concentration types are refused.
     *
     * The compliance kernel is what makes inverse(K).push_forward(F) == inverse(K.push_forward(F))
     * for non-orthogonal F.
     *
     * @param F deformation gradient (3x3)
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
     * @param metric If true (default), the corrected Kirchhoff tangent is divided by
     *        \f$ J \f$; pass false for the box tangent the solver integrates
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
     * @param F deformation gradient (3x3)
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

    /// Componentwise sum. The type tags are not checked: the result takes the left one.
    tensor4 operator+(const tensor4 &other) const;
    /// Componentwise difference. The type tags are not checked: the result takes the left one.
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
/// @p n_threads caps the threads of the loop (0 = runtime default, 1 = serial; see parallel.hpp).
arma::mat batch_rotate(const arma::mat &voigt, Tensor2Type vtype,
                       const arma::cube &rot_matrices, bool active = true,
                       unsigned int n_threads = 0);

/// Batch push-forward N tensor2. voigt:(6,N), F:(3,3,N_f).
arma::mat batch_push_forward(const arma::mat &voigt, Tensor2Type vtype,
                             const arma::cube &F, bool metric = true);

/// Batch pull-back N tensor2. voigt:(6,N), F:(3,3,N_f).
arma::mat batch_pull_back(const arma::mat &voigt, Tensor2Type vtype,
                          const arma::cube &F, bool metric = true);

/// Batch von Mises for N tensor2. voigt:(6,N) → (N).
arma::vec batch_mises(const arma::mat &voigt, Tensor2Type vtype);

/// Batch trace for N tensor2. voigt:(6,N) → (N).
arma::vec batch_trace(const arma::mat &voigt, Tensor2Type vtype);

/// Infers the output Tensor2Type of a tensor4 @ tensor2 contraction from the tensor4 type.
Tensor2Type infer_contraction_vtype(Tensor4Type t4type);
/// Batch contract tensor4 @ tensor2. t4:(6,6,N4), t2:(6,N2) → (6,N); pair with infer_contraction_vtype for the output type.
arma::mat batch_contract(const arma::cube &t4, Tensor4Type t4type,
                         const arma::mat &t2, Tensor2Type t2_vtype);

/// Batch rotate N tensor4 (overload of the tensor2 version). t4:(6,6,N), rot:(3,3,N_r).
arma::cube batch_rotate(const arma::cube &t4, Tensor4Type t4type,
                        const arma::cube &rot_matrices, bool active = true,
                        unsigned int n_threads = 0);

/// Batch push-forward N tensor4 (overload of the tensor2 version). t4:(6,6,N), F:(3,3,N_f).
arma::cube batch_push_forward(const arma::cube &t4, Tensor4Type t4type,
                              const arma::cube &F, bool metric = true);

/// Batch pull-back N tensor4 (overload of the tensor2 version). t4:(6,6,N), F:(3,3,N_f).
arma::cube batch_pull_back(const arma::cube &t4, Tensor4Type t4type,
                           const arma::cube &F, bool metric = true);

/// Infers the resulting Tensor4Type of a tensor4 inverse from the input type.
Tensor4Type infer_inverse_type(Tensor4Type t4type);
/// Batch inverse N tensor4. t4:(6,6,N) → (6,6,N); pair with infer_inverse_type for the output type.
arma::cube batch_inverse(const arma::cube &t4, Tensor4Type t4type);

} // namespace simcoon
