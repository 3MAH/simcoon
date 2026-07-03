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

///@file tensor.cpp

#include <iostream>
#include <stdexcept>
#include <exception>
#include <cmath>
#include <cstring>
#include <limits>
#include <armadillo>
#include <Fastor/Fastor.h>
#include <simcoon/parameter.hpp>
#include <simcoon/parallel.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

namespace simcoon {

// Helper: string → VoigtType
static VoigtType parse_voigt_type(const std::string &s) {
    if (s == "stress")  return VoigtType::stress;
    if (s == "strain")  return VoigtType::strain;
    if (s == "generic") return VoigtType::generic;
    if (s == "none")    return VoigtType::none;
    throw std::invalid_argument("Unknown VoigtType string: '" + s + "'. "
        "Expected: stress, strain, generic, none");
}

// Helper: string → Tensor4Type
static Tensor4Type parse_tensor4_type(const std::string &s) {
    if (s == "stiffness")             return Tensor4Type::stiffness;
    if (s == "compliance")            return Tensor4Type::compliance;
    if (s == "strain_concentration")  return Tensor4Type::strain_concentration;
    if (s == "stress_concentration")  return Tensor4Type::stress_concentration;
    if (s == "generic")               return Tensor4Type::generic;
    throw std::invalid_argument("Unknown Tensor4Type string: '" + s + "'. "
        "Expected: stiffness, compliance, strain_concentration, stress_concentration, generic");
}

// Helper: validate that a secondary operand in a batch op is either broadcast
// (n_slices == 1) or matches the primary batch size. Without this check, a
// mismatched count silently reads past the owned cube buffer in release builds
// (ARMA_NO_DEBUG), producing garbage output and, under OpenMP, non-deterministic
// crashes.
static void check_batch_broadcast(arma::uword n_secondary,
                                  arma::uword N_primary,
                                  const char *fn,
                                  const char *secondary_name,
                                  const char *primary_name) {
    if (n_secondary != 1 && n_secondary != N_primary) {
        throw std::invalid_argument(
            std::string(fn) + ": " + secondary_name + ".n_slices ("
            + std::to_string(n_secondary) + ") must be 1 or "
            + primary_name + " (" + std::to_string(N_primary) + ")");
    }
}

// Helper: compute inv of a 3x3 fixed matrix, returning fixed
static arma::mat::fixed<3,3> inv33(const arma::mat::fixed<3,3> &F) {
    arma::mat tmp;
    bool ok = arma::inv(tmp, arma::mat(F));
    if (!ok)
        throw std::runtime_error("tensor: cannot invert singular 3x3 matrix (det ~ 0)");
    arma::mat::fixed<3,3> result;
    result = tmp;
    return result;
}

// Helper: select the 3x3 contraction kernel for tensor4 push/pull on all 4 indices.
//   forward=true  (push):  stiffness/generic → F,       compliance → F^{-T}
//   forward=false (pull):  stiffness/generic → F^{-1},  compliance → F^T
// Concentration types mix covariant and contravariant indices — not implemented.
static arma::mat::fixed<3,3> tensor4_kernel(const arma::mat::fixed<3,3> &F,
                                            Tensor4Type type, bool forward) {
    switch (type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::generic:
            return forward ? F : inv33(F);
        case Tensor4Type::compliance:
            return forward ? arma::mat::fixed<3,3>(inv33(F).t()) : arma::mat::fixed<3,3>(F.t());
        case Tensor4Type::strain_concentration:
        case Tensor4Type::stress_concentration:
            throw std::runtime_error(
                std::string("tensor4::") + (forward ? "push_forward" : "pull_back")
                + " not implemented for concentration tensors "
                  "(mixed covariant/contravariant indices)");
    }
    return F; // unreachable
}

// ============================================================================
// Kelvin-Mandel convention
// ============================================================================
//
// tensor4 stores its 6x6 in Mandel. Boundaries: engineering <-> Mandel per-type
// congruence eng_to_mandel/mandel_to_eng (tensor.hpp); Mandel <-> full-index via the
// type-free mandel_to_fastor4/fastor4_to_mandel (fastor_bridge.hpp).

static const double SQ2 = std::sqrt(2.0);

// Mandel rotation operator: sigma_hat' = R6 * sigma_hat reproduces Q*sigma*Q^T for symmetric
// tensors. R6 is orthogonal, so one congruence R6 * X * R6^T rotates every Tensor4Type exactly
// (~3x cheaper than the full-index einsum route, no _fastor cache build).
// Entry: R6(I,J) = (c_I/c_J) * (Q_ik Q_jl + (1-delta_kl) Q_il Q_jk), c = (1,1,1,sqrt2,sqrt2,sqrt2).
static arma::mat::fixed<6,6> mandel_rotation(const arma::mat::fixed<3,3> &Q) {
    static const int pi[6] = {0,1,2,0,0,1};   // Voigt order [11,22,33,12,13,23]
    static const int pj[6] = {0,1,2,1,2,2};
    arma::mat::fixed<6,6> R6;
    for (int I = 0; I < 6; ++I) {
        const int i = pi[I], j = pj[I];
        const double cI = (I < 3) ? 1.0 : SQ2;
        for (int J = 0; J < 6; ++J) {
            const int k = pi[J], l = pj[J];
            double term = Q(i,k)*Q(j,l);
            if (k != l) term += Q(i,l)*Q(j,k);
            R6(I,J) = (J < 3) ? cI*term : (cI/SQ2)*term;
        }
    }
    return R6;
}

// Scalar 6x6 products for the Mandel congruences: arma Mat*Mat dispatches to BLAS dgemm,
// which at this size costs more than the flops and contends (Accelerate) when the batch
// loops run one call per item across threads.
static arma::mat::fixed<6,6> mul66(const arma::mat::fixed<6,6> &A, const arma::mat::fixed<6,6> &B) {
    arma::mat::fixed<6,6> C;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += A(i,k) * B(k,j);
            C(i,j) = s;
        }
    return C;
}

// A * B^T without forming the transpose.
static arma::mat::fixed<6,6> mul66_ABt(const arma::mat::fixed<6,6> &A, const arma::mat::fixed<6,6> &B) {
    arma::mat::fixed<6,6> C;
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int k = 0; k < 6; ++k) s += A(i,k) * B(j,k);
            C(i,j) = s;
        }
    return C;
}

// Volumetric projector (Mandel == engineering here: it has no shear entries).
static arma::mat::fixed<6,6> base_volumetric() {
    arma::mat::fixed<6,6> m;
    m.zeros();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) m(i,j) = 1.0/3.0;
    return m;
}

// ============================================================================
// tensor2 implementation
// ============================================================================

tensor2::tensor2() : _mat(arma::fill::zeros), _vtype(VoigtType::stress) {}

tensor2::tensor2(VoigtType vtype) : _mat(arma::fill::zeros), _vtype(vtype) {}

tensor2::tensor2(const arma::mat::fixed<3,3> &m, VoigtType vtype)
    : _mat(m), _vtype(vtype) {}

tensor2::tensor2(const arma::mat &m, VoigtType vtype)
    : _vtype(vtype) {
    if (m.n_rows != 3 || m.n_cols != 3)
        throw std::invalid_argument("tensor2: expected 3x3 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _mat = m;
}

tensor2 tensor2::from_voigt(const arma::vec::fixed<6> &v, VoigtType vtype) {
    arma::mat::fixed<3,3> m;
    if (vtype == VoigtType::strain) {
        m(0,0) = v(0); m(1,1) = v(1); m(2,2) = v(2);
        m(0,1) = 0.5 * v(3); m(1,0) = 0.5 * v(3);
        m(0,2) = 0.5 * v(4); m(2,0) = 0.5 * v(4);
        m(1,2) = 0.5 * v(5); m(2,1) = 0.5 * v(5);
    } else if (vtype == VoigtType::stress || vtype == VoigtType::generic) {
        m(0,0) = v(0); m(1,1) = v(1); m(2,2) = v(2);
        m(0,1) = v(3); m(1,0) = v(3);
        m(0,2) = v(4); m(2,0) = v(4);
        m(1,2) = v(5); m(2,1) = v(5);
    } else {
        throw std::runtime_error("Cannot construct tensor2 from Voigt with VoigtType::none");
    }
    tensor2 result(m, vtype);
    return result;
}

tensor2 tensor2::from_voigt(const arma::vec &v, VoigtType vtype) {
    if (v.n_elem != 6)
        throw std::invalid_argument("tensor2: expected 6-element vector, got "
            + std::to_string(v.n_elem));
    arma::vec::fixed<6> vf(v.memptr());
    return from_voigt(vf, vtype);
}

tensor2 tensor2::from_voigt(const arma::vec &v, const std::string &type_str) {
    return from_voigt(v, parse_voigt_type(type_str));
}

tensor2 tensor2::from_mandel(const arma::vec::fixed<6> &v, VoigtType vtype) {
    // Mandel shear = sqrt2 * true component, identical for stress and strain.
    arma::mat::fixed<3,3> m;
    m(0,0) = v(0); m(1,1) = v(1); m(2,2) = v(2);
    m(0,1) = m(1,0) = v(3) / SQ2;
    m(0,2) = m(2,0) = v(4) / SQ2;
    m(1,2) = m(2,1) = v(5) / SQ2;
    return tensor2(m, vtype);
}

tensor2 tensor2::zeros(VoigtType vtype) {
    return tensor2(vtype);
}

tensor2 tensor2::zeros(const std::string &type_str) {
    return zeros(parse_voigt_type(type_str));
}

tensor2 tensor2::identity(VoigtType vtype) {
    return tensor2(arma::mat::fixed<3,3>(arma::fill::eye), vtype);
}

tensor2 tensor2::identity(const std::string &type_str) {
    return identity(parse_voigt_type(type_str));
}

arma::mat::fixed<3,3>& tensor2::mat_mut() {

    return _mat;
}

void tensor2::set_mat(const arma::mat::fixed<3,3> &m) {
    _mat = m;

}

void tensor2::set_mat(const arma::mat &m) {
    if (m.n_rows != 3 || m.n_cols != 3)
        throw std::invalid_argument("tensor2: expected 3x3 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _mat = m;

}

void tensor2::set_voigt(const arma::vec::fixed<6> &v) {
    if (_vtype == VoigtType::strain) {
        _mat(0,0) = v(0); _mat(1,1) = v(1); _mat(2,2) = v(2);
        _mat(0,1) = 0.5 * v(3); _mat(1,0) = 0.5 * v(3);
        _mat(0,2) = 0.5 * v(4); _mat(2,0) = 0.5 * v(4);
        _mat(1,2) = 0.5 * v(5); _mat(2,1) = 0.5 * v(5);
    } else if (_vtype == VoigtType::stress || _vtype == VoigtType::generic) {
        _mat(0,0) = v(0); _mat(1,1) = v(1); _mat(2,2) = v(2);
        _mat(0,1) = v(3); _mat(1,0) = v(3);
        _mat(0,2) = v(4); _mat(2,0) = v(4);
        _mat(1,2) = v(5); _mat(2,1) = v(5);
    } else {
        throw std::runtime_error("Cannot set_voigt on VoigtType::none tensor");
    }
}

void tensor2::set_voigt(const arma::vec &v) {
    if (v.n_elem != 6)
        throw std::invalid_argument("tensor2::set_voigt: expected 6-element vector, got "
            + std::to_string(v.n_elem));
    arma::vec::fixed<6> vf(v.memptr());
    set_voigt(vf);
}

arma::vec::fixed<6> tensor2::voigt() const {
    if (_vtype == VoigtType::none) {
        throw std::runtime_error("Cannot compute Voigt vector for VoigtType::none (non-symmetric tensor)");
    }

    arma::vec::fixed<6> v;
    v(0) = _mat(0,0); v(1) = _mat(1,1); v(2) = _mat(2,2);

    if (_vtype == VoigtType::strain) {
        v(3) = _mat(0,1) + _mat(1,0);
        v(4) = _mat(0,2) + _mat(2,0);
        v(5) = _mat(1,2) + _mat(2,1);
    } else {
        v(3) = 0.5 * (_mat(0,1) + _mat(1,0));
        v(4) = 0.5 * (_mat(0,2) + _mat(2,0));
        v(5) = 0.5 * (_mat(1,2) + _mat(2,1));
    }
    return v;
}

arma::vec::fixed<6> tensor2::mandel() const {
    if (_vtype == VoigtType::none) {
        throw std::runtime_error("Cannot compute Mandel vector for VoigtType::none (non-symmetric tensor)");
    }
    // Mandel shear = sqrt2 * true component (= sqrt2 * symmetric part), same for stress/strain.
    arma::vec::fixed<6> v;
    v(0) = _mat(0,0); v(1) = _mat(1,1); v(2) = _mat(2,2);
    v(3) = SQ2 * 0.5 * (_mat(0,1) + _mat(1,0));
    v(4) = SQ2 * 0.5 * (_mat(0,2) + _mat(2,0));
    v(5) = SQ2 * 0.5 * (_mat(1,2) + _mat(2,1));
    return v;
}

Fastor::Tensor<double,3,3> tensor2::fastor() const {
    return arma_to_fastor2(_mat, _vtype != VoigtType::none);
}

bool tensor2::is_symmetric(double tol) const {
    return (std::abs(_mat(0,1) - _mat(1,0)) < tol) &&
           (std::abs(_mat(0,2) - _mat(2,0)) < tol) &&
           (std::abs(_mat(1,2) - _mat(2,1)) < tol);
}

tensor2 tensor2::rotate(const Rotation &R, bool active) const {
    // The 3x3 holds the true tensor, so a single full-index rotation Q*X*Q^T covers every
    // VoigtType (stress, strain, generic, non-symmetric). Q orthogonal => exact for all.
    arma::mat::fixed<3,3> Q = R.as_matrix();
    if (!active) Q = arma::mat::fixed<3,3>(Q.t());
    arma::mat::fixed<3,3> result = Q * _mat * Q.t();
    return tensor2(result, _vtype);
}

tensor2 tensor2::push_forward(const arma::mat::fixed<3,3> &F, bool metric) const {
    switch (_vtype) {
        case VoigtType::stress: {
            arma::mat::fixed<3,3> result;
            result = F * _mat * F.t();
            if (metric) {
                double J = arma::det(F);
                result *= (1.0 / J);
            }
            return tensor2(result, VoigtType::stress);
        }
        case VoigtType::strain: {
            arma::mat::fixed<3,3> invF = inv33(F);
            arma::mat::fixed<3,3> result;
            result = invF.t() * _mat * invF;
            // strain factor is 1 — no metric correction needed
            return tensor2(result, VoigtType::strain);
        }
        case VoigtType::generic:
        case VoigtType::none:
            throw std::runtime_error("push_forward requires VoigtType::stress or VoigtType::strain");
    }
    return *this;
}

tensor2 tensor2::pull_back(const arma::mat::fixed<3,3> &F, bool metric) const {
    switch (_vtype) {
        case VoigtType::stress: {
            arma::mat::fixed<3,3> invF = inv33(F);
            arma::mat::fixed<3,3> result;
            result = invF * _mat * invF.t();
            if (metric) {
                double J = arma::det(F);
                result *= J;
            }
            return tensor2(result, VoigtType::stress);
        }
        case VoigtType::strain: {
            arma::mat::fixed<3,3> result;
            result = F.t() * _mat * F;
            // strain factor is 1 — no metric correction needed
            return tensor2(result, VoigtType::strain);
        }
        case VoigtType::generic:
        case VoigtType::none:
            throw std::runtime_error("pull_back requires VoigtType::stress or VoigtType::strain");
    }
    return *this;
}

tensor2 tensor2::push_forward(const arma::mat &F, bool metric) const {
    if (F.n_rows != 3 || F.n_cols != 3)
        throw std::invalid_argument("push_forward: expected 3x3 matrix, got "
            + std::to_string(F.n_rows) + "x" + std::to_string(F.n_cols));
    return push_forward(arma::mat::fixed<3,3>(F), metric);
}

tensor2 tensor2::pull_back(const arma::mat &F, bool metric) const {
    if (F.n_rows != 3 || F.n_cols != 3)
        throw std::invalid_argument("pull_back: expected 3x3 matrix, got "
            + std::to_string(F.n_rows) + "x" + std::to_string(F.n_cols));
    return pull_back(arma::mat::fixed<3,3>(F), metric);
}

tensor2 tensor2::operator+(const tensor2 &other) const {
    arma::mat::fixed<3,3> result;
    result = _mat + other._mat;
    return tensor2(result, _vtype);
}

tensor2 tensor2::operator-(const tensor2 &other) const {
    arma::mat::fixed<3,3> result;
    result = _mat - other._mat;
    return tensor2(result, _vtype);
}

tensor2 tensor2::operator-() const {
    arma::mat::fixed<3,3> result;
    result = -_mat;
    return tensor2(result, _vtype);
}

tensor2 tensor2::operator*(double scalar) const {
    arma::mat::fixed<3,3> result;
    result = _mat * scalar;
    return tensor2(result, _vtype);
}

tensor2 tensor2::operator/(double scalar) const {
    if (scalar == 0.0)
        throw std::runtime_error("tensor2: division by zero scalar");
    arma::mat::fixed<3,3> result;
    result = _mat / scalar;
    return tensor2(result, _vtype);
}

tensor2& tensor2::operator+=(const tensor2 &other) {
    _mat += other._mat;

    return *this;
}

tensor2& tensor2::operator-=(const tensor2 &other) {
    _mat -= other._mat;

    return *this;
}

tensor2& tensor2::operator*=(double scalar) {
    _mat *= scalar;
    return *this;
}

tensor2& tensor2::operator/=(double scalar) {
    if (scalar == 0.0)
        throw std::runtime_error("tensor2: division by zero scalar");
    _mat /= scalar;
    return *this;
}

tensor2 operator*(double scalar, const tensor2 &t) {
    arma::mat::fixed<3,3> result;
    result = t._mat * scalar;
    return tensor2(result, t._vtype);
}

tensor2 tensor2::operator%(const tensor2 &other) const {
    arma::mat::fixed<3,3> result;
    result = _mat % other._mat;
    return tensor2(result, _vtype);
}

tensor2 tensor2::operator/(const tensor2 &other) const {
    arma::mat::fixed<3,3> result;
    result = _mat / other._mat;
    return tensor2(result, _vtype);
}

bool tensor2::operator==(const tensor2 &other) const {
    return _vtype == other._vtype && arma::approx_equal(_mat, other._mat, "absdiff", 1e-14);
}

bool tensor2::operator!=(const tensor2 &other) const {
    return !(*this == other);
}

// Free functions
tensor2 stress(const arma::mat::fixed<3,3> &m) {
    return tensor2(m, VoigtType::stress);
}

tensor2 stress(const arma::vec::fixed<6> &v) {
    return tensor2::from_voigt(v, VoigtType::stress);
}

tensor2 stress(const arma::mat &m) {
    return tensor2(m, VoigtType::stress);
}

tensor2 stress(const arma::vec &v) {
    return tensor2::from_voigt(v, VoigtType::stress);
}

tensor2 strain(const arma::mat::fixed<3,3> &m) {
    return tensor2(m, VoigtType::strain);
}

tensor2 strain(const arma::vec::fixed<6> &v) {
    return tensor2::from_voigt(v, VoigtType::strain);
}

tensor2 strain(const arma::mat &m) {
    return tensor2(m, VoigtType::strain);
}

tensor2 strain(const arma::vec &v) {
    return tensor2::from_voigt(v, VoigtType::strain);
}

tensor2 dev(const tensor2 &t) {
    // Deviator directly on the 3x3 (convention-free, like Mises below): no Voigt round-trip.
    arma::mat::fixed<3,3> d = t.mat();
    double tr3 = (d(0,0) + d(1,1) + d(2,2)) / 3.0;
    d(0,0) -= tr3; d(1,1) -= tr3; d(2,2) -= tr3;
    return tensor2(d, t.vtype());
}

double Mises(const tensor2 &t) {
    if (t.vtype() == VoigtType::none)
        throw std::runtime_error("Mises not defined for VoigtType::none");
    // From the 3x3 deviator (convention-free): stress -> sqrt(3/2 dev:dev),
    // strain -> sqrt(2/3 dev:dev). accu(dev % dev) is the double contraction (off-diagonals
    // counted twice via symmetry), so no engineering shear factors are needed.
    arma::mat::fixed<3,3> d = t.mat();
    double tr = (d(0,0) + d(1,1) + d(2,2)) / 3.0;
    d(0,0) -= tr; d(1,1) -= tr; d(2,2) -= tr;
    double dd = arma::accu(d % d);
    if (t.vtype() == VoigtType::strain)
        return std::sqrt((2.0/3.0) * dd);
    return std::sqrt(1.5 * dd);
}

// Gradient of a Drucker-Prager-type surface  s_eq + coeff*tr(sigma):  (3/2) dev(sigma)/s_eq + coeff*I.
// Deviator and s_eq are formed inline (no Voigt round-trip). s_eq <= iota is the cone vertex:
// the deviatoric normal is undefined there, so it is dropped and only the volumetric part remains.
static tensor2 dp_surface_normal(const tensor2 &sigma, double coeff) {
    arma::mat::fixed<3,3> n = sigma.mat();
    double tr3 = (n(0,0) + n(1,1) + n(2,2)) / 3.0;
    n(0,0) -= tr3; n(1,1) -= tr3; n(2,2) -= tr3;          // s = dev(sigma)
    double seq = std::sqrt(1.5 * arma::accu(n % n));      // s_eq = sqrt(3/2 s:s)
    if (seq > simcoon::iota) n *= (1.5 / seq);            // (3/2) s / s_eq
    else                     n.zeros();
    n(0,0) += coeff; n(1,1) += coeff; n(2,2) += coeff;
    return tensor2(n, VoigtType::strain);
}

tensor2 flow_normal(const tensor2 &sigma, double alpha) {
    return dp_surface_normal(sigma, alpha);
}

tensor2 flow(const tensor2 &sigma, double beta) {
    return dp_surface_normal(sigma, beta);
}

double trace(const tensor2 &t) {
    return t.mat()(0,0) + t.mat()(1,1) + t.mat()(2,2);
}

double sum(const tensor2 &t) {
    return arma::accu(t.mat());
}

double accu(const tensor2 &t) {
    return arma::accu(t.mat());
}

double norm(const tensor2 &t) {
    return arma::norm(t.mat(), "fro");
}

double det(const tensor2 &t) {
    return arma::det(t.mat());
}

tensor2 abs(const tensor2 &t) {
    return tensor2(arma::mat::fixed<3,3>(arma::abs(t.mat())), t.vtype());
}

tensor2 trans(const tensor2 &t) {
    return tensor2(arma::mat::fixed<3,3>(t.mat().t()), t.vtype());
}

// ============================================================================
// tensor4 implementation
// ============================================================================

tensor4::tensor4() : _mandel(arma::fill::zeros), _type(Tensor4Type::stiffness) {}

tensor4::tensor4(Tensor4Type type) : _mandel(arma::fill::zeros), _type(type) {}

// Public constructors take the ENGINEERING Voigt form and store it in Mandel.
tensor4::tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type)
    : _mandel(eng_to_mandel(m, type)), _type(type) {}

tensor4::tensor4(const arma::mat &m, Tensor4Type type)
    : _type(type) {
    if (m.n_rows != 6 || m.n_cols != 6)
        throw std::invalid_argument("tensor4: expected 6x6 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _mandel = eng_to_mandel(arma::mat::fixed<6,6>(m), type);
}

tensor4::tensor4(const arma::mat &m, const std::string &type_str)
    : tensor4(m, parse_tensor4_type(type_str)) {}

// Store an already-Mandel 6x6 directly, no congruence.
tensor4 tensor4::from_mandel(const arma::mat::fixed<6,6> &m_mandel, Tensor4Type type) {
    tensor4 t;
    t._mandel = m_mandel;
    t._type = type;
    return t;
}

tensor4 tensor4::from_voigt(const arma::mat::fixed<6,6> &m, Tensor4Type type) {
    return tensor4(m, type);
}

tensor4 tensor4::from_voigt(const arma::mat &m, Tensor4Type type) {
    return tensor4(m, type);
}

tensor4 tensor4::from_voigt(const arma::mat &m, const std::string &type_str) {
    return tensor4(m, type_str);
}

tensor4 tensor4::identity(Tensor4Type type) {
    // Mandel identity is eye(6) for every type; mat() reads it back as Ireal/Ireal2/eye(6).
    return tensor4::from_mandel(arma::mat::fixed<6,6>(arma::fill::eye), type);
}

tensor4 tensor4::volumetric(Tensor4Type type) {
    return tensor4::from_mandel(base_volumetric(), type);
}

tensor4 tensor4::deviatoric(Tensor4Type type) {
    arma::mat::fixed<6,6> Idev(arma::fill::eye);
    Idev -= base_volumetric();
    return tensor4::from_mandel(Idev, type);
}

tensor4 tensor4::zeros(Tensor4Type type) {
    return tensor4(type);
}

arma::mat::fixed<6,6> tensor4::mat() const {
    return mandel_to_eng(_mandel, _type);
}

void tensor4::set_mat(const arma::mat::fixed<6,6> &m) {
    _mandel = eng_to_mandel(m, _type);
    _invalidate_fastor();
}

void tensor4::set_mat(const arma::mat &m) {
    if (m.n_rows != 6 || m.n_cols != 6)
        throw std::invalid_argument("tensor4: expected 6x6 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _mandel = eng_to_mandel(arma::mat::fixed<6,6>(m), _type);
    _invalidate_fastor();
}

void tensor4::_ensure_fastor() const {
    if (!_fastor) {
        _fastor.emplace(mandel_to_fastor4(_mandel));
    }
}

const Fastor::Tensor<double,3,3,3,3>& tensor4::fastor() const {
    _ensure_fastor();
    return *_fastor;
}

tensor2 tensor4::contract(const tensor2 &t) const {
    // Plain Mandel matrix-vector product; both operands and result are convention-symmetric.
    arma::vec::fixed<6> v_out = _mandel * t.mandel();
    return tensor2::from_mandel(v_out, infer_contraction_vtype(_type));
}

tensor4 tensor4::push_forward(const arma::mat::fixed<3,3> &F, bool metric) const {
    _ensure_fastor();

    // Kernel branches on _type (mirrors tensor2::push_forward).
    arma::mat::fixed<3,3> kernel_mat = tensor4_kernel(F, _type, /*forward=*/true);
    auto kernel_fastor = arma_to_fastor2(kernel_mat, false);
    Fastor::Tensor<double,3,3,3,3> result = push_forward_4(*_fastor, kernel_fastor);
    arma::mat::fixed<6,6> mandel_result = fastor4_to_mandel(result);

    if (metric) {
        double J = arma::det(F);
        // Piola scaling: stiffness→1/J, compliance→J (scalar, commutes with the congruence).
        // Concentration is rejected in tensor4_kernel, so only the two cases apply here.
        mandel_result *= (_type == Tensor4Type::compliance) ? J : (1.0 / J);
    }
    return tensor4::from_mandel(mandel_result, _type);
}

tensor4 tensor4::pull_back(const arma::mat::fixed<3,3> &F, bool metric) const {
    _ensure_fastor();

    // Pull-back inverts push-forward.
    arma::mat::fixed<3,3> kernel_mat = tensor4_kernel(F, _type, /*forward=*/false);
    auto kernel_fastor = arma_to_fastor2(kernel_mat, false);
    Fastor::Tensor<double,3,3,3,3> result = push_forward_4(*_fastor, kernel_fastor);
    arma::mat::fixed<6,6> mandel_result = fastor4_to_mandel(result);

    if (metric) {
        double J = arma::det(F);
        mandel_result *= (_type == Tensor4Type::compliance) ? (1.0 / J) : J;
    }
    return tensor4::from_mandel(mandel_result, _type);
}

tensor4 tensor4::push_forward(const arma::mat &F, bool metric) const {
    if (F.n_rows != 3 || F.n_cols != 3)
        throw std::invalid_argument("push_forward: expected 3x3 matrix, got "
            + std::to_string(F.n_rows) + "x" + std::to_string(F.n_cols));
    return push_forward(arma::mat::fixed<3,3>(F), metric);
}

tensor4 tensor4::pull_back(const arma::mat &F, bool metric) const {
    if (F.n_rows != 3 || F.n_cols != 3)
        throw std::invalid_argument("pull_back: expected 3x3 matrix, got "
            + std::to_string(F.n_rows) + "x" + std::to_string(F.n_cols));
    return pull_back(arma::mat::fixed<3,3>(F), metric);
}

tensor4 tensor4::push_forward(const arma::mat::fixed<3,3> &F, CoRate rate,
                               const tensor2 &tau, bool metric) const {
    // Lie rate = plain push-forward, no correction needed
    if (rate == CoRate::lie)
        return push_forward(F, metric);

    // The B-correction recipe is defined on a Truesdell-rate Kirchhoff tangent
    // (stress-stress, contravariant on all four indices). It does not apply to
    // compliance or concentration tangents.
    if (_type != Tensor4Type::stiffness && _type != Tensor4Type::generic) {
        throw std::runtime_error(
            "tensor4::push_forward(F, CoRate, tau): corotational-rate "
            "correction is only defined for stiffness/generic tangents");
    }

    // Step 1: Truesdell push-forward (F-transport) WITHOUT metric, giving the
    // Lie-rate Kirchhoff tangent. B-correction is applied at the Kirchhoff
    // level before metric scaling.
    auto F_fastor = arma_to_fastor2(F, false);
    _ensure_fastor();
    Fastor::Tensor<double,3,3,3,3> lie_full = push_forward_4(*_fastor, F_fastor);
    // The Dtau_* corrections below operate on the ENGINEERING Voigt (solver convention).
    arma::mat Lt_v = fastor4_to_voigt(lie_full, _type);
    const arma::mat::fixed<3,3> tau_mat = tau.mat();
    arma::mat result_v;

    // Step 2: Apply the corotational rate correction. The kernel is rate-specific and
    // MUST match the solver's per-corate dispatch (DSDE_2_DtauDe_corate / DtauDe_corate_2_DSDE
    // in objective_rates.cpp), otherwise a tangent built through this API disagrees with the
    // tangent the solver integrates the stress against (different objective B^(4) kernel).
    switch (rate) {
        case CoRate::jaumann:
            result_v = Dtau_LieDD_Dtau_JaumannDD(Lt_v, tau_mat);
            break;
        case CoRate::green_naghdi:
            result_v = Dtau_LieDD_Dtau_objectiveDD(Lt_v, get_BBBB_GN(F), tau_mat);
            break;
        case CoRate::logarithmic:    // XBM logarithmic spin -> get_BBBB kernel (solver corate 2)
            result_v = Dtau_LieDD_Dtau_objectiveDD(Lt_v, get_BBBB(F), tau_mat);
            break;
        case CoRate::logarithmic_R:  // R-transport == Green-Naghdi frame -> get_BBBB_GN kernel
            // Matches solver corate 3 (DtauDe_GreenNaghdiDD_2_DSDE): log_R transports by R,
            // so it shares the Green-Naghdi spin kernel, NOT the XBM one.
            result_v = Dtau_LieDD_Dtau_objectiveDD(Lt_v, get_BBBB_GN(F), tau_mat);
            break;
        case CoRate::logarithmic_F:  // F-transport == convected/Oldroyd (B = I) -> pure Lie
            // Matches solver corate 5 (Dtau_LieDD_2_DSDE): the Lie-rate tangent already computed
            // in Step 1 needs no spin correction.
            result_v = Lt_v;
            break;
        case CoRate::lie:
            // unreachable — handled at function entry
            break;
    }

    // Step 3: Apply metric factor (Kirchhoff → Cauchy) for stiffness/generic.
    if (metric) {
        double J = arma::det(F);
        result_v *= (1.0 / J);
    }

    return tensor4(arma::mat::fixed<6,6>(result_v), _type);
}

tensor4 tensor4::push_forward(const arma::mat &F, CoRate rate,
                               const tensor2 &tau, bool metric) const {
    if (F.n_rows != 3 || F.n_cols != 3)
        throw std::invalid_argument("push_forward: expected 3x3 matrix, got "
            + std::to_string(F.n_rows) + "x" + std::to_string(F.n_cols));
    return push_forward(arma::mat::fixed<3,3>(F), rate, tau, metric);
}

tensor4 tensor4::rotate(const Rotation &R, bool active) const {
    // One orthogonal rotation applied to all four indices, valid for every type: exact Mandel
    // congruence R6 * X * R6^T (see mandel_rotation) -- no per-type Voigt kernel, no full-index
    // einsum, no _fastor cache build.
    arma::mat::fixed<3,3> Q = R.as_matrix();
    if (!active) Q = arma::mat::fixed<3,3>(Q.t());
    const arma::mat::fixed<6,6> R6 = mandel_rotation(Q);
    return tensor4::from_mandel(mul66_ABt(mul66(R6, _mandel), R6), _type);
}

tensor4 tensor4::inverse() const {
    // In Mandel the tensor inverse is a plain matrix inverse (no convention bookkeeping).
    // Hand-rolled Gauss-Jordan with partial pivoting on the fixed 6x6: arma::inv on a
    // dynamic mat goes through heap + LAPACK per call, which serializes the parallel
    // batch loops (and costs more than the ~400 flops even in serial).
    arma::mat::fixed<6,6> A = _mandel;
    arma::mat::fixed<6,6> Inv(arma::fill::eye);
    for (int col = 0; col < 6; ++col) {
        int piv = col;
        double amax = std::fabs(A(col,col));
        for (int r = col + 1; r < 6; ++r) {
            const double v = std::fabs(A(r,col));
            if (v > amax) { amax = v; piv = r; }
        }
        if (amax < std::numeric_limits<double>::min())
            throw std::runtime_error("tensor4::inverse(): matrix is singular");
        if (piv != col) { A.swap_rows(piv, col); Inv.swap_rows(piv, col); }
        const double d = 1.0 / A(col,col);
        for (int j = 0; j < 6; ++j) { A(col,j) *= d; Inv(col,j) *= d; }
        for (int r = 0; r < 6; ++r) {
            if (r == col) continue;
            const double f = A(r,col);
            if (f == 0.0) continue;
            for (int j = 0; j < 6; ++j) { A(r,j) -= f * A(col,j); Inv(r,j) -= f * Inv(col,j); }
        }
    }
    return tensor4::from_mandel(Inv, infer_inverse_type(_type));
}

tensor4 tensor4::operator+(const tensor4 &other) const {
    // Linear ops commute with the per-type congruence (same type both sides) → stay in Mandel.
    return tensor4::from_mandel(_mandel + other._mandel, _type);
}

tensor4 tensor4::operator-(const tensor4 &other) const {
    return tensor4::from_mandel(_mandel - other._mandel, _type);
}

tensor4 tensor4::operator-() const {
    return tensor4::from_mandel(-_mandel, _type);
}

tensor2 tensor4::operator*(const tensor2 &t) const {
    return contract(t);
}

tensor4 tensor4::operator*(double scalar) const {
    return tensor4::from_mandel(_mandel * scalar, _type);
}

tensor4 tensor4::operator/(double scalar) const {
    if (scalar == 0.0)
        throw std::runtime_error("tensor4: division by zero scalar");
    return tensor4::from_mandel(_mandel / scalar, _type);
}

tensor4& tensor4::operator+=(const tensor4 &other) {
    _mandel += other._mandel;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator-=(const tensor4 &other) {
    _mandel -= other._mandel;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator*=(double scalar) {
    _mandel *= scalar;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator/=(double scalar) {
    if (scalar == 0.0)
        throw std::runtime_error("tensor4: division by zero scalar");
    _mandel /= scalar;
    _invalidate_fastor();
    return *this;
}

tensor4 operator*(double scalar, const tensor4 &t) {
    return tensor4::from_mandel(t._mandel * scalar, t._type);
}

tensor4 tensor4::operator%(const tensor4 &other) const {
    // Element-wise (Schur) product is convention-dependent: do it on the ENGINEERING form so
    // the result matches the documented "element-wise on the 6x6 Voigt matrix".
    return tensor4(arma::mat::fixed<6,6>(this->mat() % other.mat()), _type);
}

bool tensor4::operator==(const tensor4 &other) const {
    return _type == other._type && arma::approx_equal(_mandel, other._mandel, "absdiff", 1e-14);
}

bool tensor4::operator!=(const tensor4 &other) const {
    return !(*this == other);
}

// Free functions — dyadic product family
// Delegates to authoritative implementations in contimech.cpp

tensor4 sym_dyadic(const tensor2 &a, const tensor2 &b) {
    arma::mat C = simcoon::sym_dyadic(arma::mat(a.mat()), arma::mat(b.mat()));
    return tensor4(arma::mat::fixed<6,6>(C), Tensor4Type::stiffness);
}

tensor4 auto_sym_dyadic(const tensor2 &a) {
    arma::mat C = simcoon::auto_sym_dyadic(arma::mat(a.mat()));
    return tensor4(arma::mat::fixed<6,6>(C), Tensor4Type::stiffness);
}

tensor4 dyadic(const tensor2 &a, const tensor2 &b) {
    // Full outer product: C_ijkl = a_ij * b_kl (not restricted to symmetric Voigt)
    arma::mat C = simcoon::dyadic(arma::mat(a.mat()), arma::mat(b.mat()));
    return tensor4(arma::mat::fixed<6,6>(C), Tensor4Type::stiffness);
}

tensor4 auto_dyadic(const tensor2 &a) {
    return dyadic(a, a);
}

// ============================================================================
// Batch operations
// ============================================================================

VoigtType infer_contraction_vtype(Tensor4Type t4type) {
    switch (t4type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::stress_concentration:
        case Tensor4Type::generic:
            return VoigtType::stress;
        case Tensor4Type::compliance:
        case Tensor4Type::strain_concentration:
            return VoigtType::strain;
    }
    return VoigtType::stress;
}

Tensor4Type infer_inverse_type(Tensor4Type t4type) {
    if (t4type == Tensor4Type::stiffness) return Tensor4Type::compliance;
    if (t4type == Tensor4Type::compliance) return Tensor4Type::stiffness;
    return t4type;
}

// Batch loops below use simcoon_parallel_for_safe (parallel.hpp): exception-safe GCD/OpenMP,
// so a singular slice raises a catchable error instead of terminating the process.
//
// IMPORTANT: never call Cube::slice(i) inside these loops. Armadillo lazily builds a
// per-slice Mat accessor behind a per-cube std::mutex + heap allocation (Cube_meat.hpp,
// create_mat_ptr): on a fresh cube that is one lock+malloc per item, which serializes the
// whole parallel loop (profiled: threads spent ~10x more time in __psynch_mutexwait than
// computing). slice_memptr() is a plain pointer lookup — copy through these helpers instead.

template<arma::uword R, arma::uword C>
static arma::mat::fixed<R,C> slice_fixed(const arma::cube &c, arma::uword s) {
    arma::mat::fixed<R,C> m;
    std::memcpy(m.memptr(), c.slice_memptr(s), R*C*sizeof(double));
    return m;
}

static arma::vec::fixed<6> col_fixed6(const arma::mat &m, arma::uword j) {
    arma::vec::fixed<6> v;
    std::memcpy(v.memptr(), m.colptr(j), 6*sizeof(double));
    return v;
}

arma::mat batch_rotate(const arma::mat &voigt, VoigtType vtype,
                       const arma::cube &rot_matrices, bool active) {
    int N = voigt.n_cols;
    check_batch_broadcast(rot_matrices.n_slices, N,
                          "batch_rotate", "rot_matrices", "voigt.n_cols");
    bool broadcast = (rot_matrices.n_slices == 1);
    arma::mat result(6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor2 t = tensor2::from_voigt(col_fixed6(voigt, i), vtype);
        Rotation R = Rotation::from_matrix(slice_fixed<3,3>(rot_matrices, broadcast ? 0 : i));
        result.col(i) = t.rotate(R, active).voigt();
    });
    return result;
}

arma::mat batch_push_forward(const arma::mat &voigt, VoigtType vtype,
                             const arma::cube &F, bool metric) {
    int N = voigt.n_cols;
    check_batch_broadcast(F.n_slices, N,
                          "batch_push_forward", "F", "voigt.n_cols");
    bool broadcast = (F.n_slices == 1);
    arma::mat result(6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor2 t = tensor2::from_voigt(col_fixed6(voigt, i), vtype);
        arma::mat::fixed<3,3> Fi = slice_fixed<3,3>(F, broadcast ? 0 : i);
        result.col(i) = t.push_forward(Fi, metric).voigt();
    });
    return result;
}

arma::mat batch_pull_back(const arma::mat &voigt, VoigtType vtype,
                          const arma::cube &F, bool metric) {
    int N = voigt.n_cols;
    check_batch_broadcast(F.n_slices, N,
                          "batch_pull_back", "F", "voigt.n_cols");
    bool broadcast = (F.n_slices == 1);
    arma::mat result(6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor2 t = tensor2::from_voigt(col_fixed6(voigt, i), vtype);
        arma::mat::fixed<3,3> Fi = slice_fixed<3,3>(F, broadcast ? 0 : i);
        result.col(i) = t.pull_back(Fi, metric).voigt();
    });
    return result;
}

arma::vec batch_mises(const arma::mat &voigt, VoigtType vtype) {
    // Stress: sqrt(3/2 * (diag^2 + 2*shear^2))
    // Strain: sqrt(2/3 * (diag^2 + 0.5*shear^2))  (shear Voigt = 2*eij)
    arma::mat dev = voigt;
    arma::rowvec hyd = (voigt.row(0) + voigt.row(1) + voigt.row(2)) / 3.0;
    dev.row(0) -= hyd;
    dev.row(1) -= hyd;
    dev.row(2) -= hyd;

    arma::rowvec diag2 = arma::sum(dev.rows(0,2) % dev.rows(0,2), 0);
    arma::rowvec shear2 = arma::sum(dev.rows(3,5) % dev.rows(3,5), 0);

    if (vtype == VoigtType::strain) {
        return arma::sqrt((2.0/3.0) * (diag2 + 0.5 * shear2)).t();
    }
    return arma::sqrt(1.5 * (diag2 + 2.0 * shear2)).t();
}

arma::vec batch_trace(const arma::mat &voigt, VoigtType /*vtype*/) {
    // trace = v(0) + v(1) + v(2) regardless of VoigtType
    return arma::vectorise(voigt.row(0) + voigt.row(1) + voigt.row(2));
}

arma::mat batch_contract(const arma::cube &t4, Tensor4Type t4type,
                         const arma::mat &t2, VoigtType t2_vtype) {
    int N4 = t4.n_slices;
    int N2 = t2.n_cols;
    int N = std::max(N4, N2);
    // Both operands must be 1 (broadcast) or N (matched). Without these guards
    // the general-path loop reads past the owned buffer when the two non-
    // broadcast sizes disagree.
    check_batch_broadcast(N4, N, "batch_contract", "t4", "max(t4.n_slices,t2.n_cols)");
    check_batch_broadcast(N2, N, "batch_contract", "t2", "max(t4.n_slices,t2.n_cols)");
    bool bc4 = (N4 == 1);
    bool bc2 = (N2 == 1);

    // Fast path: single L, N strain vectors → one BLAS dgemm
    if (bc4) {
        arma::mat::fixed<6,6> L = slice_fixed<6,6>(t4, 0);
        if (bc2) {
            // 1 tensor4, 1 tensor2 → single mat-vec
            return L * t2;
        }
        // 1 tensor4, N tensor2 → single mat-mat multiply
        return L * t2;
    }

    // General path: N tensor4, 1 or N tensor2
    arma::mat result(6, N);
    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        result.col(i) = slice_fixed<6,6>(t4, i) * t2.col(bc2 ? 0 : i);
    }
    return result;
}

arma::cube batch_rotate_t4(const arma::cube &t4, Tensor4Type t4type,
                           const arma::cube &rot_matrices, bool active) {
    int N = t4.n_slices;
    check_batch_broadcast(rot_matrices.n_slices, N,
                          "batch_rotate_t4", "rot_matrices", "t4.n_slices");
    bool broadcast = (rot_matrices.n_slices == 1);
    arma::cube result(6, 6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor4 L(slice_fixed<6,6>(t4, i), t4type);
        Rotation R = Rotation::from_matrix(slice_fixed<3,3>(rot_matrices, broadcast ? 0 : i));
        const arma::mat::fixed<6,6> out = L.rotate(R, active).mat();
        std::memcpy(result.slice_memptr(i), out.memptr(), 36*sizeof(double));
    });
    return result;
}

// Concentration tensors mix covariant/contravariant indices, so push_forward/pull_back are
// undefined for them (also rejected per-element in tensor4_kernel). The batch loops must reject
// BEFORE their OpenMP parallel region: an exception escaping an active parallel region is UB.
static void reject_concentration_transport(Tensor4Type type, const char *op) {
    if (type == Tensor4Type::strain_concentration || type == Tensor4Type::stress_concentration)
        throw std::runtime_error(std::string("tensor4::") + op + " not implemented for "
                                 "concentration tensors (mixed covariant/contravariant indices)");
}

arma::cube batch_push_forward_t4(const arma::cube &t4, Tensor4Type t4type,
                                 const arma::cube &F, bool metric) {
    int N = t4.n_slices;
    check_batch_broadcast(F.n_slices, N,
                          "batch_push_forward_t4", "F", "t4.n_slices");
    reject_concentration_transport(t4type, "push_forward");
    bool broadcast = (F.n_slices == 1);
    arma::cube result(6, 6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor4 L(slice_fixed<6,6>(t4, i), t4type);
        arma::mat::fixed<3,3> Fi = slice_fixed<3,3>(F, broadcast ? 0 : i);
        const arma::mat::fixed<6,6> out = L.push_forward(Fi, metric).mat();
        std::memcpy(result.slice_memptr(i), out.memptr(), 36*sizeof(double));
    });
    return result;
}

arma::cube batch_pull_back_t4(const arma::cube &t4, Tensor4Type t4type,
                              const arma::cube &F, bool metric) {
    int N = t4.n_slices;
    check_batch_broadcast(F.n_slices, N,
                          "batch_pull_back_t4", "F", "t4.n_slices");
    reject_concentration_transport(t4type, "pull_back");
    bool broadcast = (F.n_slices == 1);
    arma::cube result(6, 6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor4 L(slice_fixed<6,6>(t4, i), t4type);
        arma::mat::fixed<3,3> Fi = slice_fixed<3,3>(F, broadcast ? 0 : i);
        const arma::mat::fixed<6,6> out = L.pull_back(Fi, metric).mat();
        std::memcpy(result.slice_memptr(i), out.memptr(), 36*sizeof(double));
    });
    return result;
}

arma::cube batch_inverse_t4(const arma::cube &t4, Tensor4Type t4type) {
    int N = t4.n_slices;
    arma::cube result(6, 6, N);

    simcoon_parallel_for_safe(N, [&](int i) {
        tensor4 L(slice_fixed<6,6>(t4, i), t4type);
        const arma::mat::fixed<6,6> out = L.inverse().mat();
        std::memcpy(result.slice_memptr(i), out.memptr(), 36*sizeof(double));
    });
    return result;
}

} // namespace simcoon
