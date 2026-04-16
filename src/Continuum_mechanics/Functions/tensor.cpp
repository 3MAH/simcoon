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
///@brief Implementation of tensor2 and tensor4 classes with type tags and Fastor integration.
///@version 1.0

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <armadillo>
#include <Fastor/Fastor.h>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Functions/fastor_bridge.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>

namespace simcoon {

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

// ============================================================================
// Voigt factor helpers for type-aware Voigt <-> full-index conversion
// ============================================================================
//
// The Voigt 6x6 matrix V_IJ relates to the full-index C_ijkl by:
//   V_IJ = voigt_factor(I, J, type) * C_ij(I)kl(J)
//
// Factors depend on whether input/output indices are stress or strain type:
//   Left factor:  2 if I>=3 AND output is strain (compliance, strain_concentration)
//   Right factor: 2 if J>=3 AND input is stress (compliance, stress_concentration)
//
//   Type                  | I<3,J<3 | I<3,J>=3 | I>=3,J<3 | I>=3,J>=3
//   stiffness             |    1    |    1     |    1     |    1
//   compliance            |    1    |    2     |    2     |    4
//   stress_concentration  |    1    |    2     |    1     |    2
//   strain_concentration  |    1    |    1     |    2     |    2

static double voigt_factor(int I, int J, Tensor4Type type) {
    double f = 1.0;
    if (I >= 3) {
        if (type == Tensor4Type::compliance || type == Tensor4Type::strain_concentration)
            f *= 2.0;
    }
    if (J >= 3) {
        if (type == Tensor4Type::compliance || type == Tensor4Type::stress_concentration)
            f *= 2.0;
    }
    return f;
}

// Convert Voigt matrix to full-index Fastor tensor, applying type-dependent factors
static Fastor::Tensor<double,3,3,3,3> voigt_to_full(
    const arma::mat::fixed<6,6> &V, Tensor4Type type)
{
    if (type == Tensor4Type::stiffness || type == Tensor4Type::generic) {
        return voigt_to_fastor4(V);
    }

    Fastor::Tensor<double,3,3,3,3> C;
    C.zeros();
    for (int i = 0; i < 3; ++i)
    for (int j = i; j < 3; ++j) {
        int I = voigt_map[i][j];
        for (int k = 0; k < 3; ++k)
        for (int l = k; l < 3; ++l) {
            int J = voigt_map[k][l];
            double val = V(I, J) / voigt_factor(I, J, type);
            C(i,j,k,l) = val;
            C(i,j,l,k) = val;
            C(j,i,k,l) = val;
            C(j,i,l,k) = val;
        }
    }
    return C;
}

// Convert full-index Fastor tensor to Voigt matrix, applying type-dependent factors
static arma::mat::fixed<6,6> full_to_voigt(
    const Fastor::Tensor<double,3,3,3,3> &C, Tensor4Type type)
{
    arma::mat::fixed<6,6> V = fastor4_to_voigt(C);

    if (type == Tensor4Type::stiffness || type == Tensor4Type::generic) {
        return V;
    }

    for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J)
            V(I, J) *= voigt_factor(I, J, type);
    return V;
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

tensor2 tensor2::zeros(VoigtType vtype) {
    return tensor2(vtype);
}

tensor2 tensor2::identity(VoigtType vtype) {
    return tensor2(arma::mat::fixed<3,3>(arma::fill::eye), vtype);
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

Fastor::Tensor<double,3,3> tensor2::fastor() const {
    return arma_to_fastor2(_mat, _vtype != VoigtType::none);
}

bool tensor2::is_symmetric(double tol) const {
    return (std::abs(_mat(0,1) - _mat(1,0)) < tol) &&
           (std::abs(_mat(0,2) - _mat(2,0)) < tol) &&
           (std::abs(_mat(1,2) - _mat(2,1)) < tol);
}

tensor2 tensor2::rotate(const Rotation &R, bool active) const {
    switch (_vtype) {
        case VoigtType::stress: {
            arma::vec::fixed<6> v_rot = R.apply_stress(voigt(), active);
            return tensor2::from_voigt(v_rot, VoigtType::stress);
        }
        case VoigtType::strain: {
            arma::vec::fixed<6> v_rot = R.apply_strain(voigt(), active);
            return tensor2::from_voigt(v_rot, VoigtType::strain);
        }
        case VoigtType::generic:
        case VoigtType::none: {
            arma::mat::fixed<3,3> R_mat = R.as_matrix();
            arma::mat::fixed<3,3> result;
            if (active) {
                result = R_mat * _mat * R_mat.t();
            } else {
                result = R_mat.t() * _mat * R_mat;
            }
            return tensor2(result, _vtype);
        }
    }
    return *this;
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
    arma::vec dv = simcoon::dev(arma::vec(t.voigt()));
    return tensor2::from_voigt(arma::vec::fixed<6>(dv.memptr()), t.vtype());
}

double Mises(const tensor2 &t) {
    arma::vec v(t.voigt());
    if (t.vtype() == VoigtType::stress || t.vtype() == VoigtType::generic)
        return Mises_stress(v);
    if (t.vtype() == VoigtType::strain)
        return Mises_strain(v);
    throw std::runtime_error("Mises not defined for VoigtType::none");
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

tensor4::tensor4() : _voigt(arma::fill::zeros), _type(Tensor4Type::stiffness) {}

tensor4::tensor4(Tensor4Type type) : _voigt(arma::fill::zeros), _type(type) {}

tensor4::tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type)
    : _voigt(m), _type(type) {}

tensor4::tensor4(const arma::mat &m, Tensor4Type type)
    : _type(type) {
    if (m.n_rows != 6 || m.n_cols != 6)
        throw std::invalid_argument("tensor4: expected 6x6 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _voigt = m;
}

tensor4::tensor4(const tensor4 &other)
    : _voigt(other._voigt), _type(other._type),
      _fastor(other._fastor) {}

tensor4::tensor4(tensor4 &&other) noexcept
    : _voigt(std::move(other._voigt)), _type(other._type),
      _fastor(std::move(other._fastor)) {}

tensor4& tensor4::operator=(const tensor4 &other) {
    if (this != &other) {
        _voigt = other._voigt;
        _type = other._type;
        _fastor = other._fastor;
    }
    return *this;
}

tensor4& tensor4::operator=(tensor4 &&other) noexcept {
    if (this != &other) {
        _voigt = std::move(other._voigt);
        _type = other._type;
        _fastor = std::move(other._fastor);
    }
    return *this;
}

tensor4 tensor4::identity(Tensor4Type type) {
    arma::mat::fixed<6,6> m = Ireal();
    return tensor4(m, type);
}

tensor4 tensor4::volumetric(Tensor4Type type) {
    arma::mat::fixed<6,6> m = Ivol();
    return tensor4(m, type);
}

tensor4 tensor4::deviatoric(Tensor4Type type) {
    arma::mat::fixed<6,6> m = Idev();
    return tensor4(m, type);
}

tensor4 tensor4::identity2(Tensor4Type type) {
    arma::mat::fixed<6,6> m = Ireal2();
    return tensor4(m, type);
}

tensor4 tensor4::deviatoric2(Tensor4Type type) {
    arma::mat::fixed<6,6> m = Idev2();
    return tensor4(m, type);
}

tensor4 tensor4::zeros(Tensor4Type type) {
    return tensor4(type);
}

arma::mat::fixed<6,6>& tensor4::mat_mut() {
    _invalidate_fastor();
    return _voigt;
}

void tensor4::set_mat(const arma::mat::fixed<6,6> &m) {
    _voigt = m;
    _invalidate_fastor();
}

void tensor4::set_mat(const arma::mat &m) {
    if (m.n_rows != 6 || m.n_cols != 6)
        throw std::invalid_argument("tensor4: expected 6x6 matrix, got "
            + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));
    _voigt = m;
    _invalidate_fastor();
}

void tensor4::_ensure_fastor() const {
    if (!_fastor) {
        _fastor.emplace(voigt_to_full(_voigt, _type));
    }
}

const Fastor::Tensor<double,3,3,3,3>& tensor4::fastor() const {
    _ensure_fastor();
    return *_fastor;
}

tensor2 tensor4::contract(const tensor2 &t) const {
    arma::vec::fixed<6> v_out = _voigt * t.voigt();

    VoigtType out_vtype;
    switch (_type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::stress_concentration:
            out_vtype = VoigtType::stress;
            break;
        case Tensor4Type::compliance:
        case Tensor4Type::strain_concentration:
            out_vtype = VoigtType::strain;
            break;
        case Tensor4Type::generic:
        default:
            out_vtype = VoigtType::stress;
            break;
    }

    return tensor2::from_voigt(v_out, out_vtype);
}

tensor4 tensor4::push_forward(const arma::mat::fixed<3,3> &F, bool metric) const {
    _ensure_fastor();

    auto F_fastor = arma_to_fastor2(F, false);  // F is non-symmetric

    Fastor::Tensor<double,3,3,3,3> result = push_forward_4(*_fastor, F_fastor);
    arma::mat::fixed<6,6> voigt_result = full_to_voigt(result, _type);

    if (metric) {
        double J = arma::det(F);
        switch (_type) {
            case Tensor4Type::stiffness:
            case Tensor4Type::generic:
                voigt_result *= (1.0 / J);
                break;
            case Tensor4Type::compliance:
                voigt_result *= J;
                break;
            case Tensor4Type::strain_concentration:
            case Tensor4Type::stress_concentration:
                break; // factor is 1
        }
    }
    return tensor4(voigt_result, _type);
}

tensor4 tensor4::pull_back(const arma::mat::fixed<3,3> &F, bool metric) const {
    _ensure_fastor();

    arma::mat::fixed<3,3> invF = inv33(F);

    auto invF_fastor = arma_to_fastor2(invF, false);  // invF is non-symmetric

    Fastor::Tensor<double,3,3,3,3> result = push_forward_4(*_fastor, invF_fastor);
    arma::mat::fixed<6,6> voigt_result = full_to_voigt(result, _type);

    if (metric) {
        double J = arma::det(F);
        switch (_type) {
            case Tensor4Type::stiffness:
            case Tensor4Type::generic:
                voigt_result *= J;
                break;
            case Tensor4Type::compliance:
                voigt_result *= (1.0 / J);
                break;
            case Tensor4Type::strain_concentration:
            case Tensor4Type::stress_concentration:
                break; // factor is 1
        }
    }
    return tensor4(voigt_result, _type);
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

    // Step 1: Lie push-forward WITHOUT metric (Kirchhoff-level tangent).
    // B-correction must be applied at the Kirchhoff level before metric scaling.
    arma::mat Lt_v = push_forward(F, /*metric=*/false).mat();
    const arma::mat::fixed<3,3> &tau_mat = tau.mat();
    arma::mat result_v;

    // Step 2: Apply corotational rate correction
    switch (rate) {
        case CoRate::jaumann:
            result_v = Dtau_LieDD_Dtau_JaumannDD(Lt_v, tau_mat);
            break;
        case CoRate::green_naghdi:
            result_v = Dtau_LieDD_Dtau_objectiveDD(Lt_v, get_BBBB_GN(F), tau_mat);
            break;
        case CoRate::logarithmic:
        case CoRate::logarithmic_R:
        case CoRate::logarithmic_F:
            result_v = Dtau_LieDD_Dtau_objectiveDD(Lt_v, get_BBBB(F), tau_mat);
            break;
        default:
            throw std::runtime_error("Unknown CoRate");
    }

    // Step 3: Apply metric factor (Kirchhoff → Cauchy)
    if (metric) {
        double J = arma::det(F);
        switch (_type) {
            case Tensor4Type::stiffness:
            case Tensor4Type::generic:
                result_v *= (1.0 / J);
                break;
            case Tensor4Type::compliance:
                result_v *= J;
                break;
            case Tensor4Type::strain_concentration:
            case Tensor4Type::stress_concentration:
                break;
        }
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
    switch (_type) {
        case Tensor4Type::stiffness:
        case Tensor4Type::generic: {
            arma::mat::fixed<6,6> result = R.apply_stiffness(_voigt, active);
            return tensor4(result, _type);
        }
        case Tensor4Type::compliance: {
            arma::mat::fixed<6,6> result = R.apply_compliance(_voigt, active);
            return tensor4(result, _type);
        }
        case Tensor4Type::strain_concentration: {
            arma::mat::fixed<6,6> result = R.apply_strain_concentration(_voigt, active);
            return tensor4(result, _type);
        }
        case Tensor4Type::stress_concentration: {
            arma::mat::fixed<6,6> result = R.apply_stress_concentration(_voigt, active);
            return tensor4(result, _type);
        }
    }
    return *this;
}

tensor4 tensor4::inverse() const {
    arma::mat inv_mat;
    bool ok = arma::inv(inv_mat, arma::mat(_voigt));
    if (!ok)
        throw std::runtime_error("tensor4::inverse(): Voigt matrix is singular");
    arma::mat::fixed<6,6> inv_voigt;
    inv_voigt = inv_mat;

    Tensor4Type inv_type;
    switch (_type) {
        case Tensor4Type::stiffness:
            inv_type = Tensor4Type::compliance;
            break;
        case Tensor4Type::compliance:
            inv_type = Tensor4Type::stiffness;
            break;
        default:
            inv_type = _type;
            break;
    }
    return tensor4(inv_voigt, inv_type);
}

tensor4 tensor4::operator+(const tensor4 &other) const {
    arma::mat::fixed<6,6> result;
    result = _voigt + other._voigt;
    return tensor4(result, _type);
}

tensor4 tensor4::operator-(const tensor4 &other) const {
    arma::mat::fixed<6,6> result;
    result = _voigt - other._voigt;
    return tensor4(result, _type);
}

tensor4 tensor4::operator-() const {
    arma::mat::fixed<6,6> result;
    result = -_voigt;
    return tensor4(result, _type);
}

tensor2 tensor4::operator*(const tensor2 &t) const {
    return contract(t);
}

tensor4 tensor4::operator*(double scalar) const {
    arma::mat::fixed<6,6> result;
    result = _voigt * scalar;
    return tensor4(result, _type);
}

tensor4 tensor4::operator/(double scalar) const {
    if (scalar == 0.0)
        throw std::runtime_error("tensor4: division by zero scalar");
    arma::mat::fixed<6,6> result;
    result = _voigt / scalar;
    return tensor4(result, _type);
}

tensor4& tensor4::operator+=(const tensor4 &other) {
    _voigt += other._voigt;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator-=(const tensor4 &other) {
    _voigt -= other._voigt;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator*=(double scalar) {
    _voigt *= scalar;
    _invalidate_fastor();
    return *this;
}

tensor4& tensor4::operator/=(double scalar) {
    if (scalar == 0.0)
        throw std::runtime_error("tensor4: division by zero scalar");
    _voigt /= scalar;
    _invalidate_fastor();
    return *this;
}

tensor4 operator*(double scalar, const tensor4 &t) {
    arma::mat::fixed<6,6> result;
    result = t._voigt * scalar;
    return tensor4(result, t._type);
}

tensor4 tensor4::operator%(const tensor4 &other) const {
    arma::mat::fixed<6,6> result;
    result = _voigt % other._voigt;
    return tensor4(result, _type);
}

bool tensor4::operator==(const tensor4 &other) const {
    return _type == other._type && arma::approx_equal(_voigt, other._voigt, "absdiff", 1e-14);
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

arma::mat batch_rotate(const arma::mat &voigt, VoigtType vtype,
                       const arma::cube &rot_matrices, bool active) {
    int N = voigt.n_cols;
    bool broadcast = (rot_matrices.n_slices == 1);
    arma::mat result(6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor2 t = tensor2::from_voigt(arma::vec(voigt.col(i)), vtype);
        Rotation R = Rotation::from_matrix(
            arma::mat::fixed<3,3>(rot_matrices.slice(broadcast ? 0 : i)));
        tensor2 t_rot = t.rotate(R, active);
        result.col(i) = t_rot.voigt();
    }
    return result;
}

arma::mat batch_push_forward(const arma::mat &voigt, VoigtType vtype,
                             const arma::cube &F, bool metric) {
    int N = voigt.n_cols;
    bool broadcast = (F.n_slices == 1);
    arma::mat result(6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor2 t = tensor2::from_voigt(arma::vec(voigt.col(i)), vtype);
        arma::mat::fixed<3,3> Fi(F.slice(broadcast ? 0 : i));
        result.col(i) = t.push_forward(Fi, metric).voigt();
    }
    return result;
}

arma::mat batch_pull_back(const arma::mat &voigt, VoigtType vtype,
                          const arma::cube &F, bool metric) {
    int N = voigt.n_cols;
    bool broadcast = (F.n_slices == 1);
    arma::mat result(6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor2 t = tensor2::from_voigt(arma::vec(voigt.col(i)), vtype);
        arma::mat::fixed<3,3> Fi(F.slice(broadcast ? 0 : i));
        result.col(i) = t.pull_back(Fi, metric).voigt();
    }
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
    bool bc4 = (N4 == 1);
    bool bc2 = (N2 == 1);

    // Fast path: single L, N strain vectors → one BLAS dgemm
    if (bc4) {
        arma::mat::fixed<6,6> L(t4.slice(0));
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
        result.col(i) = arma::mat::fixed<6,6>(t4.slice(i)) * t2.col(bc2 ? 0 : i);
    }
    return result;
}

arma::cube batch_rotate_t4(const arma::cube &t4, Tensor4Type t4type,
                           const arma::cube &rot_matrices, bool active) {
    int N = t4.n_slices;
    bool broadcast = (rot_matrices.n_slices == 1);
    arma::cube result(6, 6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor4 L(arma::mat::fixed<6,6>(t4.slice(i)), t4type);
        Rotation R = Rotation::from_matrix(
            arma::mat::fixed<3,3>(rot_matrices.slice(broadcast ? 0 : i)));
        result.slice(i) = L.rotate(R, active).mat();
    }
    return result;
}

arma::cube batch_push_forward_t4(const arma::cube &t4, Tensor4Type t4type,
                                 const arma::cube &F, bool metric) {
    int N = t4.n_slices;
    bool broadcast = (F.n_slices == 1);
    arma::cube result(6, 6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor4 L(arma::mat::fixed<6,6>(t4.slice(i)), t4type);
        arma::mat::fixed<3,3> Fi(F.slice(broadcast ? 0 : i));
        result.slice(i) = L.push_forward(Fi, metric).mat();
    }
    return result;
}

arma::cube batch_pull_back_t4(const arma::cube &t4, Tensor4Type t4type,
                              const arma::cube &F, bool metric) {
    int N = t4.n_slices;
    bool broadcast = (F.n_slices == 1);
    arma::cube result(6, 6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor4 L(arma::mat::fixed<6,6>(t4.slice(i)), t4type);
        arma::mat::fixed<3,3> Fi(F.slice(broadcast ? 0 : i));
        result.slice(i) = L.pull_back(Fi, metric).mat();
    }
    return result;
}

arma::cube batch_inverse_t4(const arma::cube &t4, Tensor4Type t4type) {
    int N = t4.n_slices;
    arma::cube result(6, 6, N);

    #pragma omp parallel for schedule(static) if(N > 100)
    for (int i = 0; i < N; i++) {
        tensor4 L(arma::mat::fixed<6,6>(t4.slice(i)), t4type);
        result.slice(i) = L.inverse().mat();
    }
    return result;
}

} // namespace simcoon
