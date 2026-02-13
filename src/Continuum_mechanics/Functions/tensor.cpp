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
#include <simcoon/Simulation/Maths/rotation.hpp>

namespace simcoon {

// Helper: compute inv of a 3x3 fixed matrix, returning fixed
static arma::mat::fixed<3,3> inv33(const arma::mat::fixed<3,3> &F) {
    arma::mat tmp = arma::inv(arma::mat(F));
    arma::mat::fixed<3,3> result;
    result = tmp;
    return result;
}

// ============================================================================
// tensor2 implementation
// ============================================================================

tensor2::tensor2() : _mat(arma::fill::zeros), _vtype(VoigtType::stress),
                     _symmetric(true), _symmetry_checked(true) {}

tensor2::tensor2(VoigtType vtype) : _mat(arma::fill::zeros), _vtype(vtype),
                                     _symmetric(true), _symmetry_checked(true) {}

tensor2::tensor2(const arma::mat::fixed<3,3> &m, VoigtType vtype)
    : _mat(m), _vtype(vtype), _symmetric(false), _symmetry_checked(false) {}

tensor2::tensor2(const arma::mat &m, VoigtType vtype)
    : _vtype(vtype), _symmetric(false), _symmetry_checked(false) {
    assert(m.n_rows == 3 && m.n_cols == 3);
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
    result._symmetric = true;
    result._symmetry_checked = true;
    return result;
}

tensor2 tensor2::from_voigt(const arma::vec &v, VoigtType vtype) {
    assert(v.n_elem == 6);
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
    _symmetry_checked = false;
    return _mat;
}

void tensor2::set_mat(const arma::mat::fixed<3,3> &m) {
    _mat = m;
    _symmetry_checked = false;
}

void tensor2::set_mat(const arma::mat &m) {
    assert(m.n_rows == 3 && m.n_cols == 3);
    _mat = m;
    _symmetry_checked = false;
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
    _symmetric = true;
    _symmetry_checked = true;
}

void tensor2::set_voigt(const arma::vec &v) {
    assert(v.n_elem == 6);
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

Fastor::TensorMap<double,3,3> tensor2::fastor() {
    return arma_to_fastor2(_mat);
}

Fastor::TensorMap<const double,3,3> tensor2::fastor() const {
    return arma_to_fastor2(_mat);
}

bool tensor2::is_symmetric(double tol) const {
    if (_symmetry_checked) return _symmetric;

    _symmetric = (std::abs(_mat(0,1) - _mat(1,0)) < tol) &&
                 (std::abs(_mat(0,2) - _mat(2,0)) < tol) &&
                 (std::abs(_mat(1,2) - _mat(2,1)) < tol);
    _symmetry_checked = true;
    return _symmetric;
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

tensor2 tensor2::push_forward(const arma::mat::fixed<3,3> &F) const {
    switch (_vtype) {
        case VoigtType::stress: {
            arma::mat::fixed<3,3> result;
            result = F * _mat * F.t();
            return tensor2(result, VoigtType::stress);
        }
        case VoigtType::strain: {
            arma::mat::fixed<3,3> invF = inv33(F);
            arma::mat::fixed<3,3> result;
            result = invF.t() * _mat * invF;
            return tensor2(result, VoigtType::strain);
        }
        case VoigtType::generic:
        case VoigtType::none:
            throw std::runtime_error("push_forward requires VoigtType::stress or VoigtType::strain");
    }
    return *this;
}

tensor2 tensor2::pull_back(const arma::mat::fixed<3,3> &F) const {
    switch (_vtype) {
        case VoigtType::stress: {
            arma::mat::fixed<3,3> invF = inv33(F);
            arma::mat::fixed<3,3> result;
            result = invF * _mat * invF.t();
            return tensor2(result, VoigtType::stress);
        }
        case VoigtType::strain: {
            arma::mat::fixed<3,3> result;
            result = F.t() * _mat * F;
            return tensor2(result, VoigtType::strain);
        }
        case VoigtType::generic:
        case VoigtType::none:
            throw std::runtime_error("pull_back requires VoigtType::stress or VoigtType::strain");
    }
    return *this;
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

tensor2 tensor2::operator*(double scalar) const {
    arma::mat::fixed<3,3> result;
    result = _mat * scalar;
    return tensor2(result, _vtype);
}

tensor2& tensor2::operator+=(const tensor2 &other) {
    _mat += other._mat;
    _symmetry_checked = false;
    return *this;
}

tensor2& tensor2::operator-=(const tensor2 &other) {
    _mat -= other._mat;
    _symmetry_checked = false;
    return *this;
}

tensor2& tensor2::operator*=(double scalar) {
    _mat *= scalar;
    return *this;
}

tensor2 operator*(double scalar, const tensor2 &t) {
    arma::mat::fixed<3,3> result;
    result = t._mat * scalar;
    return tensor2(result, t._vtype);
}

bool tensor2::operator==(const tensor2 &other) const {
    return _vtype == other._vtype && arma::approx_equal(_mat, other._mat, "absdiff", 1e-14);
}

// Free functions
tensor2 stress(const arma::mat::fixed<3,3> &m) {
    return tensor2(m, VoigtType::stress);
}

tensor2 stress(const arma::vec::fixed<6> &v) {
    return tensor2::from_voigt(v, VoigtType::stress);
}

tensor2 strain(const arma::mat::fixed<3,3> &m) {
    return tensor2(m, VoigtType::strain);
}

tensor2 strain(const arma::vec::fixed<6> &v) {
    return tensor2::from_voigt(v, VoigtType::strain);
}

arma::vec::fixed<6> dev(const tensor2 &t) {
    arma::vec::fixed<6> v = t.voigt();
    double tr_val = v(0) + v(1) + v(2);
    arma::vec::fixed<6> result = v;
    result(0) -= tr_val / 3.0;
    result(1) -= tr_val / 3.0;
    result(2) -= tr_val / 3.0;
    return result;
}

double Mises(const tensor2 &t) {
    arma::vec::fixed<6> d = dev(t);
    if (t.vtype() == VoigtType::stress || t.vtype() == VoigtType::generic) {
        return std::sqrt(1.5 * (d(0)*d(0) + d(1)*d(1) + d(2)*d(2) +
                                2.0*(d(3)*d(3) + d(4)*d(4) + d(5)*d(5))));
    } else if (t.vtype() == VoigtType::strain) {
        return std::sqrt(2.0/3.0 * (d(0)*d(0) + d(1)*d(1) + d(2)*d(2) +
                                    0.5*(d(3)*d(3) + d(4)*d(4) + d(5)*d(5))));
    }
    throw std::runtime_error("Mises not defined for VoigtType::none");
}

double trace(const tensor2 &t) {
    return t.mat()(0,0) + t.mat()(1,1) + t.mat()(2,2);
}

// ============================================================================
// tensor4 implementation
// ============================================================================

tensor4::tensor4() : _voigt(arma::fill::zeros), _type(Tensor4Type::stiffness),
                     _fastor_valid(false) {}

tensor4::tensor4(Tensor4Type type) : _voigt(arma::fill::zeros), _type(type),
                                      _fastor_valid(false) {}

tensor4::tensor4(const arma::mat::fixed<6,6> &m, Tensor4Type type)
    : _voigt(m), _type(type), _fastor_valid(false) {}

tensor4::tensor4(const arma::mat &m, Tensor4Type type)
    : _type(type), _fastor_valid(false) {
    assert(m.n_rows == 6 && m.n_cols == 6);
    _voigt = m;
}

tensor4::tensor4(const tensor4 &other)
    : _voigt(other._voigt), _type(other._type),
      _fastor(other._fastor), _fastor_valid(other._fastor_valid) {}

tensor4::tensor4(tensor4 &&other) noexcept
    : _voigt(std::move(other._voigt)), _type(other._type),
      _fastor(std::move(other._fastor)), _fastor_valid(other._fastor_valid) {}

tensor4& tensor4::operator=(const tensor4 &other) {
    if (this != &other) {
        _voigt = other._voigt;
        _type = other._type;
        _fastor = other._fastor;
        _fastor_valid = other._fastor_valid;
    }
    return *this;
}

tensor4& tensor4::operator=(tensor4 &&other) noexcept {
    if (this != &other) {
        _voigt = std::move(other._voigt);
        _type = other._type;
        _fastor = std::move(other._fastor);
        _fastor_valid = other._fastor_valid;
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
    assert(m.n_rows == 6 && m.n_cols == 6);
    _voigt = m;
    _invalidate_fastor();
}

void tensor4::_ensure_fastor() const {
    if (!_fastor_valid) {
        _fastor = voigt_to_fastor4(_voigt);
        _fastor_valid = true;
    }
}

const Fastor::Tensor<double,3,3,3,3>& tensor4::fastor() const {
    _ensure_fastor();
    return _fastor;
}

tensor2 tensor4::contract(const tensor2 &t) const {
    arma::vec::fixed<6> v_in = t.voigt();
    arma::vec::fixed<6> v_out;

    for (int i = 0; i < 6; ++i) {
        double sum = 0.0;
        for (int j = 0; j < 6; ++j) {
            sum += _voigt(i,j) * v_in(j);
        }
        v_out(i) = sum;
    }

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

tensor4 tensor4::push_forward(const arma::mat::fixed<3,3> &F) const {
    _ensure_fastor();

    Fastor::Tensor<double,3,3> F_fastor;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            F_fastor(i,j) = F(i,j);

    Fastor::Tensor<double,3,3,3,3> result = push_forward_4(_fastor, F_fastor);
    arma::mat::fixed<6,6> voigt_result = fastor4_to_voigt(result);
    return tensor4(voigt_result, _type);
}

tensor4 tensor4::pull_back(const arma::mat::fixed<3,3> &F) const {
    _ensure_fastor();

    arma::mat::fixed<3,3> invF = inv33(F);

    Fastor::Tensor<double,3,3> invF_fastor;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            invF_fastor(i,j) = invF(i,j);

    Fastor::Tensor<double,3,3,3,3> result = pull_back_4(_fastor, invF_fastor);
    arma::mat::fixed<6,6> voigt_result = fastor4_to_voigt(result);
    return tensor4(voigt_result, _type);
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

tensor4 tensor4::operator*(double scalar) const {
    arma::mat::fixed<6,6> result;
    result = _voigt * scalar;
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

tensor4 operator*(double scalar, const tensor4 &t) {
    arma::mat::fixed<6,6> result;
    result = t._voigt * scalar;
    return tensor4(result, t._type);
}

bool tensor4::operator==(const tensor4 &other) const {
    return _type == other._type && arma::approx_equal(_voigt, other._voigt, "absdiff", 1e-14);
}

// Free functions
tensor4 dyadic(const tensor2 &a, const tensor2 &b) {
    arma::vec::fixed<6> a_v, b_v;

    const arma::mat::fixed<3,3> &am = a.mat();
    const arma::mat::fixed<3,3> &bm = b.mat();

    a_v(0) = am(0,0); a_v(1) = am(1,1); a_v(2) = am(2,2);
    a_v(3) = 0.5*(am(0,1) + am(1,0));
    a_v(4) = 0.5*(am(0,2) + am(2,0));
    a_v(5) = 0.5*(am(1,2) + am(2,1));

    b_v(0) = bm(0,0); b_v(1) = bm(1,1); b_v(2) = bm(2,2);
    b_v(3) = 0.5*(bm(0,1) + bm(1,0));
    b_v(4) = 0.5*(bm(0,2) + bm(2,0));
    b_v(5) = 0.5*(bm(1,2) + bm(2,1));

    arma::mat::fixed<6,6> C;
    C = a_v * b_v.t();
    return tensor4(C, Tensor4Type::stiffness);
}

tensor4 auto_dyadic(const tensor2 &a) {
    return dyadic(a, a);
}

} // namespace simcoon
