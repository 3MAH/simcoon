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

///@file parameter.hpp
///@brief parameters of simcoon
///@version 1.0

#pragma once
#include <armadillo>
#if defined(__cpp_lib_math_constants)
    #include <numbers>
    #define HAS_STD_NUMBERS
#endif

#define UNUSED(x) [&x]{}()

namespace simcoon{

#ifdef HAS_STD_NUMBERS
constexpr double pi = std::numbers::pi;
#else
constexpr double pi = 3.14159265358979323846;
#endif

/// Convert degrees to radians (scalar)
constexpr double deg2rad(double deg) { return deg * pi / 180.0; }
/// Convert radians to degrees (scalar)
constexpr double rad2deg(double rad) { return rad * 180.0 / pi; }

/// Convert degrees to radians (vector)
inline arma::vec deg2rad(const arma::vec& deg) { return deg * pi / 180.0; }
/// Convert radians to degrees (vector)
inline arma::vec rad2deg(const arma::vec& rad) { return rad * 180.0 / pi; }

constexpr int axis_psi = 3;
constexpr int axis_theta = 1;
constexpr int axis_phi = 3;

constexpr double limit = 1.E-9;
constexpr double iota = 1.E-12;
constexpr int miniter_umat = 10;

constexpr int maxiter_umat = 100;
constexpr double precision_umat = 1E-9;

// Tangent-operator mode (state_variables::tangent_mode, solver(), every UMAT).
// Renumbered for 2.0 — pre-2.0 values: 0 = continuum, 1 = algorithmic.
constexpr int tangent_none = 0;          ///< skip tangent assembly, Lt = elastic L (explicit integration)
constexpr int tangent_continuum = 1;     ///< continuum operator (pre-2.0 mode 0)
constexpr int tangent_algorithmic = 2;   ///< Simo-Hughes consistent operator (pre-2.0 mode 1)
constexpr int tangent_closest_point = 3; ///< closest-point-projection integrator with its exact operator (return_mapping.hpp)
constexpr int tangent_default = tangent_algorithmic;

/// Finite-strain formulation carried in the corate_type slot. Rate formulation: corate_type 0..5
/// (corate_kinematics), the law works on \f$ \ln\mathbf{V} \f$ and the Kirchhoff stress
/// \f$ \boldsymbol{\tau} \f$. Total formulation (no rate): the law works on \f$ \ln\mathbf{U} \f$,
/// a state function of \f$ \mathbf{F} \f$, and on its dual stress \f$ \mathbf{T} \f$.
constexpr int corate_type_total = 6;

/// Kinematic form of the finite-strain constitutive update, for the interfaces that name it:
/// @c rate (objective rates, \f$ \ln\mathbf{V} \f$ and \f$ \boldsymbol{\tau} \f$) or @c total
/// (\f$ \ln\mathbf{U} \f$ and its dual stress \f$ \mathbf{T} \f$, no rate). In a corate_type
/// slot the total form is corate_type_total.
enum class KinematicForm { rate, total };

constexpr double div_tnew_dt_umat = 0.2;
constexpr double mul_tnew_dt_umat  = 2.0;
constexpr int maxiter_micro = 100;

constexpr double precision_micro = 1E-6;

} //namespace simcoon
