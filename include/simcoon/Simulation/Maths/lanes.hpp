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

///@file lanes.hpp
///@brief Two-lane double arithmetic: a SIMD vector type where the compiler offers one, a plain pair elsewhere
///@version 1.0

#pragma once

namespace simcoon{

/** @addtogroup maths
 *  @{
 */

/**
 * @brief Two doubles operated on together.
 *
 * On GCC and Clang `lane2` is a 16-byte vector type: `lanes`, `mul` and `mul_add` compile to one
 * SIMD instruction each, and a loop written on `lane2` serves two scalar sums per operation while
 * each lane follows the scalar expressions exactly (same rounding as the scalar code). On other
 * compilers `lane2` is a pair of doubles and the same functions run lane by lane: same results,
 * no speed-up. This is the only compiler-specific definition of the library; the kernels that
 * use it (the Eshelby quadrature) stay free of preprocessor conditions.
 */
#if defined(__GNUC__) || defined(__clang__)
typedef double lane2 __attribute__((vector_size(16)));
/// \f$ (a, b) \f$
inline lane2 lanes(const double &a, const double &b) { lane2 r = {a, b}; return r; }
/// lane by lane \f$ a_k b_k \f$
inline lane2 mul(const lane2 &a, const lane2 &b) { return a*b; }
/// lane by lane \f$ a_k b_k + c_k \f$ (the compiler may contract it into one fma, as for the scalar expression)
inline lane2 mul_add(const lane2 &a, const lane2 &b, const lane2 &c) { return a*b + c; }
#else
struct lane2 {
    double v[2];
    double operator[](int i) const { return v[i]; }
};
inline lane2 lanes(const double &a, const double &b) { lane2 r = {{a, b}}; return r; }
inline lane2 mul(const lane2 &a, const lane2 &b) { lane2 r = {{a.v[0]*b.v[0], a.v[1]*b.v[1]}}; return r; }
inline lane2 mul_add(const lane2 &a, const lane2 &b, const lane2 &c) { lane2 r = {{a.v[0]*b.v[0] + c.v[0], a.v[1]*b.v[1] + c.v[1]}}; return r; }
#endif

/** @}*/

} //namespace simcoon
