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

#if __has_include(<numbers>)
    #include <numbers>
    #define HAS_STD_NUMBERS
#endif

#define UNUSED(x) [&x]{}()

#ifndef version_full
#define version_full 1
#endif

namespace simcoon{

#ifdef HAS_STD_NUMBERS
    constexpr double pi = std::numbers::pi;
#else
    constexpr double pi = 3.14159265358979323846;
#endif
    
constexpr axis_psi 3
constexpr axis_theta 1
constexpr axis_phi 3

constexpr limit 1.E-9
constexpr iota 1.E-12
constexpr miniter_umat 10

constexpr maxiter_umat 100
constexpr precision_umat 1E-9

constexpr div_tnew_dt_umat 0.2
constexpr mul_tnew_dt_umat 2
constexpr maxiter_micro 100

constexpr precision_micro 1E-6

} //namespace simcoon
