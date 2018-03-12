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

#include <boost/math/constants/constants.hpp>
#define UNUSED(x) [&x]{}()

#ifndef version_full
#define version_full 1
#endif

namespace simcoon{

#ifndef pi
#define pi boost::math::constants::pi<double>()
#endif
    
#ifndef axis_psi
#define axis_psi 3
#endif

#ifndef axis_theta
#define axis_theta 1
#endif

#ifndef axis_phi
#define axis_phi 3
#endif

#ifndef limit
#define limit 1E-9
#endif

#ifndef iota
#define iota 1E-12
#endif

#ifndef miniter_umat
#define miniter_umat 10
#endif

#ifndef maxiter_umat
#define maxiter_umat 100
#endif

#ifndef precision_umat
#define precision_umat 1E-9
#endif

#ifndef div_tnew_dt_umat
#define div_tnew_dt_umat 0.2
#endif

#ifndef mul_tnew_dt_umat
#define mul_tnew_dt_umat 2
#endif

#ifndef maxiter_micro
#define maxiter_micro 100
#endif

#ifndef precision_micro
#define precision_micro 1E-6
#endif

} //namespace simcoon
