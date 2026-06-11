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

///@file fibre_Fp.hpp
///@brief Single-fibre transversely-isotropic plasticity with an explicit
///plastic deformation gradient Fp (multiplicative decomposition F = Fe Fp,
///exponential-map update) and a hyperelastic Hencky stress evaluated on the
///true elastic logarithmic strain ln(Ve), Ve from Fe = F Fp^{-1}.
///Requires the natural basis carrying F (solver corate_type = 5).

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

namespace simcoon {

/**
 * @brief UMAT "EPTRF": fibre plasticity with explicit Fp.
 *
 * props (12): E, nu, alpha, sigmaY, k, m, c1, c2, a0_x, a0_y, a0_z, anchor
 *   - yield: Phi^2 = 3/2 s:s + c1 (a.s.a)^2 + c2 [(s.a).(s.a) - (a.s.a)^2],
 *     isotropic hardening sigmaY + k p^m (same criterion as EPTRI).
 *   - anchor = 0: embedded material line, a = F a0 / |F a0|  (rides Fe Fp).
 *     anchor = 1: lattice direction,      a = Fe a0 / |Fe a0| (rides Fe only;
 *     plastic flow leaves the direction invariant).
 *
 * statev (11): T_init, p, Fp(3x3 row-major).
 *
 * The stress is an exact function of the state (F, Fp): the elastic branch is
 * the Hencky potential on ln(Ve), so closed elastic cycles return identically
 * (no O(Dt) hypoelastic residual). Plastic flow updates Fp by the exponential
 * map in the relaxed configuration (zero plastic spin), preserving det Fp = 1.
 */
void umat_plasticity_fibre_Fp(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR, const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const natural_basis &nb, const int &tangent_mode = 0);

} //namespace simcoon
