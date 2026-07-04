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

///@file chaboche_cpp.hpp
///@brief Shared closest-point-projection (tangent_mode = 2) local solve for the
///Chaboche family: Voce isotropic hardening + two Armstrong-Frederick backstresses,
///parameterised by the yield criterion / flow-direction pair. One implementation
///serves EPCHA (Mises), EPHAC (Hill), EPANI (Ani) and EPDFA (DFA) — the UMATs differ
///only in the injected `crit`/`flow` callables and in their success/failure blocks.
///@version 1.0

#pragma once

#include <functional>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Umat/return_mapping.hpp>

namespace simcoon {

///@brief CPP local solve for a 1-mechanism Chaboche model (doc §cpp_return_mapping).
///
///Backward-Euler closed forms resolve the state at every iterate:
///  a_i = (a_{n,i} + Dl*flow(xi)) / (1 + Dl*D_i)   (Armstrong-Frederick recall)
///  Hp  = (Hp_n + b*Q*Dl) / (1 + b*Dl)             (Voce isotropic hardening)
///iterated on xi = sigma - X. ALL derivative callbacks handed to the return-mapping
///helper are TOTAL derivatives of the inner-consistent composed map (central FD with
///save/restore of the referenced state) — partial-only inputs would leave an
///O(dX/dsigma) error in the consistent tangent.
///
///@param sigma_tr  trial stress (Voigt 6)
///@param L         elastic stiffness
///@param crit      equivalent-stress criterion evaluated on xi = sigma - X
///@param flow      flow direction (criterion gradient) evaluated on xi = sigma - X
///@param sigmaY,Q,b,C_1,D_1,C_2,D_2  Chaboche material parameters
///@param p_n,a_1n,a_2n,X_1n,X_2n,Hp_n  start-of-increment state (frozen copies)
///@param p,Hp,dHpdp,a_1,a_2,X_1,X_2,X  live state references — on a CONVERGED return
///        they hold the inner-consistent end-of-increment state; on failure they hold
///        whatever the last iterate produced (the caller's failure block restores).
ReturnMappingResult chaboche_cpp_return_mapping(
    const arma::vec &sigma_tr,
    const arma::mat &L,
    const std::function<double(const arma::vec &)> &crit,
    const std::function<arma::vec(const arma::vec &)> &flow,
    const double &sigmaY, const double &Q, const double &b,
    const double &C_1, const double &D_1, const double &C_2, const double &D_2,
    const double &p_n, const arma::vec &a_1n, const arma::vec &a_2n,
    const arma::vec &X_1n, const arma::vec &X_2n, const double &Hp_n,
    double &p, double &Hp, double &dHpdp,
    arma::vec &a_1, arma::vec &a_2, arma::vec &X_1, arma::vec &X_2, arma::vec &X);

} // namespace simcoon
