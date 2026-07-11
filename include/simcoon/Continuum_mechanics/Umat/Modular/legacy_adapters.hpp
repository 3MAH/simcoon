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

/**
 * @file legacy_adapters.hpp
 * @brief Legacy UMAT names served by the modular engine.
 *
 * The legacy UMATs whose physics ModularUMAT replicates exactly (proven by
 * the equivalence tests in test_modular.py and gated by
 * bench/bench_legacy_vs_modular.py) keep their NAME and PROPS ABI: a
 * per-name translator turns the legacy props vector into the MODUL props
 * stream and the call is forwarded to umat_modular. statev: the modular
 * engine claims the FIRST required slots of the caller's array (legacy
 * allocations are always large enough); trailing legacy slots are left
 * untouched — column MEANING may differ from the deleted implementation
 * (see the migration table in the docs).
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <armadillo>

namespace simcoon {

/// Translate a legacy props vector into the MODUL props stream.
using LegacyPropsTranslator = arma::vec (*)(const arma::vec& legacy_props);

/// True if @p umat_name is served by the modular engine via an adapter.
bool has_legacy_adapter(const std::string& umat_name);

/**
 * @brief UMAT-shaped entry for adapter-served legacy names: translate the
 * props for @p umat_name, then run umat_modular. All other arguments are
 * forwarded untouched (signature identical to umat_modular).
 *
 * @throws std::runtime_error (with the legacy name prefixed) if the name has
 * no registered translator or if the caller's nstatev is smaller than the
 * modular configuration requires.
 */
void umat_legacy_modular(
    const std::string& umat_name,
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    arma::mat& Lt,
    arma::mat& L,
    const arma::mat& DR,
    const int& nprops,
    const arma::vec& props,
    const int& nstatev,
    arma::vec& statev,
    const double& T,
    const double& DT,
    const double& Time,
    const double& DTime,
    double& Wm,
    double& Wm_r,
    double& Wm_ir,
    double& Wm_d,
    const int& ndi,
    const int& nshr,
    const bool& start,
    double& tnew_dt,
    const int& tangent_mode = 0
);

} // namespace simcoon
