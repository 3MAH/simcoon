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
 * The legacy UMATs whose physics ModularUMAT replicates (proven by the
 * equivalence tests in test_modular.py / test_modular_thermo.py and gated by
 * benchmark/bench_legacy_vs_modular.py) keep their NAME and PROPS ABI: a
 * per-name translator turns the legacy props vector into the MODUL props
 * stream and the call is forwarded to umat_modular (mechanical face) or
 * umat_modular_T (thermomechanical face, the legacy props being
 * `[rho, c_p, mechanical legacy props]`). statev: the modular engine claims
 * the FIRST required slots of the caller's array (legacy allocations are
 * always large enough); trailing legacy slots are left untouched — column
 * MEANING may differ from the dedicated implementation (table below). The
 * dedicated kernels that still exist answer the `_closed` names (EPICP_closed,
 * EPCHA_closed, PRONK_closed, ZENER_closed, EPKCP_closed on the thermomechanical face,
 * NEOHC_closed and NEOHI_closed built on F).
 *
 * The hyperelastic names run the HYPER_INVARIANTS block of the modular engine, a potential of
 * the logarithmic strain: the finite dispatch serves them under corate 3 (log_R) and the total
 * formulation only (adapter_is_hyperelastic); NEOHC_closed / NEOHI_closed, built on F, serve
 * every corate and the native material tangent.
 *
 * Name  | legacy props                               | MODUL stream (translate_NAME)                            | modular statev
 * ------|--------------------------------------------|----------------------------------------------------------|---------------
 * EPICP | E nu alpha sigmaY k m                      | isotropic block; plasticity, Mises, power law, no kin    | [T_init, p, EP] (8)
 * EPKCP | E nu alpha sigmaY k m kX                   | ..., Prager C = 1.5 kX                                   | [T_init, p, EP, a] (14)
 * EPCHA | E nu alpha sigmaY Q b C1 D1 C2 D2          | ..., Voce, Chaboche (C1, D1), (C2, D2)                   | [T_init, p, EP, a_1, a_2] (20)
 * PRONK | E0 nu0 alpha N (E_i nu_i etaB_i etaS_i)xN  | isotropic block; one viscoelastic mechanism, N branches  | [T_init, (v_i, EV_i) x N] (1 + 7N)
 * ZENER | E0 nu0 alpha E1 nu1 etaB1 etaS1            | isotropic block; one Maxwell branch, the Kelvin unit mapped per isotropic mode X in {K, mu}: X_M = X0^2/(X0 + X1), eta_M = eta_K X0^2/(X0 + X1)^2 | [T_init, v, EV] (8), EV the Maxwell branch strain (not the Kelvin viscous strain of ZENER_closed)
 * NEOHC, MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA | the potential's own parameters [+ volumetric selector] (HyperPotential) | HYPER_INVARIANTS block, the same HyperPotential, alpha = 0 | [T_init] (1)
 * NEOHI | E nu alpha                                 | HYPER_INVARIANTS block, NEOHC potential [mu = E/2(1+nu), kappa = E/3(1-2nu)], quadratic U(J), alpha | [T_init] (1)
 *
 * The other families (ELISO/ELIST/ELORT, EPHIL/EPTRI/EPHAC/EPANI/EPDFA/EPCHG/EPHIN) are listed in
 * docs/simulation/umat_catalog.rst; a_i are back-strains (engineering shear), EV_i branch strains.
 *
 * @version 1.0
 */

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon {

/// Translate a legacy props vector into the MODUL props stream.
using LegacyPropsTranslator = arma::vec (*)(const arma::vec& legacy_props);

/// True if @p umat_name is served by the modular engine via an adapter.
bool has_legacy_adapter(const std::string& umat_name);

/// True if the adapter of @p umat_name is a hyperelastic potential of the logarithmic strain
/// (NEOHC, NEOHI, MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA): the finite dispatch serves it under
/// corate 3 (log_R) and the total formulation only, as MODUL.
bool adapter_is_hyperelastic(const std::string& umat_name);

/**
 * @brief UMAT-shaped entry for adapter-served legacy names: translate the
 * props for @p umat_name, then run umat_modular. All other arguments are
 * forwarded untouched (umat_modular's signature). @p ncp, the NCPMethod code of the run, is the
 * complementarity function of the cutting plane (EPKCP); @p kinematic_form is the form of the
 * finite dispatch (the hyperelastic block returns the stress dual to H under the total form).
 *
 * @throws std::invalid_argument (with the legacy name prefixed) if the name has
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
    const int& tangent_mode = tangent_default,
    const int& ncp = 0
,
    const KinematicForm& kinematic_form = KinematicForm::rate);

/**
 * @brief The thermomechanical counterpart of umat_legacy_modular: @p props is
 * `[rho, c_p, legacy mechanical props]`, translated to `[rho, c_p, MODUL stream]` and
 * forwarded to umat_modular_T (its signature, @p ncp included). Serves EPICP, EPKCP and
 * PRONK. Same exceptions as the mechanical face.
 */
void umat_legacy_modular_T(
    const std::string& umat_name,
    const arma::vec& Etot,
    const arma::vec& DEtot,
    arma::vec& sigma,
    double& r,
    arma::mat& dSdE,
    arma::vec& dSdT,
    arma::vec& drdE,
    double& drdT,
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
    double& Wt,
    double& Wt_r,
    double& Wt_ir,
    const int& ndi,
    const int& nshr,
    const bool& start,
    double& tnew_dt,
    const int& tangent_mode,
    const KinematicForm& kinematic_form,
    const int& ncp);

} // namespace simcoon
