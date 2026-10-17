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

///@file condensed_state.hpp
///@brief Condensed stress states (ndi, nshr): truss, beam fibre, plane stress.
///@version 1.0

#pragma once

#include <armadillo>
#include <functional>

namespace simcoon {

/**
 * @brief Components of a stress state in the Abaqus (ndi, nshr) convention.
 *
 * The kept set \f$ K \f$ is the first ndi direct components followed by the first nshr shears
 * (12, 13, 23), in the order an element passes them (Abaqus STRESS / STRAN / DDSDDE). With
 * ndi < 3 every other component \f$ C \f$ is stress-free and its strain is an unknown of the
 * material point. With ndi = 3 nothing is condensed: the components outside \f$ K \f$ (e.g. the
 * transverse shears of plane strain, ndi = 3, nshr = 1) have zero strain, and \f$ C \f$ is empty.
 */
struct condensed_state {
    arma::uvec K;   ///< kept components, Abaqus order (the element's ntens components)
    arma::uvec C;   ///< stress-free components (empty for ndi = 3)
    bool active() const { return C.n_elem > 0; }
};

/**
 * @brief Abaqus' number of shears for ndi: 0 (truss), 1 (plane stress), 3 (3D).
 * @param ndi Number of direct stress components (1, 2 or 3)
 * @return The default nshr
 */
int default_nshr(const int &ndi);

/**
 * @brief Builds the kept and stress-free sets of the state (ndi, nshr).
 * @param ndi Number of direct stress components (1, 2 or 3)
 * @param nshr Number of shear components (0 to 3); negative picks default_nshr(ndi)
 * @return The condensed state
 * @details Supported condensed states: (1, 0) truss, (1, 2) beam fibre (11, 12, 13),
 * (2, 1) plane stress (11, 22, 12), (2, 3) plane stress with transverse shears. Throws
 * std::invalid_argument on any other ndi or nshr.
 */
condensed_state make_condensed_state(const int &ndi, const int &nshr = -1);

/**
 * @brief Solves the stress-free condition of a condensed state around a 3D law.
 *
 * The material point is the 3D law with the end-of-increment strain on \f$ C \f$ as unknown:
 * \f[
 *    \boldsymbol{\sigma}_C \left( \Delta\boldsymbol{\varepsilon}_K, \Delta\boldsymbol{\varepsilon}_C \right) = \mathbf{0},
 *    \qquad \delta \Delta\boldsymbol{\varepsilon}_C = - \mathbf{L}_{t,CC}^{-1} \, \boldsymbol{\sigma}_C ,
 * \f]
 * a Newton iteration on the law's own tangent. Every evaluation integrates the whole increment
 * from the start state, so the converged point is the 3D integrator itself (its backward-Euler
 * map, its tangent) restricted to the condensed state, whatever the law: the inelastic and
 * viscous strains on \f$ C \f$ are free, and the law never sees an incompatible strain. The
 * Newton is quadratic when \f$ \mathbf{L}_t \f$ is the law's consistent tangent.
 *
 * @param cs The condensed state (must be active)
 * @param DEtot Strain increment (6). On entry its \f$ C \f$ components are the starting guess
 *        (the caller's, zero when it has none); on exit they are the converged ones
 * @param sigma Stress written by evaluate (6), read for the residual
 * @param Lt Tangent written by evaluate (6x6), read for the Newton step
 * @param evaluate Integrates the increment with the current DEtot from the start state (the
 *        callee restores the start state); returns false when the law requests a step cut
 * @param maxiter Maximal number of evaluations
 * @return true on convergence; false on a step-cut request, a singular \f$ \mathbf{L}_{t,CC} \f$,
 *         a non-finite residual or maxiter evaluations (the caller requests a step cut)
 */
bool solve_condensed(const condensed_state &cs, arma::vec &DEtot, const arma::vec &sigma, const arma::mat &Lt,
                     const std::function<bool()> &evaluate, const int &maxiter = 40);

/**
 * @brief Static condensation of a tangent on the kept set.
 * @param cs The condensed state (no-op when inactive)
 * @param Lt Tangent (6x6), replaced by \f$ \mathbf{L}_{KK} - \mathbf{L}_{KC} \mathbf{L}_{CC}^{-1} \mathbf{L}_{CK} \f$
 *        on \f$ K \times K \f$ and zero on every row and column of \f$ C \f$
 * @details The derivative of the kept stress with respect to the kept strain at a stress-free
 * \f$ C \f$: exact when Lt is the law's consistent tangent at the converged point.
 */
void condense_tangent(const condensed_state &cs, arma::mat &Lt);

/**
 * @brief Static condensation of the thermomechanical derivatives on the kept set.
 * @param cs The condensed state (no-op when inactive)
 * @param dSdE Stress-strain tangent (6x6), condensed as by condense_tangent
 * @param dSdT Stress-temperature derivative (6 entries, 6x1 or 1x6)
 * @param drdE Heat-strain derivative (6 entries, 6x1 or 1x6)
 * @param drdT Heat-temperature derivative (1x1)
 * @details With the stress-free condition \f$ \mathrm{d}\boldsymbol{\sigma}_C = \mathbf{0} \f$,
 * \f$ \mathrm{d}\boldsymbol{\varepsilon}_C = - \mathbf{L}_{CC}^{-1} \left( \mathbf{L}_{CK} \mathrm{d}\boldsymbol{\varepsilon}_K + \partial_T \boldsymbol{\sigma}_C \, \mathrm{d}T \right) \f$, so
 \f[
    \partial_T \boldsymbol{\sigma}_K \leftarrow \partial_T \boldsymbol{\sigma}_K - \mathbf{L}_{KC} \mathbf{L}_{CC}^{-1} \partial_T \boldsymbol{\sigma}_C, \quad
    \partial_{\varepsilon_K} r \leftarrow \partial_{\varepsilon_K} r - \partial_{\varepsilon_C} r \cdot \mathbf{L}_{CC}^{-1} \mathbf{L}_{CK}, \quad
    \partial_T r \leftarrow \partial_T r - \partial_{\varepsilon_C} r \cdot \mathbf{L}_{CC}^{-1} \partial_T \boldsymbol{\sigma}_C
 \f]
 * The \f$ C \f$ entries are set to zero.
 */
void condense_tangent_T(const condensed_state &cs, arma::mat &dSdE, arma::mat &dSdT, arma::mat &drdE, arma::mat &drdT);

} //namespace simcoon
