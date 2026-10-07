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
 * @file umat_tutorial_J2.hpp
 * @brief TUTORIAL UMAT — J2 plasticity with linear isotropic hardening,
 * written with the typed Tensor2/Tensor4 API.
 *
 * This is the companion code of the "writing a UMAT" tutorial
 * (docs/simulation/umat_tutorial.rst). The model is chosen so that EVERY
 * quantity has a closed form:
 *
 * - Yield: \f$ f = \sigma_{eq} - (\sigma_Y + H\,p) \f$ with the von Mises
 *   equivalent stress \f$ \sigma_{eq} = \sqrt{\tfrac{3}{2}\,\mathbf{s}:\mathbf{s}} \f$.
 * - Radial return (associated flow, linear hardening):
 *   \f$ \Delta p = f^{trial} / (3\mu + H) \f$ — no local iteration needed.
 * - Continuum tangent, derived analytically:
 *   \f[ \mathbf{L}_t = \mathbf{L}
 *       - \frac{(\mathbf{L}:\boldsymbol{\Lambda}) \otimes (\mathbf{L}:\boldsymbol{\Lambda})}
 *              {\boldsymbol{\Lambda}:\mathbf{L}:\boldsymbol{\Lambda} + H}
 *       \;=\; \mathbf{L} - \frac{4\mu^2}{3\mu + H}\,
 *             \boldsymbol{\Lambda}\otimes\boldsymbol{\Lambda} \f]
 *   since \f$ \mathbf{L}:\boldsymbol{\Lambda} = 2\mu\boldsymbol{\Lambda} \f$
 *   (deviatoric direction) and
 *   \f$ \boldsymbol{\Lambda}:\mathbf{L}:\boldsymbol{\Lambda} = 3\mu \f$.
 *
 * Props (5): E, nu, alpha, sigma_Y, H.
 * statev (8): T_init, p, EP (6, strain Voigt) — identical layout to EPICP, so
 * the two are directly comparable (EPICP with m = 1, k = H is the same model).
 *
 * The tutorial gtest (Ttutorial_umat) validates this implementation three
 * ways: against the EPICP reference, against the shared
 * assemble_continuum_tangent operator, and against a finite-difference
 * Jacobian of the discrete map.
 */

#pragma once

#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>

namespace simcoon {

void umat_tutorial_J2(const std::string &umat_name, const arma::vec &Etot,
                      const arma::vec &DEtot, arma::vec &sigma, arma::mat &Lt,
                      arma::mat &L, const arma::mat &DR, const int &nprops,
                      const arma::vec &props, const int &nstatev,
                      arma::vec &statev, const double &T, const double &DT,
                      const double &Time, const double &DTime, double &Wm,
                      double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi,
                      const int &nshr, const bool &start, double &tnew_dt,
                      const int &tangent_mode = tangent_default);

} // namespace simcoon
