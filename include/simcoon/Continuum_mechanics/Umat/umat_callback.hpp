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

///@file umat_callback.hpp
///@brief Process-wide callback slot served by the "PYEXT" UMAT name: a small-strain
///       constitutive law implemented outside libsimcoon (typically in Python) and
///       driven by the simcoon solver or by the batch UMAT entry point.
///@version 2.0

#pragma once

#include <functional>
#include <string>
#include <armadillo>

namespace simcoon {

/**
 * @brief Signature of a small-strain mechanical UMAT (the 24-argument simcoon convention).
 *
 * Identical to the built-in kernels (e.g. umat_plasticity_iso_CCP): the callee receives the
 * strain at the beginning of the increment \f$ \boldsymbol{\varepsilon}_n \f$ (@p Etot), the
 * increment \f$ \Delta\boldsymbol{\varepsilon} \f$ (@p DEtot), the stress at the beginning of
 * the increment (@p sigma, in/out), the rotation increment @p DR, the material properties,
 * the internal variables (@p statev, in/out), the temperature and time data, the accumulated
 * energies (in/out) and must return the updated stress, the tangent operator
 * \f$ \mathbf{L}_t = \partial \boldsymbol{\sigma} / \partial \boldsymbol{\varepsilon} \f$
 * (@p Lt), the elastic operator (@p L), the updated internal variables and energies.
 *
 * Voigt order is 11, 22, 33, 12, 13, 23 with engineering shear strains, all quantities
 * expressed in the material (local) frame. Under finite strain the dispatcher feeds the
 * logarithmic strain and expects the Kirchhoff stress on output ("PYEXT" belongs to the
 * kirchhoff_box set of stress_output_is_kirchhoff()), exactly like the built-in small-strain
 * kernels.
 *
 * Contract (see the solver): the callee must be a pure function of its arguments — the solver
 * restores @p statev, @p sigma and the energies to their start-of-increment values before
 * every Newton retrial of the same increment. A recurrent model must therefore keep its
 * hidden state inside @p statev. @p start is true on the first call of a block
 * (initialise @p statev there). Writing @p tnew_dt < 1 requests a smaller increment.
 */
using umat_M_callback = std::function<void(
    const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot,
    arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR,
    const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev,
    const double &T, const double &DT, const double &Time, const double &DTime,
    double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
    const int &ndi, const int &nshr, const bool &start, double &tnew_dt,
    const int &tangent_mode)>;

/**
 * @brief Install the process-wide callback served by the "PYEXT" UMAT name.
 * @param cb Callable following umat_M_callback. An empty function clears the slot.
 *
 * Thread-safe (mutex). The callback itself is invoked outside the lock, serially, by the
 * solver (one material point) or by the serial branch of the batch UMAT entry point. The
 * Python bindings register their bridge here; a C++ program may register a plain lambda.
 */
void set_umat_callback(umat_M_callback cb);

/// @brief Clear the process-wide callback (subsequent "PYEXT" calls throw exception_solver).
void clear_umat_callback();

/// @brief True when a callback is installed.
bool has_umat_callback();

/// @brief Snapshot (copy) of the installed callback; empty when none is installed.
umat_M_callback get_umat_callback();

/**
 * @brief UMAT entry point of the "PYEXT" name: forwards the 24-argument call to the installed
 *        callback.
 * @throws exception_solver when no callback is installed.
 *
 * Has exactly the small-strain UMAT signature so that it can be stored in the same function
 * pointer as the built-in kernels by the dispatchers.
 */
void umat_callback_M(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot,
                     arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR,
                     const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev,
                     const double &T, const double &DT, const double &Time, const double &DTime,
                     double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
                     const int &ndi, const int &nshr, const bool &start, double &tnew_dt,
                     const int &tangent_mode);

} // namespace simcoon
