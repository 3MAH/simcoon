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

///@file umat_callback.cpp
///@brief Process-wide callback slot served by the "PYEXT" UMAT name.

#include <mutex>
#include <utility>

#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_callback.hpp>

namespace simcoon {

namespace {

// Function-local statics: no static-initialization-order issue, and the slot is never an
// exported symbol (CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS only exports the accessor functions).
std::mutex &cb_mutex() {
    static std::mutex m;
    return m;
}

umat_M_callback &cb_slot() {
    static umat_M_callback cb;
    return cb;
}

} // namespace

void set_umat_callback(umat_M_callback cb) {
    std::lock_guard<std::mutex> lock(cb_mutex());
    cb_slot() = std::move(cb);
}

void clear_umat_callback() {
    set_umat_callback(umat_M_callback{});
}

bool has_umat_callback() {
    std::lock_guard<std::mutex> lock(cb_mutex());
    return static_cast<bool>(cb_slot());
}

umat_M_callback get_umat_callback() {
    std::lock_guard<std::mutex> lock(cb_mutex());
    return cb_slot();
}

void umat_callback_M(const std::string &umat_name, const arma::vec &Etot, const arma::vec &DEtot,
                     arma::vec &sigma, arma::mat &Lt, arma::mat &L, const arma::mat &DR,
                     const int &nprops, const arma::vec &props, const int &nstatev, arma::vec &statev,
                     const double &T, const double &DT, const double &Time, const double &DTime,
                     double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d,
                     const int &ndi, const int &nshr, const bool &start, double &tnew_dt,
                     const int &tangent_mode)
{
    // Snapshot under the lock, call outside it (the callback may take a long time).
    umat_M_callback cb = get_umat_callback();
    if (!cb) {
        throw exception_solver(
            "PYEXT: no UMAT callback registered. From Python, pass the law object to "
            "simcoon.solver.solve(blocks, umat, ...) or wrap the call in "
            "`with simcoon.registered(umat): ...`.");
    }
    cb(umat_name, Etot, DEtot, sigma, Lt, L, DR, nprops, props, nstatev, statev,
       T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, ndi, nshr, start, tnew_dt, tangent_mode);
}

} // namespace simcoon
