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
///@file unified_TR.cpp
///@brief Thermomechanical face of the unified SMA model with martensite reorientation

#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/sma_flow.hpp>
#include <simcoon/Continuum_mechanics/Umat/Mechanical/SMA/unified_TR.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_T.hpp>
#include <simcoon/Continuum_mechanics/Umat/Thermomechanical/SMA/unified_TR.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

void umat_sma_unified_TR_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &r, mat &dSdE, vec &dSdT, vec &drdE, double &drdT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, double &Wt, double &Wt_r, double &Wt_ir, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &tangent_mode, const int &ncp) {
    const NCPMethod ncp_m = sma_ncp_method(ncp, tangent_mode, "unified_TR_T");

    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(ndi);
    UNUSED(nshr);

    const double rho = props(0);
    const double c_pA = props(1);
    const double c_pM = props(2);

    if (start) {
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    mat L;
    const SMAIncrement inc = sma_unified_TR_increment(umat_name, Etot, DEtot, sigma, L, DR, props.subvec(3, props.n_elem - 1), statev, T, DT, start, tnew_dt, tangent_mode, ncp_m);
    sma_thermal_face(inc, sigma, DEtot, rho, c_pA, c_pM, T, DT, tangent_mode, r, dSdE, dSdT, drdE, drdT, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir);
}

} //namespace simcoon
