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

///@file Tumat_single.hpp
///@brief The Abaqus entries software/umat_singleM.cpp and umat_singleT.cpp as Abaqus calls them:
///       the STATEV layout guard and the failures that must come back as PNEWDT, never as an
///       exception. Included by Tumat_singleM.cpp and Tumat_singleT.cpp (TUMAT_SINGLE_T), one
///       executable per entry since each defines umat_
///@version 1.0

#pragma once

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>
#include "fea_entries.hpp"

using namespace std;
using namespace arma;

namespace {

// EPICP; the thermomechanical entry takes rho and c_p first and keeps Wt after Wm
#if TUMAT_SINGLE_T
const vec props_EPICP = {7.8e-9, 4.6e8, 210000., 0.3, 1.2e-5, 300., 500., 0.4};
const int nheader = 7;
#else
const vec props_EPICP = {210000., 0.3, 1.2e-5, 300., 500., 0.4};
const int nheader = 4;
#endif
const int nstatev_EPICP = 8;

// One 3D call on (DFGRD0, DFGRD1), CMNAME space padded to 80 characters
struct call_result {
    vec stress;
    vector<double> statev;
    double pnewdt;
};

call_result call(const string &name, const int &nstatev, const mat &F0, const mat &F1) {
    call_result r{zeros(6), vector<double>(nstatev, 0.), 1.};
    string cmname = name;
    cmname.resize(80, ' ');
    double ddsdde[36] = {0.}, ddsddt[6] = {0.}, drplde[6] = {0.};
    double sse = 0., spd = 0., scd = 0., rpl = 0., drpldt = 0.;
    const double stran[6] = {0.}, dstran[6] = {0.}, time[2] = {0., 0.};
    const double drot[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    const double predef = 0., dpred = 0., coords = 0., celent = 1., layer = 0.;
    umat_(r.stress.memptr(), r.statev.data(), ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran,
          time, 1., 293.15, 0., predef, dpred, &cmname[0], 3, 3, 6, nstatev, props_EPICP.memptr(),
          static_cast<int>(props_EPICP.n_elem), coords, drot, r.pnewdt, celent, F0.memptr(), F1.memptr(),
          1, 1, layer, 1, 1, 1);
    return r;
}

const mat F_stretch = {{1.002, 0.001, 0.}, {0.001, 0.999, 0.}, {0., 0., 1.}};

}  // namespace

TEST(Tumat_single, nlgeom_rate_form_writes_its_frame)
{
    const int nstatev = nheader + nstatev_EPICP + simcoon::abaqus_kinematic_statev(3);
    const call_result r = call("EPICP_R", nstatev, eye(3, 3), F_stretch);
    EXPECT_EQ(r.pnewdt, 1.);
    EXPECT_TRUE(r.stress.is_finite());
    EXPECT_GT(r.stress(0), 0.);
    // the log_R frame closes STATEV: a rotation, no longer the zero entries of the start
    const mat R(r.statev.data() + nstatev - 9, 3, 3);
    EXPECT_LT(norm(R.t()*R - eye(3, 3), "fro"), 1.e-12);
}

TEST(Tumat_single, numerical_failure_is_a_step_cut)
{
    // a singular trial DFGRD1: its ln U cannot be formed, a smaller increment may cure it
    const mat F1_singular = diagmat(vec({1., 1., 0.}));
    const call_result r = call("EPICP", nheader + nstatev_EPICP, eye(3, 3), F1_singular);
    EXPECT_LT(r.pnewdt, 1.);
}

TEST(Tumat_single, cmname_is_read_whole)
{
    // the dedicated kernel under its full name, rate form: CMNAME is the whole token before the blanks
    const int nstatev = nheader + nstatev_EPICP + simcoon::abaqus_kinematic_statev(3);
    const call_result r = call("EPICP_closed_R", nstatev, eye(3, 3), F_stretch);
    EXPECT_EQ(r.pnewdt, 1.);
    EXPECT_TRUE(r.stress.is_finite());
    EXPECT_GT(r.stress(0), 0.);
    const call_result m = call("EPICP_R", nstatev, eye(3, 3), F_stretch);   // the modular engine, same law
    EXPECT_LT(norm(r.stress - m.stress, 2), 1.e-8 * norm(m.stress, 2));
    const call_result u = call("EPICP_CLOSED_R", nstatev, eye(3, 3), F_stretch);   // as Abaqus upper-cases it
    EXPECT_EQ(norm(u.stress - r.stress, 2), 0.);
}

TEST(Tumat_single, descriptive_tail_falls_back_to_five_letters)
{
    // a deck naming its material EPICP_STEEL runs EPICP (a warning once); the _closed suffix is kept whole
    const call_result r = call("EPICP", nheader + nstatev_EPICP, eye(3, 3), F_stretch);
    const call_result s = call("EPICP_STEEL", nheader + nstatev_EPICP, eye(3, 3), F_stretch);
    EXPECT_EQ(s.pnewdt, 1.);
    EXPECT_EQ(norm(s.stress - r.stress, 2), 0.);
    const int nstatev_R = nheader + nstatev_EPICP + simcoon::abaqus_kinematic_statev(3);
    const call_result rr = call("EPICP_R", nstatev_R, eye(3, 3), F_stretch);
    const call_result sr = call("EPICP_STEEL_R", nstatev_R, eye(3, 3), F_stretch);
    EXPECT_EQ(norm(sr.stress - rr.stress, 2), 0.);
}

TEST(Tumat_singleDeathTest, input_errors_stop_the_analysis)
{
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    // DEPVAR sized for EPICP without the 15 entries _R adds
    EXPECT_DEATH(call("EPICP_R", nheader + nstatev_EPICP, eye(3, 3), F_stretch), "DEPVAR");
    EXPECT_DEATH(call("XXXXX", nheader + nstatev_EPICP, eye(3, 3), F_stretch), "the analysis stops");
    // a singular DFGRD0 is the converged state of the last increment: no step cut cures it
    EXPECT_DEATH(call("EPICP", nheader + nstatev_EPICP, diagmat(vec({1., 1., 0.})), F_stretch), "DFGRD0 is singular");
    EXPECT_DEATH(call("EPICP_R", nheader + nstatev_EPICP + 15, diagmat(vec({1., 1., 0.})), F_stretch), "DFGRD0 is singular");
}
