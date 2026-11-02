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

///@file Tumat_singleT.cpp
///@brief The Abaqus entry software/umat_singleT.cpp as Abaqus calls it: the STATEV layout guard
///       and the failures that must come back as PNEWDT, never as an exception
///@version 1.0

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <armadillo>

#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

using namespace std;
using namespace arma;

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double &sse, double &spd, double &scd, double &rpl, double *ddsddt, double *drplde, double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const double &predef, const double &dpred, char *cmname, const int &ndi, const int &nshr, const int &ntens, const int &nstatev, const double *props, const int &nprops, const double &coords, const double *drot, double &pnewdt, const double &celent, const double *dfgrd0, const double *dfgrd1, const int &noel, const int &npt, const double &layer, const int &kspt, const int &kstep, const int &kinc);

namespace {

const vec props_EPICP = {7.8e-9, 4.6e8, 210000., 0.3, 1.2e-5, 300., 500., 0.4};
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

TEST(Tumat_singleT, nlgeom_rate_form_writes_its_frame)
{
    const int nstatev = 7 + nstatev_EPICP + simcoon::abaqus_kinematic_statev(3);
    const call_result r = call("EPICP_R", nstatev, eye(3, 3), F_stretch);
    EXPECT_EQ(r.pnewdt, 1.);
    EXPECT_TRUE(r.stress.is_finite());
    EXPECT_GT(r.stress(0), 0.);
    // the log_R frame closes STATEV: a rotation, no longer the zero entries of the start
    const mat R(r.statev.data() + nstatev - 9, 3, 3);
    EXPECT_LT(norm(R.t()*R - eye(3, 3), "fro"), 1.e-12);
}

TEST(Tumat_singleT, numerical_failure_is_a_step_cut)
{
    // a singular DFGRD0: ln U_0 of the total form cannot be formed
    const mat F0_singular = diagmat(vec({1., 1., 0.}));
    const call_result r = call("EPICP", 7 + nstatev_EPICP, F0_singular, F_stretch);
    EXPECT_LT(r.pnewdt, 1.);
}

TEST(Tumat_singleTDeathTest, short_depvar_stops_the_analysis)
{
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    // DEPVAR sized for EPICP without the 15 entries _R adds
    EXPECT_DEATH(call("EPICP_R", 7 + nstatev_EPICP, eye(3, 3), F_stretch), "DEPVAR");
    EXPECT_DEATH(call("XXXXX", 7 + nstatev_EPICP, eye(3, 3), F_stretch), "the analysis stops");
}
