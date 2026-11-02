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

///@file Tusermat_singleM.cpp
///@brief The Ansys entry software/usermat_singleM.cpp as Ansys calls it, against the Abaqus entry
///       software/umat_singleM.cpp on the same path: same stress, tangent and state variables in
///       3D, plane stress and 1D; the input errors stop the analysis
///@version 1.0

#include <gtest/gtest.h>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

#include "fea_entries.hpp"

namespace {

// The simcoon component of each array entry: Abaqus (11, 22, 33, 12, 13, 23), Ansys (11, 22, 33,
// 12, 23, 13) in 3D; both (11, 22, 12) in plane stress and (11) in 1D
struct layout {
    int ndi, nshr, nC;
    vector<int> abaqus, ansys;
};
const layout L3D = {3, 3, 0, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 5, 4}};
const layout LPS = {2, 1, 3, {0, 1, 3}, {0, 1, 3}};
const layout L1D = {1, 0, 5, {0}, {0}};

struct law {
    string name;
    int code;
    vec props;
    int nstatev;
};
const law EPICP = {"EPICP", 5, {210000., 0.3, 1.2e-5, 300., 500., 0.4}, 8};
const law ZENER = {"ZENER", 11, {3000., 0.35, 0., 1500., 0.35, 3000., 1200.}, 8};
const law ELISO = {"ELISO", 2, {210000., 0.3, 0.}, 1};

// One material point, carried over the increments by each entry as the FE code carries it
struct point {
    layout lay;
    law mat;
    int ntens;
    vector<double> stress_aba, statev_aba, stran_aba, stress_ans, statev_ans, strain_ans;
    vector<double> ddsdde, dsdePl;
    double time = 0.;
    arma::mat DR = eye(3, 3);   // rotation increment: DROT (column-major) for Abaqus, rotateM (row-major) for Ansys
    int keycut = 0;
    double pnewdt = 1.;

    point(const layout &l, const law &m)
        : lay(l), mat(m), ntens(static_cast<int>(l.abaqus.size())),
          stress_aba(ntens, 0.), statev_aba(4 + l.nC + m.nstatev, 0.), stran_aba(ntens, 0.),
          stress_ans(ntens, 0.), statev_ans(4 + l.nC + m.nstatev, 0.), strain_ans(ntens, 0.),
          ddsdde(ntens*ntens, 0.), dsdePl(ntens*ntens, 0.) {}

    // the increment dE (simcoon Voigt, engineering shears) over dtime, on both entries (Ansys
    // alone with abaqus = false)
    void step(const vec &dE, const double &dtime, const bool &abaqus = true) {
        vector<double> dstran(ntens), dstrain(ntens);
        for (int i = 0; i < ntens; i++) {
            dstran[i] = dE(lay.abaqus[i]);
            dstrain[i] = dE(lay.ansys[i]);
        }
        // Abaqus, NLGEOM = NO
        const int nstatev = static_cast<int>(statev_ans.size()), nprops = static_cast<int>(mat.props.n_elem);
        if (abaqus) {
            string cmname = mat.name;
            cmname.resize(80, ' ');
            double sse = 0., spd = 0., scd = 0., rpl = 0., drpldt = 0., ddsddt[6] = {0.}, drplde[6] = {0.};
            const double time_a[2] = {time, time}, I[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
            umat_(stress_aba.data(), statev_aba.data(), ddsdde.data(), sse, spd, scd, rpl, ddsddt, drplde, drpldt,
                  stran_aba.data(), dstran.data(), time_a, dtime, 293.15, 0., 0., 0., &cmname[0], lay.ndi, lay.nshr,
                  ntens, nstatev, mat.props.memptr(), nprops, 0., DR.memptr(), pnewdt, 1., I, I, 1, 1, 0., 1, 1, 1);
        }
        // Ansys: the model code first
        vector<double> prop(1 + nprops);
        prop[0] = mat.code;
        std::copy(mat.props.begin(), mat.props.end(), prop.begin() + 1);
        int one = 1, nDirect = lay.ndi, nShear = lay.nshr, ncomp = ntens, nStatev = nstatev, nProp = 1 + nprops;
        double Time = time, dTime = dtime, Temp = 293.15, dTemp = 0., sedEl = 0., sedPl = 0., epseq = 0.;
        vector<double> epsPl(ntens, 0.);
        const arma::mat rot_rows = DR.t();
        double coords[3] = {0.}, rotateM[9], F[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.}, tsstif[2] = {0.};
        std::copy(rot_rows.begin(), rot_rows.end(), rotateM);
        usermat_(&one, &one, &one, &one, &one, &one, &one, &keycut, &nDirect, &nShear, &ncomp, &nStatev, &nProp,
                 &one, &Time, &dTime, &Temp, &dTemp, stress_ans.data(), statev_ans.data(), dsdePl.data(), &sedEl,
                 &sedPl, &epseq, strain_ans.data(), dstrain.data(), epsPl.data(), prop.data(), coords, rotateM,
                 F, F, tsstif);
        for (int i = 0; i < ntens; i++) {
            stran_aba[i] += dstran[i];
            strain_ans[i] += dstrain[i];
        }
        time += dtime;
    }

    // Ansys against Abaqus, entry by simcoon component
    void expect_same(const double &tol) const {
        ASSERT_EQ(pnewdt, 1.);
        ASSERT_EQ(keycut, 0);
        auto ans_of = [&](const int &k) {
            for (int i = 0; i < ntens; i++) if (lay.ansys[i] == k) return i;
            return -1;
        };
        double scale = 0.;
        for (double v : ddsdde) scale = std::max(scale, fabs(v));
        for (int i = 0; i < ntens; i++) {
            const int ia = ans_of(lay.abaqus[i]);
            EXPECT_NEAR(stress_ans[ia], stress_aba[i], tol*(1. + fabs(stress_aba[i])));
            for (int j = 0; j < ntens; j++) {
                const int ja = ans_of(lay.abaqus[j]);
                EXPECT_NEAR(dsdePl[ja*ntens + ia], ddsdde[j*ntens + i], tol*scale);
            }
        }
        for (size_t i = 0; i < statev_aba.size(); i++)
            EXPECT_NEAR(statev_ans[i], statev_aba[i], tol*(1. + fabs(statev_aba[i]))) << "statev " << i;
    }
};

// A loading and unloading path with shears, past the yield stress of EPICP
vec dE_path(const int &inc) {
    const vec dE = {1.5e-3, -0.4e-3, 0.3e-3, 0.8e-3, -0.5e-3, 0.6e-3};
    return (inc < 4) ? vec(dE) : vec(-0.5*dE);
}

}  // namespace

TEST(Tusermat_singleM, matches_the_abaqus_entry)
{
    for (const layout &lay : {L3D, LPS, L1D}) {
        for (const law &mat : {EPICP, ZENER}) {
            point p(lay, mat);
            for (int inc = 0; inc < 6; inc++) {
                p.step(dE_path(inc), 0.5);
                SCOPED_TRACE(mat.name + " ndi " + to_string(lay.ndi) + " inc " + to_string(inc));
                p.expect_same(1.e-12);
            }
        }
    }
}

TEST(Tusermat_singleM, rotation_increment_and_unfilled_rotateM)
{
    // a rotation increment turns the state variables on both entries alike
    point p(L3D, EPICP);
    for (int inc = 0; inc < 6; inc++) {
        p.DR = {{cos(0.05), -sin(0.05), 0.}, {sin(0.05), cos(0.05), 0.}, {0., 0., 1.}};
        p.step(dE_path(inc), 0.5);
        SCOPED_TRACE("inc " + to_string(inc));
        p.expect_same(1.e-12);
    }
    // a zero rotateM (an Ansys release that leaves it unfilled) is no rotation
    point q(L3D, EPICP), r(L3D, EPICP);
    q.DR.zeros();
    for (int inc = 0; inc < 6; inc++) {
        q.step(dE_path(inc), 0.5, false);
        r.step(dE_path(inc), 0.5, false);
    }
    double moved = 0.;
    for (size_t i = 0; i < q.statev_ans.size(); i++) {
        EXPECT_NEAR(q.statev_ans[i], r.statev_ans[i], 1.e-12*(1. + fabs(r.statev_ans[i])));
        moved = std::max(moved, fabs(p.statev_ans[i] - r.statev_ans[i]));
    }
    EXPECT_GT(moved, 1.e-6);   // the rotation did turn the plastic strain
}

TEST(Tusermat_singleM, plane_stress_elastic_tangent)
{
    // the tangent goes out once condensed: E/(1 - nu^2) on 11, 11, never 0/0
    point p(LPS, ELISO);
    p.step(vec({1.e-3, 0., 0., 0., 0., 0.}), 1.);
    const double E = 210000., nu = 0.3;
    EXPECT_NEAR(p.dsdePl[0], E/(1. - nu*nu), 1.e-6*E);
    EXPECT_NEAR(p.dsdePl[3*1 + 0], nu*E/(1. - nu*nu), 1.e-6*E);
    EXPECT_NEAR(p.dsdePl[3*2 + 2], E/(2.*(1. + nu)), 1.e-6*E);
    EXPECT_NEAR(p.stress_ans[0], E/(1. - nu*nu)*1.e-3, 1.e-9*E);
    // the solved 33 strain is kept in USTATEV after the works
    EXPECT_NEAR(p.statev_ans[4], -nu/(1. - nu)*1.e-3, 1.e-12);
}

TEST(Tusermat_singleMDeathTest, input_errors_stop_the_analysis)
{
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    // unknown model code
    law unknown = EPICP;
    unknown.code = 999;
    EXPECT_DEATH({ point p(L3D, unknown); p.step(dE_path(0), 1., false); }, "unknown model code");
    // USTATEV without room for the plane-stress free strain: the Wm header only
    EXPECT_DEATH({
        point p(LPS, EPICP);
        p.statev_ans.assign(4, 0.);
        p.step(dE_path(0), 1., false);
    }, "below the");
    // a beam layout (nDirect = 1, nShear = 2)
    const layout beam = {1, 2, 3, {0, 3, 4}, {0, 3, 4}};
    EXPECT_DEATH({ point p(beam, EPICP); p.step(dE_path(0), 1., false); }, "not supported");
}
