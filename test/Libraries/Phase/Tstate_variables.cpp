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

///@file Tstate_variables.cpp
///@brief Test the copies of the state variables of a phase
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

namespace {

//Distinct, non-default values in every field of the common part, statev included
template <class SV>
void fill_common(SV &sv)
{
    sv.tangent_mode = 3;
    sv.ncp = 1;
    sv.Etot = {1., 2., 3., 4., 5., 6.};
    sv.DEtot = sv.Etot + 0.1;
    sv.etot = sv.Etot + 0.2;
    sv.Detot = sv.Etot + 0.3;
    sv.sigma = {10., 20., 30., 40., 50., 60.};
    sv.sigma_start = sv.sigma + 1.;
    sv.PKII = sv.sigma + 2.;
    sv.PKII_start = sv.sigma + 3.;
    sv.tau = sv.sigma + 4.;
    sv.tau_start = sv.sigma + 5.;
    sv.F0 = {{1.1, 0.1, 0.}, {0., 1., 0.2}, {0., 0., 0.9}};
    sv.F1 = sv.F0 + 0.01;
    sv.U0 = {{1.2, 0., 0.}, {0., 1.1, 0.}, {0., 0., 0.8}};
    sv.U1 = sv.U0 + 0.02;
    sv.R = {{0., -1., 0.}, {1., 0., 0.}, {0., 0., 1.}};
    sv.DR = sv.R.t();
    sv.T = 300.;
    sv.DT = 5.;
    sv.resize(3, true, 1.5);
    sv.statev_start = {7., 8., 9.};
}

//full = copy constructor or operator= (tangent_mode and statev included); otherwise copy_fields
template <class SV>
void expect_common(const SV &a, const SV &b, bool full)
{
    if (full) {
        EXPECT_EQ(a.tangent_mode, b.tangent_mode);
        EXPECT_EQ(a.ncp, b.ncp);
        EXPECT_EQ(a.nstatev, b.nstatev);
        EXPECT_TRUE(approx_equal(a.statev, b.statev, "absdiff", 0.));
        EXPECT_TRUE(approx_equal(a.statev_start, b.statev_start, "absdiff", 0.));
    }
    EXPECT_TRUE(approx_equal(a.Etot, b.Etot, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.DEtot, b.DEtot, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.etot, b.etot, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.Detot, b.Detot, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.sigma, b.sigma, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.sigma_start, b.sigma_start, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.PKII, b.PKII, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.PKII_start, b.PKII_start, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.tau, b.tau, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.tau_start, b.tau_start, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.F0, b.F0, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.F1, b.F1, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.U0, b.U0, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.U1, b.U1, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.R, b.R, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.DR, b.DR, "absdiff", 0.));
    EXPECT_EQ(a.T, b.T);
    EXPECT_EQ(a.DT, b.DT);
}

}

// The copy constructor and the assignment operators copy every common field, U0, U1, R, DR,
// tangent_mode, ncp and statev included; copy_fields_M/T copy the same fields as copy_fields.
TEST(Tstate_variables, copies_keep_every_common_field)
{
    state_variables s;
    fill_common(s);
    state_variables s_copy(s);
    expect_common(s_copy, s, true);

    state_variables_M m;
    fill_common(m);
    m.L = 2.*eye(6,6);
    state_variables_M m_copy(m);
    expect_common(m_copy, m, true);
    state_variables_M m_assigned;
    m_assigned = m;
    expect_common(m_assigned, m, true);
    EXPECT_TRUE(approx_equal(m_assigned.L, m.L, "absdiff", 0.));
    state_variables_M m_fields;
    m_fields.copy_fields_M(m);
    expect_common(m_fields, m, false);

    state_variables_T t;
    fill_common(t);
    t.Q = 7.;
    state_variables_T t_copy(t);
    expect_common(t_copy, t, true);
    state_variables_T t_assigned;
    t_assigned = t;
    expect_common(t_assigned, t, true);
    EXPECT_EQ(t_assigned.Q, t.Q);
    state_variables_T t_fields;
    t_fields.copy_fields_T(t);
    expect_common(t_fields, t, false);
}

TEST(Tstate_variables, tensor_members_have_fixed_sizes)
{
    state_variables_M m;
    EXPECT_EQ(m.Etot.n_elem, 6u);
    EXPECT_EQ(m.Lt.n_rows, 6u);
    EXPECT_EQ(m.Wm.n_elem, 4u);
    EXPECT_EQ(m.DR.n_rows, 3u);
    EXPECT_THROW(m.Etot = zeros(7), std::logic_error);
    EXPECT_THROW(m.Lt = zeros(3, 3), std::logic_error);
    EXPECT_THROW(m.Wm = zeros(3), std::logic_error);
    m.statev = zeros(11);   // sized by the law
    EXPECT_EQ(m.statev.n_elem, 11u);
}
