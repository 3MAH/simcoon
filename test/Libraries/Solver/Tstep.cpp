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


///@file Tstep.cpp
///@brief Copies of step_meca / step_thermomeca reproduce every field (memberwise, defaulted)
///@version 1.0

#include <gtest/gtest.h>
#include <armadillo>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>

using namespace arma;
using namespace simcoon;

namespace {

template <typename S>
void fill_step(S &s) {
    s.number = 3; s.Dn_init = 0.5; s.Dn_mini = 0.01; s.Dn_inc = 0.2; s.ninc = 4; s.mode = 3;
    s.control_type = 3; s.times = linspace(0., 1., 4); s.BC_Time = 1.; s.tab_data = randu(4, 7);
    s.cBC_meca = Col<int>({0, 1, 1, 1, 1, 1}); s.BC_meca = randu(6); s.mecas = randu(4, 6);
    s.BC_mecas = randu(4, 6); s.BC_w = randu(3, 3); s.BC_R = randu(3, 3); s.BC_T = 320.;
    s.cBC_T = 1; s.Ts = randu(4); s.BC_Ts = randu(4);
}

template <typename S>
void expect_same(const S &a, const S &b) {
    EXPECT_EQ(a.number, b.number); EXPECT_EQ(a.Dn_init, b.Dn_init); EXPECT_EQ(a.Dn_mini, b.Dn_mini);
    EXPECT_EQ(a.Dn_inc, b.Dn_inc); EXPECT_EQ(a.ninc, b.ninc); EXPECT_EQ(a.mode, b.mode);
    EXPECT_EQ(a.control_type, b.control_type); EXPECT_EQ(a.BC_Time, b.BC_Time);
    EXPECT_TRUE(approx_equal(a.times, b.times, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.tab_data, b.tab_data, "absdiff", 0.));
    EXPECT_TRUE(all(a.cBC_meca == b.cBC_meca));
    EXPECT_TRUE(approx_equal(a.BC_meca, b.BC_meca, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.mecas, b.mecas, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.BC_mecas, b.BC_mecas, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.BC_w, b.BC_w, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.BC_R, b.BC_R, "absdiff", 0.));
    EXPECT_EQ(a.BC_T, b.BC_T); EXPECT_EQ(a.cBC_T, b.cBC_T);
    EXPECT_TRUE(approx_equal(a.Ts, b.Ts, "absdiff", 0.));
    EXPECT_TRUE(approx_equal(a.BC_Ts, b.BC_Ts, "absdiff", 0.));
}

}  // namespace

TEST(Tstep, step_meca_copy_and_assignment)
{
    step_meca s;
    fill_step(s);
    const step_meca copied(s);
    expect_same(s, copied);
    step_meca assigned;
    assigned = s;
    expect_same(s, assigned);
}

TEST(Tstep, step_thermomeca_copy_and_assignment)
{
    step_thermomeca s;
    fill_step(s);
    const step_thermomeca copied(s);
    expect_same(s, copied);
    step_thermomeca assigned;
    assigned = s;
    expect_same(s, assigned);
}
