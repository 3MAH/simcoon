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

///@file Tmaterial_characteristics.cpp
///@brief Test for the material frame of a phase
///@version 1.0

#include <gtest/gtest.h>
#include <algorithm>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

// The frame follows the angles through every writer, is copied with the object, and refuses
// to serve a frame that no longer matches directly assigned angles.
TEST(Tmaterial_characteristics, frame_follows_the_angle_writers)
{
    vec props = {70000., 0.3, 1.e-5};
    material_characteristics m(1, "ELISO", 1, 0.5, 0.2, -0.1, 3, props);

    frame_rotation ref(Rotation::from_euler(0.5, 0.2, -0.1, "zxz"));
    vec v = {100., -40., 25., 30., -12., 8.};
    vec a = v, b = v;
    m.frame().rotate_stress(a);
    ref.rotate_stress(b);
    EXPECT_FALSE(m.frame().is_identity());
    EXPECT_TRUE(std::equal(a.begin(), a.end(), b.begin()));

    // update() rebuilds
    m.update(1, "ELISO", 1, 0., 0., 0., 3, props);
    EXPECT_TRUE(m.frame().is_identity());

    // copy constructor and assignment carry the frame
    m.update(1, "ELISO", 1, 0.5, 0.2, -0.1, 3, props);
    material_characteristics c(m);
    material_characteristics d;
    d = m;
    vec ca = v, da = v;
    c.frame().rotate_stress(ca);
    d.frame().rotate_stress(da);
    EXPECT_TRUE(std::equal(ca.begin(), ca.end(), b.begin()));
    EXPECT_TRUE(std::equal(da.begin(), da.end(), b.begin()));

    // a direct assignment to the public angles is detected instead of serving a stale frame
    m.psi_mat = 1.0;
    EXPECT_THROW(m.frame(), std::logic_error);

    // the default object is unoriented
    EXPECT_TRUE(material_characteristics().frame().is_identity());
}
