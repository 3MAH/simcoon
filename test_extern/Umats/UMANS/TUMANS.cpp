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

///@file TUMANS.cpp
///@brief The UMANS plugin's data layout against a USERMAT that exposes what it receives
///       (usermat_layout.cpp): component order, column-major tangent and rotation, the law's
///       statev without header, energies, condensed (plane stress) components
///@version 1.0

#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <armadillo>
#include <dylib.hpp>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_plugin_api.hpp>

using namespace arma;
using namespace simcoon;

namespace {

struct ans_plugin {
    dylib::library lib{"external/umat_plugin_ans", dylib::decorations::os_default()};
    std::unique_ptr<umat_plugin_ans_api, void(*)(umat_plugin_ans_api*)> api{
        lib.get_function<umat_plugin_ans_api*()>("create_api")(),
        lib.get_function<void(umat_plugin_ans_api*)>("destroy_api")};
};

// one call of the plugin on a fresh point with nstatev law variables; returns the local state
std::shared_ptr<state_variables_M> call(ans_plugin &p, const mat &DR, const int &ndi, const int &nshr) {
    const int nstatev = 10;
    phase_characteristics rve;
    rve.construct(0, 1);
    rve.sptr_matprops->update(0, "UMANS", 1, 0., 0., 0., 3, vec{70000., 0.3, 0.});
    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6),
                               zeros(6), zeros(6), eye(3, 3), eye(3, 3), eye(3, 3), eye(3, 3), eye(3, 3), eye(3, 3), 290., 0.,
                               nstatev, zeros(nstatev), zeros(nstatev));
    rve.global2local();
    auto sv = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_local);
    sv->DEtot = vec{1e-3, 0., 0., 0., 0., 0.};
    bool start = false;
    double tnew_dt = 1.;
    p.api->umat_ansys(rve, DR, 0., 1., ndi, nshr, start, 0, tnew_dt);
    return sv;
}

}  // namespace

TEST(TUMANS, data_layout_3D)
{
    ans_plugin p;
    const mat DR = Rotation::from_euler(0.3, 0.5, 0.7, "zxz").as_matrix();
    auto sv = call(p, DR, 3, 3);
    const int K[6] = {0, 1, 2, 3, 5, 4};   // Ansys 11, 22, 33, 12, 23, 13
    for (int i = 0; i < 6; i++) {
        EXPECT_EQ(sv->sigma(K[i]), 10.*(i + 1));
        for (int j = 0; j < 6; j++) EXPECT_EQ(sv->Lt(K[i], K[j]), 100.*(i + 1) + (j + 1)) << i << "," << j;
    }
    for (int k = 0; k < 9; k++) EXPECT_EQ(sv->statev(k), DR(k % 3, k / 3));   // rotateM(r,c) = DR(r,c)
    EXPECT_EQ(sv->statev(9), 10.);       // the law's own count, no simcoon header
    EXPECT_EQ(sv->Wm(1), 1.5);
    EXPECT_EQ(sv->Wm(3), 2.5);
}

TEST(TUMANS, data_layout_plane_stress)
{
    ans_plugin p;
    auto sv = call(p, eye(3, 3), 2, 1);
    const int K[3] = {0, 1, 3};          // Ansys plane stress 11, 22, 12
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) EXPECT_EQ(sv->Lt(K[i], K[j]), 100.*(i + 1) + (j + 1));
    EXPECT_EQ(sv->Lt(2, 2), 0.);         // absent components stay zero
}
