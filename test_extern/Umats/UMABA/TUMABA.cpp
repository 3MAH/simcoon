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

///@file TUMABA.cpp
///@brief The UMABA plugin driven by the solver: the historical reference case, its
///loading programme built in code (nothing is read from a file since the 2.0 JSON-only
///migration; the plugin itself is what the test exercises).
///@version 1.0

#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Solver/solver_sink.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

TEST(TUMABA, TUMABA_solver)
{
    //The material and the programme of the reference case: an elastic plugin (E = 70000 MPa,
    //nu = 0.3), loaded in uniaxial tension to 2 % strain and unloaded, 100 increments each
    const vec props = {70000., 0.3, 0.};
    const double T_init = 290.;
    const double targets[2] = {0.02, 0.};

    std::vector<block> blocks(1);
    blocks[0].number = 1;
    blocks[0].type = 1;
    blocks[0].control_type = 1;
    blocks[0].ncycle = 1;
    blocks[0].nstep = 2;
    blocks[0].generate();
    for (unsigned int j = 0; j < 2; j++) {
        auto s = std::dynamic_pointer_cast<step_meca>(blocks[0].steps[j]);
        s->number = j + 1;
        s->control_type = 1;
        s->mode = 1;
        s->Dn_init = 1.;
        s->Dn_mini = 1.;
        s->Dn_inc = 0.01;
        s->BC_Time = 1.;
        s->cBC_meca(0) = 0;                 //strain-driven 11
        for (unsigned int k = 1; k < 6; k++) {
            s->cBC_meca(k) = 1;             //stress-free otherwise
        }
        s->BC_meca = zeros(6);
        s->BC_meca(0) = targets[j];
        s->cBC_T = 0;
        s->BC_T = T_init;
    }

    solver_output so(1);
    so.o_type(0) = 1;
    so.o_nfreq(0) = 1;
    solver_params ctrl;
    solver_memory_sink sink;
    sink.record_tangent = false;

    const int status = solver_run(blocks, T_init, so, "UMABA", props, 1, 0., 0., 0., 0, 2, ctrl, sink);
    ASSERT_EQ(status, 0);

    //The committed reference: strain in columns 8:14, stress in 14:20
    mat C;
    ASSERT_TRUE(C.load("comparison/results_job_global-0.txt")) << "missing reference";
    ASSERT_EQ(C.n_rows, sink.sigma.size());
    for (unsigned int i = 0; i < C.n_rows; i++) {
        for (unsigned int k = 0; k < 6; k++) {
            EXPECT_LT(fabs(C(i, 8 + k) - sink.Etot[i](k)), 1.E-6);
            EXPECT_LT(fabs(C(i, 14 + k) - sink.sigma[i](k)), 1.E-6);
        }
    }
}
