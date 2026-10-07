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

///@file solver_sink.cpp
///@brief Results sinks for solver_run (see solver_sink.hpp)

#include <string>
#include <memory>
#include <armadillo>
#include <simcoon/Simulation/Solver/solver_sink.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

solver_file_sink::solver_file_sink(const std::string &mpath_results, const std::string &moutputfile_global, const std::string &moutputfile_local)
    : path_results(mpath_results), outputfile_global(moutputfile_global), outputfile_local(moutputfile_local) {}

void solver_file_sink::init(phase_characteristics &rve) {
    rve.define_output(path_results, outputfile_global, "global");
    rve.define_output(path_results, outputfile_local, "local");
}

void solver_file_sink::record(phase_characteristics &rve, const solver_output &so, const int &kblock, const int &kcycle, const int &kstep, const int &kinc, const double &Time) {
    rve.output(so, kblock, kcycle, kstep, kinc, Time, "global");
    rve.output(so, kblock, kcycle, kstep, kinc, Time, "local");
}

void solver_memory_sink::init(phase_characteristics &rve) {
    if (std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global) != nullptr)
        sv_type = 2;
    else if (std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global) != nullptr)
        sv_type = 1;
}

void solver_memory_sink::record(phase_characteristics &rve, const solver_output &so, const int &kblock, const int &kcycle, const int &kstep, const int &kinc, const double &Time) {
    (void)so;
    auto sv = rve.sptr_sv_global;
    // defensive re-detection of the state-variable type (runs are single-type:
    // the solver_run binding rejects mixed mechanical/thermomechanical blocks)
    auto sv_T = std::dynamic_pointer_cast<state_variables_T>(sv);
    std::shared_ptr<state_variables_M> sv_M;
    if (sv_T == nullptr) {
        sv_M = std::dynamic_pointer_cast<state_variables_M>(sv);
    }
    sv_type = (sv_T != nullptr) ? 2 : ((sv_M != nullptr) ? 1 : 0);

    blocks_i.push_back(kblock);
    cycles_i.push_back(kcycle);
    steps_i.push_back(kstep);
    incs_i.push_back(kinc);
    time.push_back(Time);

    Etot.push_back(sv->Etot);
    etot.push_back(sv->etot);
    PKII.push_back(sv->PKII);
    tau.push_back(sv->tau);
    // consistent with phase_characteristics::output (o_stress_type 4); the thermomechanical
    // route sets tau = sigma (F = I) in select_umat_T, so this holds for both block types
    sigma.push_back(sv->tau / det(sv->F1));
    F1.push_back(sv->F1);
    R.push_back(sv->R);
    DR.push_back(sv->DR);
    T.push_back(sv->T);
    statev.push_back(sv->statev);

    if (sv_T != nullptr) {
        Wm.push_back(sv_T->Wm);
        Wt.push_back(sv_T->Wt);
        Q.push_back(sv_T->Q);
        r.push_back(sv_T->r);
        if (record_tangent) {
            dSdE.push_back(sv_T->dSdE);
            dSdT.push_back(sv_T->dSdT);
            drdE.push_back(sv_T->drdE);
            drdT.push_back(sv_T->drdT);
        }
    }
    else if (sv_M != nullptr) {
        Wm.push_back(sv_M->Wm);
        if (record_tangent) {
            Lt.push_back(sv_M->Lt);
        }
    }
}

} //namespace simcoon
