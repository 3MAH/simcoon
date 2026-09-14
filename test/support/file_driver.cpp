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

///@file file_driver.cpp
///@brief The historical file-driven solver entry point, kept for the C++ test suite
///@version 1.0

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>
#include <armadillo>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Simulation/Solver/solver_assembly.hpp>
#include <simcoon/Simulation/Solver/solver_sink.hpp>
#include <simcoon/Continuum_mechanics/Micromechanics/multiphase.hpp>

#include "file_readers.hpp"
#include "file_driver.hpp"

using namespace std;
using namespace arma;

namespace simcoon{

namespace {

//Recursive, a sub-phase may itself be a mean-field model (the deleted get_phase_charateristics).
void read_phase_tree(phase_characteristics &phase, const std::string &path_data) {
    const int shape = sub_phase_shape(phase.sptr_matprops->umat_name);
    if (shape == 0) {
        return;
    }
    const std::string nfile = (shape == 2 ? "Nellipsoids" : "Nlayers")
                            + std::to_string(int(phase.sptr_matprops->props(1))) + ".dat";
    if (shape == 2) {
        read_ellipsoid(phase, path_data, nfile);
    }
    else {
        read_layer(phase, path_data, nfile);
    }
    for (auto &sub : phase.sub_phases) {
        read_phase_tree(sub, path_data);
    }
}

} //namespace

void solver(const string &umat_name, const vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const double &div_tnew_dt_solver, const double &mul_tnew_dt_solver, const int &miniter_solver, const int &maxiter_solver, const int &inforce_solver, const double &precision_solver, const double &lambda_solver, const std::string &path_data, const std::string &path_results, const std::string &pathfile, const std::string &outputfile, const int &tangent_mode) {
    if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_algorithmic) {
        throw std::invalid_argument("solver: tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
                                    + std::to_string(tangent_mode) + " (3 = closest-point is reserved)");
    }

    //Check if the required directories exist:
    if(!filesystem::is_directory(path_data)) {
        cout << "error: the folder for the data, " << path_data << ", is not present" << endl;
        return;
    }
    if(!filesystem::is_directory(path_results)) {
        cout << "The folder for the results, " << path_results << ", is not present and has been created" << endl;
        filesystem::create_directory(path_results);
    }

    std::string ext_filename = outputfile.substr(outputfile.length()-4,outputfile.length());
    std::string filename = outputfile.substr(0,outputfile.length()-4); //to remove the extension

    std::string outputfile_global = filename + "_global" + ext_filename;
    std::string outputfile_local = filename + "_local" + ext_filename;

    std::string output_info_file = "output.dat";

    std::vector<block> blocks;  //loading blocks
    double T_init = 0.;

    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);

    solver_output so(blocks.size());
    read_output(so, blocks.size(), nstatev, path_data, output_info_file);

    //Check output and step files
    check_path_output(blocks, so);

    solver_params ctrl;
    ctrl.div_tnew_dt = div_tnew_dt_solver;
    ctrl.mul_tnew_dt = mul_tnew_dt_solver;
    ctrl.miniter = miniter_solver;
    ctrl.maxiter = maxiter_solver;
    ctrl.inforce = inforce_solver;
    ctrl.precision = precision_solver;
    ctrl.lambda = lambda_solver;
    ctrl.tangent_mode = tangent_mode;

    //The phase tree is read here, with the rest of the file semantics: solver_run takes its
    //sub-phases in memory.
    std::vector<phase_characteristics> sub_phases;
    {
        phase_characteristics rve_phases;
        rve_phases.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
        rve_phases.construct(0, 1);
        rve_phases.sptr_sv_global->T = T_init;   //the sub-phases start at the path's initial temperature
        read_phase_tree(rve_phases, path_data);
        sub_phases = std::move(rve_phases.sub_phases);
    }

    solver_file_sink sink(path_results, outputfile_global, outputfile_local);
    //status intentionally ignored: the historical file-driven solver() returned void on early aborts
    solver_run(blocks, T_init, so, umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, corate_type, ctrl, sink, sub_phases);
}

} //namespace simcoon
