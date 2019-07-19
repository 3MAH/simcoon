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

///@file L_eff.cpp
///@brief solver: Determination of the effective elastic properties	of a composite
///@version 1.9

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

int main()
{

  string umat_name;
  string path_data = "data";
  string path_results = "results";
  string pathfile = "path.txt";
  string materialfile = "material.dat";
  string output_info_file = "output.dat";
  string resultslfile = "results.dat";
  string inputfile;

  string outputfile = "results_aba.txt";
  string ext_filename = outputfile.substr(outputfile.length()-4,outputfile.length());
  string filename = outputfile.substr(0,outputfile.length()-4); //to remove the extension
  string outputfile_global = filename + "_global" + ext_filename;

  unsigned int nprops = 0;
  unsigned int nstatev = 0;

  double psi_rve = 0.;
  double theta_rve = 0.;
  double phi_rve = 0.;
  vec props, E, S;

  double T_init = 273.15;
  std::vector<block> blocks;

  solver_output so(blocks.size());
  read_output(so, blocks.size(), nstatev, path_data, output_info_file);

  read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
  read_path(blocks, T_init, path_data, pathfile);
  phase_characteristics rve;

  rve.construct(0,1);
  rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);
  rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(3,3), zeros(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev));

  /*
  inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
  read_layer(rve, path_data, inputfile);
  */

  inputfile = "Nphases" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
  read_phase(rve, path_data, inputfile);

  rve.define_output(path_results, outputfile_global, "global");

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    for(unsigned int n = 0 ; n < blocks[i].ncycle ; n++){
      if (i != 0){
        for(unsigned int j = 0 ; j < blocks[i].nstep ; j++){
            read_results(E,S,path_data,"results-"+to_string(j+blocks[i-1].nstep)+".dat");
            auto sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
            sv_M->Etot = E;
            sv_M->sigma = S;
            for (unsigned int k=0; k<rve.sub_phases.size(); k++) {
                read_subphases_results(E,S,k,path_data,"results-"+to_string(j+blocks[i-1].nstep)+".dat");
                auto sv_M_k = std::dynamic_pointer_cast<state_variables_M>(rve.sub_phases[k].sptr_sv_global);
                sv_M_k->Etot = E;
                sv_M_k->sigma = S;
            }
            double inc = round(1./blocks[i].steps[j]->Dn_inc);
            rve.output(so, i, n, j, inc-1, 1, "global");
            }
      }
      else {
        for(unsigned int j = 0 ; j < blocks[i].nstep ; j++){
            read_results(E,S,path_data,"results-"+to_string(j)+".dat");
            auto sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
            sv_M->Etot = E;
            sv_M->sigma = S;
            for (unsigned int k=0; k<rve.sub_phases.size(); k++) {
                read_subphases_results(E,S,k,path_data,"results-"+to_string(j)+".dat");
                auto sv_M_k = std::dynamic_pointer_cast<state_variables_M>(rve.sub_phases[k].sptr_sv_global);
                sv_M_k->Etot = E;
                sv_M_k->sigma = S;
            }
            double inc = round(1./blocks[i].steps[j]->Dn_inc);
            rve.output(so, i, n, j, inc-1, 1., "global");
            }
      }

        }
      }

  return 0;
}