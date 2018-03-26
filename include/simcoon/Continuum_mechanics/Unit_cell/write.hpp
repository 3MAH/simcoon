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

///@file write
///@brief To write the necessary Abaqus input files
///@version 1.0

#pragma once
#include <armadillo>
#include <string>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/steps.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Simulation/Solver/block.hpp>

namespace simcoon{
    
/// Function that reads the output parameters
void update_sections(section_characteristics &, const int &, const int &, const int &, const std::string & = "data");
    
void write_section(section_characteristics &, const std::string & = "data", const std::string & = "Nmat_0.inp");

void write_sections(std::vector<section_characteristics> &, const std::string & = "data", const std::string & = "Nmat_0.inp");
    
void write_sections(section_characteristics &, const std::string & = "data", const std::string & = "Nmat_0.inp");    

void update_materials(std::vector<aba_material> &, const int &, const int &, const int &, const std::string & = "data");

void write_materials(std::vector<aba_material> &, const std::string & = "data", const std::string & = "Nmat_0.inp");

void update_steps(std::vector<aba_step_meca> &, const std::vector<block> &, const bool & = false);

void write_steps(std::vector<aba_step_meca> &, const double &, const std::string & = "data", const std::string & = "Nstep_0.inp");
    
void write_nodes_file(const Node &, std::ofstream &);
void write_node_set(const std::string &, const Node &, std::ofstream &);
void write_nodes_set(const std::string &, const std::vector<Node> &, std::ofstream &);

void append_perio_nodes(const cubic_mesh &, const std::string &, const std::string &);
    
void write_PBC(const cubic_mesh &, const unsigned int &, const std::string & = "data", const std::string & = "PBC_0.inp");

void write_TIE(const cubic_mesh &, const cubic_mesh &, const std::string &, const std::string &);
    
void write_CDN(const cubic_mesh &, const std::string & = "data", const std::string & = "CDN_0.inp");

void write_eq(std::ostream &, const equation &);
    
void write_NonPerio2_CDN(const cubic_mesh &, const cubic_mesh &, const unsigned int &, const std::string & = "data", const std::string & = "CDN_0.inp");
    
//void write_NonPerio_CDN(const cubic_mesh &, const cubic_mesh &, const std::string & = "data", const std::string & = "CDN_0.inp");
    
} //namespace simcoon
