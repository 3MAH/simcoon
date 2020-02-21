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
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Simulation/Solver/block.hpp>

namespace simcoon{
    
/// Function that reads the output parameters
void update_sections(section_characteristics &, const int &, const int &, const int &, const std::string & = "data");
    
void write_section(section_characteristics &, const unsigned int &, const std::string & = "data", const std::string & = "Nmat_0.inp");

void write_sections(std::vector<section_characteristics> &, const unsigned int &, const std::string & = "data", const std::string & = "Nmat_0.inp");
    
void write_sections(section_characteristics &, const unsigned int &,const std::string & = "data", const std::string & = "Nmat_0.inp");

void update_materials(std::vector<aba_material> &, const int &, const int &, const int &, const std::string & = "data");

void write_materials(std::vector<aba_material> &, const unsigned &, const std::string & = "data", const std::string & = "Nmat_0.inp");

void update_steps(std::vector<std::shared_ptr<step> > &, const std::vector<block> &, const bool &, const int &, const int &);

void write_steps(const std::vector<std::shared_ptr<step> > &, const int &, const double &, const std::string & = "data", const std::string & = "Nstep_0.inp");
    
void write_nodes_file(const std::vector<Node> &, const std::string &, const std::string &);
void write_elements_file(const std::vector<Element> &, const std::string &, const std::string &);
void write_sets_file(const std::vector<section_characteristics> &, const std::vector<Node> &, const std::string &, const std::string &);

void write_node_set(const std::string &, const Node &, std::ofstream &);
void write_nodes_set(const std::string &, const std::vector<Node> &, std::ofstream &);

void write_element_set(const std::string &, const Element &, std::ofstream &);
void write_elements_set(const std::string &, const std::vector<Element> &, std::ofstream &);

void append_perio_nodes(const cubic_mesh &, const std::string &, const std::string &);
    
void write_PBC(const cubic_mesh &, const std::string & = "data", const std::string & = "PBC_0.inp");

void write_TIE(const cubic_mesh &, const cubic_mesh &, const std::string &, const std::string &);
    
void write_CDN(const cubic_mesh &, const std::string & = "data", const std::string & = "CDN_0.inp");

void write_eq(std::ostream &, const equation &);
    
void write_NonPerio_CDN(const cubic_mesh &, const cubic_mesh &, const unsigned int &, const std::string & = "data", const std::string & = "CDN_0.inp");
    
void write_run_perturbation_file(const std::string &, const std::string &);
} //namespace simcoon
