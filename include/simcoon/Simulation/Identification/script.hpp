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

///@file script.hpp
///@brief Scripts that allows to run identification algorithms based on Smart+ Control functions
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "constants.hpp"
#include "parameters.hpp"
#include "opti_data.hpp"
#include "individual.hpp"
#include "generation.hpp"

namespace simcoon{

/**
 * @file script.hpp
 * @brief Parameter identification functions.
 */

/** @addtogroup identification
 *  @{
 */


//This function will copy the parameters files
void copy_parameters(const std::vector<parameters> &, const std::string &, const std::string &);

//This function will copy the parameters files
void copy_constants(const std::vector<constants> &, const std::string &, const std::string &);
    
//This function will replace the keys by the parameters
void apply_parameters(const std::vector<parameters> &, const std::string &);

//This function will replace the keys by the parameters
void apply_constants(const std::vector<constants> &, const std::string &);
    
//Read the control parameters of the optimization algorithm
void launch_solver(const generation &, const int &, std::vector<parameters> &, std::vector<constants> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);
    
//Read the control parameters of the optimization algorithm
void launch_odf(const generation &, std::vector<parameters> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);
    
//Read the control parameters of the optimization algorithm
    void launch_func_N(const generation &, const int &, std::vector<parameters> &, std::vector<constants> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);
    
void run_simulation(const std::string &, const individual &, const int &, std::vector<parameters> &, std::vector<constants> &, std::vector<opti_data> &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string&);
    
double calc_cost(const arma::vec &, arma::vec &, const arma::vec &, const std::vector<opti_data> &, const std::vector<opti_data> &, const int &, const int &);

arma::mat calc_sensi(const individual &, generation &, const std::string &, const int &, const int &, std::vector<parameters> &, std::vector<constants> &, arma::vec &, std::vector<opti_data> &, std::vector<opti_data> &, const std::string &, const std::string &, const std::string &, const std::string &, const int &, const arma::vec &, const std::string&);

    

/** @} */ // end of identification group

} //namespace simcoon
