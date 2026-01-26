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

///@file individual.hpp
///@brief individual for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <fstream>
#include <armadillo>
#include "generation.hpp"

namespace simcoon{

/**
 * @file methods.hpp
 * @brief Parameter identification functions.
 */

/** @addtogroup identification
 *  @{
 */

    
//Genetic method
void genetic(generation &, generation &, int &, const double &, const double &, const std::vector<parameters> &);

///Genrun creation
void to_run(generation &, generation &, generation &, const double &, const std::vector<parameters> &);

//Find the bests from the gensons and previous generation, considering the gboys
//Define the new gen_cur and gboys_cur accordingly
void find_best(generation &, generation &, const generation &, const generation &, const generation &, const int &, const int &, int &);

//Write the results in an output file
void write_results(std::ofstream &, const std::string &outputfile, const generation &, const int &, const int &, const int &);


/** @} */ // end of identification group

} //namespace simcoon