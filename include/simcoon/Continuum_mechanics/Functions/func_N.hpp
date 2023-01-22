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

/**
* @file func_N.hpp
* @author Yves Chemisky 
* @section Functions that computes cumulative  laws
*/

#include <armadillo>
#include <string>

namespace simcoon{

/**
 * @brief Read the parameters and variables for a user-defined parametric function
 * 
 * This function is utilized to use the simcoon optimization algorithm to identifiy the parameters of a parametric function
 * 
 * @warning : This function should be deprecated in a next release of simcoon
 * 
 * @param[out] params (arma::vec) vector of parameters
 * @param[out] variables (arma::vec) vector of variables
 * @param[in] path_data (std::string) path to find the file to be read 
 * @param[in] inputfile (std::string) name of the input file to be read
 * 
*/
void read_func_N(arma::vec&params, arma::vec&variables, const std::string &path_data, const std::string &inputfile);

/**
 * @brief Compute the result a user-defined parametric function
 * 
 *  This function computes the result of a single function \f$ \underline{y} = f(\underline{x}) \f$, providing the vector of input values \f$ \underline{x} \f$
 *  The list of values shall be stored in a file named \b inputfile, in the data folder \b path_data.
 * 
 * @note This is a temproary function and necessitate to modify the code with your own function. You shall define your own function along with the definition of the vector, for example by adding:
 *  @code
 *      vec y = p_cumulative(N, variables(0), variables(1), params);
 *      Insert here the fonction you want
 *  @endcode
 * 
 * 
 * in the file func_N.cpp. You should then reinstall the library
 * 
 * @warning This function should be deprecated in a next release of simcoon
 * 
 * The x and y values are written in a file named \b outputfile, in the data folder \b path_results
 * 
 * @param[out] params (arma::vec) vector of parameters
 * @param[out] variables (arma::vec) vector of variables
 * @param[in] inputfile (std::string) name of the input file to be read
 * @param[in] outputfile Name of the output file
 * @param[in] path_data The path of the data folder
 * @param[in] path_results The path of the results folder
 * 
 * @details Example
 *  @code
 *      string outputfile = "results.txt";
 *      string inputfile = "list_inputs.txt";
 *      string path_data = "data";
 *      string path_results = "results";
 *
 *      vec props = {1., 2.};
 *      vec sigma = randu(6);
 *      double sigma_eq = Mises_stress(sigma);
 *      vec variables = {sigma_eq};
 *
 *      func_N(props, variables, inputfile, outputfile, path_data, path_results);
 *  @endcode
 * 
*/
void func_N(const arma::vec &params, const arma::vec &variables, const std::string &inputfile, const std::string &outputfile, const std::string &path_data, const std::string &path_results);

} //namespace simcoon
