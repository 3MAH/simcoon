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
#include <iostream>
#include <armadillo>

namespace simcoon{

/**
 * @file individual.hpp
 * @brief Parameter identification functions.
 */

/** @addtogroup identification
 *  @{
 */



/**
 * @brief Class representing an individual in a genetic algorithm for parameter identification.
 *
 * An individual encodes a candidate solution (set of parameters) and its associated cost (fitness).
 * Used in evolutionary algorithms for model calibration or optimization.
 *
 * @details The class stores:
 * - Parameter vector (p)
 * - Cost function value (cout)
 * - Unique identifier (id) and rank in the population
 * - Lambda: step size or regularization parameter
 */
class individual
{
private:
protected:
public:
	int np; ///< Number of parameters
	double cout; ///< Cost function value (fitness)
	int id; ///< Unique identifier
	int rank; ///< Rank in the population
	arma::vec p; ///< Parameter vector
	double lambda; ///< Step or regularization parameter

	/**
	 * @brief Default constructor.
	 */
	individual();

	/**
	 * @brief Constructor with number of parameters, id, and lambda.
	 * @param np Number of parameters
	 * @param id Unique identifier
	 * @param lambda Step or regularization parameter
	 */
	individual(const int &np, const int &id, const double &lambda);

	/**
	 * @brief Copy constructor.
	 * @param ind Individual to copy
	 */
	individual(const individual &ind);

	/**
	 * @brief Destructor.
	 */
	~individual();

	/**
	 * @brief Allocate and initialize parameter vector.
	 */
	void construct();

	/**
	 * @brief Assignment operator.
	 * @param ind Individual to assign
	 * @return Reference to this object
	 */
	virtual individual& operator = (const individual &ind);

	/**
	 * @brief Stream output operator.
	 * @param os Output stream
	 * @param ind Individual to output
	 * @return Output stream
	 */
	friend std::ostream& operator << (std::ostream& os, const individual &ind);
};
    

/** @} */ // end of identification group

} //namespace simcoon
