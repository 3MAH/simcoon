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

///@file generation.hpp
///@brief generation for genetic algorithm (among others)
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>
#include "individual.hpp"

namespace simcoon{

/**
 * @file generation.hpp
 * @brief Parameter identification functions.
 */

/** @addtogroup identification
 *  @{
 */


//======================================
class generation
//======================================
{
	private:

	protected:

	public :
    
	/**
	 * @brief Class representing a generation (population) in a genetic algorithm.
	 *
	 * A generation is a collection of individuals (candidate solutions) used in evolutionary
	 * parameter identification. It provides methods for population management, classification,
	 * and unique id assignment.
	 *
	 * @details The class stores:
	 * - Vector of individuals (pop)
	 * - Methods for construction, classification, and destruction
	 */
	class generation
	{
	private:
	protected:
	public:
		std::vector<individual> pop; ///< Population of individuals

		/**
		 * @brief Default constructor.
		 */
		generation();

		/**
		 * @brief Constructor with population size and parameter count.
		 * @param npop Number of individuals
		 * @param nparam Number of parameters per individual
		 * @param id_start Starting id for individuals
		 * @param lambda Initial lambda value (default: 0.0)
		 */
		generation(const int&, const int&, int&, const double & = 0.);

		/**
		 * @brief Copy constructor.
		 * @param gen Generation to copy
		 */
		generation(const generation&);

		/**
		 * @brief Destructor.
		 */
		~generation();

		/**
		 * @brief Get the number of individuals in the population.
		 * @return Population size
		 */
		int size() const {return pop.size();}

		/**
		 * @brief Construct the population with given size and parameters.
		 * @param npop Number of individuals
		 * @param nparam Number of parameters
		 * @param id_start Starting id
		 * @param lambda Initial lambda value
		 */
		void construct(const int&, const int&, int&, const double & = 0.);

		/**
		 * @brief Classify individuals by cost (fitness).
		 */
		void classify();

		/**
		 * @brief Assign new unique ids to individuals.
		 * @param id_start Starting id
		 */
		void newid(int &);

		/**
		 * @brief Destroy the population (clear individuals).
		 */
		void destruct();

		/**
		 * @brief Assignment operator.
		 * @param gen Generation to assign
		 * @return Reference to this object
		 */
		virtual generation& operator = (const generation&);

		/**
		 * @brief Stream output operator.
		 * @param os Output stream
		 * @param gen Generation to output
		 * @return Output stream
		 */
		friend std::ostream& operator << (std::ostream&, const generation&);
	};
		void destruct();

		virtual generation& operator = (const generation&);
		
		friend std::ostream& operator << (std::ostream&, const generation&);
};


/** @} */ // end of identification group

} //namespace simcoon
