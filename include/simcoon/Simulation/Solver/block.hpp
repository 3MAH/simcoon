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

///@file block.hpp
///@brief object that defines a block
///@version 1.0

#pragma once
#include <iostream>
#include <memory>
#include "step.hpp"

namespace simcoon{

/**
 * @file block.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


/**
 * @brief Class representing a loading block containing multiple steps.
 * 
 * A block groups several loading steps that can be repeated (cycled) multiple times.
 * This is useful for simulating cyclic loading, fatigue tests, or repeated thermomechanical cycles.
 * 
 * @details The block structure allows organizing complex loading histories:
 * - Multiple steps within a block can represent different loading phases
 * - The block can be cycled to simulate repetitive loading
 * - Different control types (stress, strain, mixed) can be specified
 */
class block
{
	private:

	protected:

	public :
		unsigned int number; ///< Block identification number
        unsigned int nstep; ///< Number of steps in the block
        unsigned int ncycle; ///< Number of cycles to repeat the block
        unsigned int type; ///< Type of block (loading type identifier)
        unsigned int control_type; ///< Control type (0: strain, 1: stress, 2: mixed)
    
        std::vector<std::shared_ptr<step> > steps; ///< Vector of step pointers
    
        /**
         * @brief Default constructor.
         */
        block();
        
        /**
         * @brief Constructor with basic parameters.
         * @param number Block identification number
         * @param nstep Number of steps
         * @param ncycle Number of cycles
         * @param type Block type
         * @param control_type Control type
         */
        block(const unsigned int &number, const unsigned int &nstep, const unsigned int &ncycle, const unsigned int &type, const unsigned int &control_type);
        
        /**
         * @brief Full constructor with step vector.
         * @param number Block identification number
         * @param nstep Number of steps
         * @param ncycle Number of cycles
         * @param type Block type
         * @param control_type Control type
         * @param steps Vector of step pointers
         */
        block(const unsigned int &number, const unsigned int &nstep, const unsigned int &ncycle, const unsigned int &type, const unsigned int &control_type, const std::vector<std::shared_ptr<step> > &steps);
        
        /**
         * @brief Copy constructor.
         * @param b Block to copy
         */
        block(const block &b);
        
        /**
         * @brief Virtual destructor.
         */
		virtual ~block();
		
        /**
         * @brief Generate the steps within the block.
         */
		void generate();
    
        /**
         * @brief Assignment operator.
         * @param b Block to assign
         * @return Reference to this object
         */
        virtual block& operator = (const block &b);
		
        /**
         * @brief Stream output operator.
         * @param os Output stream
         * @param b Block to output
         * @return Output stream
         */
        friend std::ostream& operator << (std::ostream& os, const block &b);
};


/** @} */ // end of solver group

} //namespace simcoon
