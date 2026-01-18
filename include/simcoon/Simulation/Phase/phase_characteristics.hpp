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

///@file phase_characteristics.hpp
///@brief Characteristics of a phase, which hereditates from:
///-material characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <armadillo>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Solver/output.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file phase_characteristics.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */


/**
 * @brief Class representing a material phase with full mechanical state.
 * 
 * This class combines geometry, material properties, and state variables to represent
 * a complete phase in a heterogeneous material. It supports hierarchical structures
 * for multi-scale modeling, where a phase can contain sub-phases.
 * 
 * @details The phase_characteristics class manages:
 * - Geometric description (shape type and parameters)
 * - Material properties (UMAT name, elastic constants, etc.)
 * - State variables in both local and global coordinate systems
 * - Output streams for results
 * - Sub-phases for multi-scale homogenization
 */
class phase_characteristics
{
	private:

	protected:

	public :

        int shape_type; ///< Type of geometric shape (0: sphere, 1: cylinder, 2: ellipsoid, etc.)
        int sv_type; ///< Type of state variables (0: mechanical, 1: thermomechanical)
        std::shared_ptr<geometry> sptr_shape; ///< Pointer to the geometric shape of the phase
        std::shared_ptr<phase_multi> sptr_multi; ///< Pointer to multiscale information
        std::shared_ptr<material_characteristics> sptr_matprops; ///< Pointer to material properties
        std::shared_ptr<state_variables> sptr_sv_global; ///< Pointer to state variables in global frame
        std::shared_ptr<state_variables> sptr_sv_local; ///< Pointer to state variables in local frame
        std::shared_ptr<std::ofstream> sptr_out_global; ///< Output stream for global results
        std::shared_ptr<std::ofstream> sptr_out_local; ///< Output stream for local results
    
        std::vector<phase_characteristics> sub_phases; ///< Vector of sub-phases for multi-scale modeling
        std::string sub_phases_file; ///< Filename containing sub-phase definitions
    
        /**
         * @brief Default constructor.
         */
		phase_characteristics();
    
        /**
         * @brief Full constructor with all parameters.
         * @param shape_type Type of geometric shape
         * @param sv_type Type of state variables
         * @param sptr_shape Pointer to geometry
         * @param sptr_multi Pointer to multiscale information
         * @param sptr_matprops Pointer to material properties
         * @param sptr_sv_global Pointer to global state variables
         * @param sptr_sv_local Pointer to local state variables
         * @param sptr_out_global Global output stream
         * @param sptr_out_local Local output stream
         * @param sub_phases_file File with sub-phase definitions
         */
        phase_characteristics(const int &shape_type, const int &sv_type, const std::shared_ptr<geometry> &sptr_shape, const std::shared_ptr<phase_multi> &sptr_multi, const std::shared_ptr<material_characteristics> &sptr_matprops, const std::shared_ptr<state_variables> &sptr_sv_global, const std::shared_ptr<state_variables> &sptr_sv_local, const std::shared_ptr<std::ofstream> &sptr_out_global, const std::shared_ptr<std::ofstream> &sptr_out_local, const std::string &sub_phases_file);

        /**
         * @brief Copy constructor.
         * @param pc Phase characteristics to copy
         */
        phase_characteristics(const phase_characteristics &pc);
        
        /**
         * @brief Virtual destructor.
         */
        virtual ~phase_characteristics();
    
        /**
         * @brief Construct the phase with given types.
         * @param shape_type Type of geometric shape
         * @param sv_type Type of state variables
         */
        virtual void construct(const int &shape_type, const int &sv_type);
        
        /**
         * @brief Construct sub-phases recursively.
         * @param shape_type Type of geometric shape
         * @param sv_type Type of state variables
         * @param n_sub Number of sub-phases
         */
        virtual void sub_phases_construct(const int &shape_type, const int &sv_type, const int &n_sub);
        
        /**
         * @brief Copy current values to start-of-increment values.
         */
        virtual void to_start();
        
        /**
         * @brief Set current values from start-of-increment values.
         * @param control Control flag for selective update
         */
        virtual void set_start(const int &control);
        
        /**
         * @brief Transform state variables from local to global frame.
         */
        virtual void local2global();
        
        /**
         * @brief Transform state variables from global to local frame.
         */
        virtual void global2local();
        
        /**
         * @brief Copy from another phase_characteristics (without ofstreams).
         * @param pc Phase characteristics to copy
         * @warning Output streams are NOT copied
         */
        virtual void copy(const phase_characteristics& pc);

        /**
         * @brief Assignment operator.
         * @param pc Phase characteristics to assign
         * @return Reference to this object
         */
        virtual phase_characteristics& operator = (const phase_characteristics& pc);
    
        /**
         * @brief Define output file for results.
         * @param path Output directory path
         * @param filename Base filename (default: "results")
         * @param outputtype Type of output "global" or "local" (default: "global")
         */
        virtual void define_output(const std::string &path, const std::string &filename = "results", const std::string &outputtype = "global");
        
        /**
         * @brief Write output to file.
         * @param so Solver output configuration
         * @param block Current block number
         * @param step Current step number
         * @param inc Current increment number
         * @param sub_inc Current sub-increment number
         * @param time Current time
         * @param outputtype Type of output "global" or "local" (default: "global")
         */
        virtual void output(const solver_output &so, const int &block, const int &step, const int &inc, const int &sub_inc, const double &time, const std::string &outputtype = "global");
    
    
        /**
         * @brief Stream output operator.
         * @param os Output stream
         * @param pc Phase characteristics to output
         * @return Output stream
         */
        friend std::ostream& operator << (std::ostream& os, const phase_characteristics &pc);
};


/** @} */ // end of phase group

} //namespace simcoon
