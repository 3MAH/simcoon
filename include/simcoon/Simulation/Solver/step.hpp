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

///@file step.hpp
///@brief object that defines an step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "output.hpp"
#include "../Phase/state_variables.hpp"
#include "../Phase/phase_characteristics.hpp"

namespace simcoon{

/**
 * @file step.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


/**
 * @brief Class representing a loading step within a simulation block.
 * 
 * A step defines a portion of the loading path with specified boundary conditions
 * and time discretization parameters. It controls how the load is applied and
 * how time increments are managed.
 * 
 * @details The step class manages:
 * - Time discretization (initial, minimum increments)
 * - Boundary condition mode (stress, strain, or mixed control)
 * - Loading path definition through external files
 */
class step
{
private:
    
protected:
    
	public :
    int number; ///< Step identification number
    double Dn_init; ///< Initial fraction of the step (initial time increment ratio)
    double Dn_mini; ///< Minimal fraction of the step (minimum time increment ratio)
    double Dn_inc; ///< Maximum fraction of the step (maximum time increment ratio)
    int ninc; ///< Number of milestones/increments in the step
    int mode; ///< Loading mode identifier
    unsigned int control_type; ///< Control type (0: strain, 1: stress, 2: mixed)
    
    arma::vec times; ///< Vector of time values for the step
    double BC_Time; ///< Boundary condition application time
    
    std::string file; ///< Input/output file for loading path values
    
    /**
     * @brief Default constructor.
     */
    step();
    
    /**
     * @brief Constructor with parameters.
     * @param number Step identification number
     * @param Dn_init Initial time increment fraction
     * @param Dn_mini Minimum time increment fraction
     * @param Dn_inc Maximum time increment fraction
     * @param mode Loading mode
     * @param control_type Control type
     */
    step(const int &number, const double &Dn_init, const double &Dn_mini, const double &Dn_inc, const int &mode, const unsigned int &control_type);
    
    /**
     * @brief Copy constructor.
     * @param s Step to copy
     */
    step(const step &s);
    
    /**
     * @brief Virtual destructor.
     */
    virtual ~step();
   
    /**
     * @brief Generate the time discretization for the step.
     */
    virtual void generate();
    
    /**
     * @brief Compute the next increment parameters.
     * @param tnew_dt Suggested new time increment ratio (output)
     * @param inc Current increment number
     * @param tinc Time increment (output)
     * @param Dtinc Delta time increment (output)
     * @param Dn Increment fraction (output)
     * @param control Increment control flag
     */
    virtual void compute_inc(double &tnew_dt, const int &inc, double &tinc, double &Dtinc, double &Dn, const int &control);
    
    /**
     * @brief Assignment operator.
     * @param s Step to assign
     * @return Reference to this object
     */
    virtual step& operator = (const step& s);
    
    /**
     * @brief Stream output operator.
     * @param os Output stream
     * @param s Step to output
     * @return Output stream
     */
    friend  std::ostream& operator << (std::ostream& os, const step &s);
};


/** @} */ // end of solver group

} //namespace simcoon
