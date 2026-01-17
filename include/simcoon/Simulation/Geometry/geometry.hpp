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

///@file geometry.hpp
///@brief Class that define the geometry of a phase
///@version 1.0

#pragma once

#include <iostream>
#include <string>

namespace simcoon{

/**
 * @file geometry.hpp
 * @brief Inclusion geometry functions.
 */

/** @addtogroup geometry
 *  @{
 */

/**
 * @brief Base class for defining the geometry of a phase.
 * 
 * This class provides a foundation for describing the geometric properties
 * of material phases in homogenization schemes. Derived classes include
 * ellipsoid, cylinder, and layer geometries.
 */
class geometry
{
	private:

	protected:

	public :

        double concentration;   ///< Volume fraction of the phase
    
        /**
         * @brief Default constructor.
         * @details Initializes an empty geometry object.
         */
		geometry();

        /**
         * @brief Constructor with concentration parameter.
         * @param c The volume fraction of the phase (double)
         */
		geometry(const double &c);

        /**
         * @brief Copy constructor.
         * @param geo The geometry object to copy
         */
		geometry(const geometry &geo);

        /**
         * @brief Virtual destructor.
         */
        virtual ~geometry();
    
        /**
         * @brief Assignment operator.
         * @param geo The geometry object to assign
         * @return Reference to this geometry object
         */
		virtual geometry& operator = (const geometry &geo);
		
        /**
         * @brief Output stream operator.
         * @param os Output stream
         * @param geo Geometry object to output
         * @return Reference to the output stream
         */
        friend std::ostream& operator << (std::ostream &os, const geometry &geo);
};


/** @} */ // end of geometry group

} //namespace simcoon
