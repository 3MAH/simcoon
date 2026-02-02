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

///@file layer.hpp
///@brief Characteristics of a layer geometry, which inherits from geometry
///@version 2.0
///
///@details This class represents planar layers for laminate composite
///      homogenization. Layers are defined by their orientation using
///      Euler angles and links to adjacent layers.
///
///@example Python interface with JSON I/O:
///      @code{.py}
///      from simcoon.solver.micromechanics import Layer, GeometryOrientation
///      layer = Layer(number=0, concentration=0.5,
///                    geometry_orientation=GeometryOrientation(psi=0, theta=90, phi=-90))
///      @endcode

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace simcoon{

/**
 * @file layer.hpp
 * @brief Layer geometry for laminate homogenization.
 */

/** @addtogroup geometry
 *  @{
 */

/**
 * @brief Class representing a layer geometry for laminate homogenization.
 * 
 * This class extends the geometry base class to describe planar layers
 * for laminate composite homogenization. Each layer has a defined
 * orientation using Euler angles and links to adjacent layers.
 */
class layer : public geometry
{
	private:

	protected:

	public :
    
        int layerup;        ///< Index of the layer above this one (-1 if none)
        int layerdown;      ///< Index of the layer below this one (-1 if none)
    
        double psi_geom;    ///< First Euler angle for layer orientation (radians)
        double theta_geom;  ///< Second Euler angle for layer orientation (radians)
        double phi_geom;    ///< Third Euler angle for layer orientation (radians)
    
        /**
         * @brief Default constructor.
         */
		layer();

        /**
         * @brief Constructor with full parameters.
         * @param concentration Volume fraction of the layer
         * @param layerup Index of the layer above
         * @param layerdown Index of the layer below
         * @param psi First Euler angle (radians)
         * @param theta Second Euler angle (radians)
         * @param phi Third Euler angle (radians)
         */
        layer(const double &concentration, const int &layerup, const int &layerdown, 
              const double &psi, const double &theta, const double &phi);

        /**
         * @brief Copy constructor.
         * @param lay The layer to copy
         */
		layer(const layer &lay);

        /**
         * @brief Virtual destructor.
         */
        virtual ~layer();
    
        /**
         * @brief Assignment operator.
         * @param lay The layer to assign
         * @return Reference to this layer
         */
		virtual layer& operator = (const layer &lay);
		
        /**
         * @brief Output stream operator.
         */
        friend std::ostream& operator << (std::ostream &os, const layer &lay);
};


/** @} */ // end of geometry group

} //namespace simcoon
