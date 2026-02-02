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

///@file cylinder.hpp
///@brief Characteristics of a cylindrical phase, which inherits from geometry
///@version 2.0
///
///@details This class represents cylindrical inclusions for micromechanical
///      homogenization schemes. Cylinders are defined by length L and radius R.
///
///@example Python interface with JSON I/O:
///      @code{.py}
///      from simcoon.solver.micromechanics import Cylinder
///      cyl = Cylinder(number=0, concentration=0.3, L=50, R=1)
///      print(cyl.aspect_ratio)  # 50.0
///      @endcode

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace simcoon{

/**
 * @file cylinder.hpp
 * @brief Cylindrical inclusion geometry for micromechanics.
 */

/** @addtogroup geometry
 *  @{
 */

/**
 * @brief Class representing a cylindrical inclusion geometry.
 *
 * This class extends the geometry base class to describe cylindrical inclusions
 * for micromechanical homogenization schemes. The cylinder is defined by its
 * length L and radius R, with orientation specified by Euler angles.
 *
 * The aspect ratio \f$ L/R \f$ determines fiber characteristics:
 * - High aspect ratio: Long fibers
 * - Low aspect ratio: Short fibers or disc-like inclusions
 */
class cylinder : public geometry
{
	private:

	protected:

	public :

        int coatingof;      ///< Index of the phase this cylinder is coating (0 if none)
        int coatedby;       ///< Index of the phase coating this cylinder (0 if none)

        double L;           ///< Length of the cylinder
        double R;           ///< Radius of the cylinder

        double psi_geom;    ///< First Euler angle for orientation (radians)
        double theta_geom;  ///< Second Euler angle for orientation (radians)
        double phi_geom;    ///< Third Euler angle for orientation (radians)

        /**
         * @brief Default constructor.
         */
		cylinder();

        /**
         * @brief Constructor with full parameters.
         * @param Lval Length of the cylinder
         * @param coatingof Index of the phase this cylinder coats
         * @param coatedby Index of the phase coating this cylinder
         * @param Rval Radius of the cylinder
         * @param psi_geom First Euler angle (radians)
         * @param theta_geom Second Euler angle (radians)
         * @param phi_geom Third Euler angle (radians)
         * @param dummy Unused parameter (for compatibility)
         */
    	cylinder(const double &Lval, const int &coatingof, const int &coatedby, const double &Rval, const double &psi_geom, const double &theta_geom, const double &phi_geom, const double &dummy);

        /**
         * @brief Copy constructor.
         * @param cyl The cylinder to copy
         */
		cylinder(const cylinder&);

        /**
         * @brief Virtual destructor.
         */
        virtual ~cylinder();

        /**
         * @brief Assignment operator.
         * @param cyl The cylinder to assign
         * @return Reference to this cylinder
         */
		virtual cylinder& operator = (const cylinder&);

        /**
         * @brief Output stream operator.
         */
        friend std::ostream& operator << (std::ostream& os, const cylinder &cyl);
};


/** @} */ // end of geometry group

} //namespace simcoon
