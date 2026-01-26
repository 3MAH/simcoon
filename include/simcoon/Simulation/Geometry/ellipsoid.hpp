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

///@file ellipsoid.hpp
///@brief Characteristics of an ellipsoidal phase, which hereditates from:
///-phase characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include "geometry.hpp"

namespace simcoon{

/**
 * @file ellipsoid.hpp
 * @brief Inclusion geometry functions.
 */

/** @addtogroup geometry
 *  @{
 */

/**
 * @brief Class representing an ellipsoidal inclusion geometry.
 * 
 * This class extends the geometry base class to describe ellipsoidal inclusions
 * for micromechanical homogenization schemes. The ellipsoid is defined by three
 * semi-axes and three Euler angles for orientation.
 * 
 * The semi-axes ratios \f$ a_1, a_2, a_3 \f$ define the shape:
 * - Sphere: \f$ a_1 = a_2 = a_3 \f$
 * - Prolate spheroid: \f$ a_1 > a_2 = a_3 \f$
 * - Oblate spheroid: \f$ a_1 = a_2 > a_3 \f$
 * - General ellipsoid: \f$ a_1 \neq a_2 \neq a_3 \f$
 */
class ellipsoid : public geometry
{
	private:

	protected:

	public :

        int coatingof;      ///< Index of the phase this ellipsoid is coating (0 if none)
        int coatedby;       ///< Index of the phase coating this ellipsoid (0 if none)
    
        double a1;          ///< First semi-axis (relative)
        double a2;          ///< Second semi-axis (relative)
        double a3;          ///< Third semi-axis (relative)
        
        double psi_geom;    ///< First Euler angle for orientation (radians)
        double theta_geom;  ///< Second Euler angle for orientation (radians)
        double phi_geom;    ///< Third Euler angle for orientation (radians)
    
        /**
         * @brief Default constructor.
         */
		ellipsoid();

        /**
         * @brief Constructor with full parameters.
         * @param concentration Volume fraction of the phase
         * @param coatingof Index of the phase this ellipsoid coats
         * @param coatedby Index of the phase coating this ellipsoid
         * @param a1 First semi-axis
         * @param a2 Second semi-axis
         * @param a3 Third semi-axis
         * @param psi First Euler angle (radians)
         * @param theta Second Euler angle (radians)
         * @param phi Third Euler angle (radians)
         */
        ellipsoid(const double &concentration, const int &coatingof, const int &coatedby, 
                  const double &a1, const double &a2, const double &a3, 
                  const double &psi, const double &theta, const double &phi);

        /**
         * @brief Copy constructor.
         * @param ell The ellipsoid to copy
         */
		ellipsoid(const ellipsoid &ell);

        /**
         * @brief Virtual destructor.
         */
        virtual ~ellipsoid();
    
        /**
         * @brief Assignment operator.
         * @param ell The ellipsoid to assign
         * @return Reference to this ellipsoid
         */
		virtual ellipsoid& operator = (const ellipsoid &ell);
		
        /**
         * @brief Output stream operator.
         */
        friend std::ostream& operator << (std::ostream &os, const ellipsoid &ell);
};


/** @} */ // end of geometry group

} //namespace simcoon
