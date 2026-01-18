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

///@file material_characteristics.hpp
///@brief Characteristics of a material, the parent class of:
// - phase_characteristics
// - ellipsoid_characteristics
// - layer_characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file material_characteristics.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */


/**
 * @brief Base class for material characteristics.
 * 
 * This class stores material properties and identification information used by constitutive models.
 * It serves as the parent class for phase_characteristics and other specialized material classes.
 * 
 * @details The material is characterized by:
 * - A UMAT name identifying the constitutive model
 * - Material orientation angles (Euler angles \f$ \psi, \theta, \phi \f$)
 * - Material properties vector (elastic constants, hardening parameters, etc.)
 */
class material_characteristics
{
	private:

	protected:

	public :

		int number; ///< Material identification number
        std::string umat_name; ///< Name of the constitutive model (UMAT)
        int save; ///< Flag indicating if results should be saved (1) or not (0)
        double psi_mat; ///< First Euler angle for material orientation (rad)
        double theta_mat; ///< Second Euler angle for material orientation (rad)
        double phi_mat; ///< Third Euler angle for material orientation (rad)
        
		int nprops; ///< Number of material properties
		arma::vec props; ///< Vector of material properties
    
        /**
         * @brief Default constructor.
         */
        material_characteristics();
        
        /**
         * @brief Constructor allocating memory for material properties.
         * @param nprops Number of material properties
         * @param init If true, initialize arrays to zero (default: true)
         * @param value Initial value for properties (default: 0.0)
         */
        material_characteristics(const int &nprops, const bool &init = true, const double &value = 0.);
    
        /**
         * @brief Full constructor with all parameters.
         * @param number Material identification number
         * @param umat_name Name of the constitutive model
         * @param save Flag for saving results
         * @param psi_mat First Euler angle (rad)
         * @param theta_mat Second Euler angle (rad)
         * @param phi_mat Third Euler angle (rad)
         * @param nprops Number of properties
         * @param props Properties vector
         */
        material_characteristics(const int &number, const std::string &umat_name, const int &save, const double &psi_mat, const double &theta_mat, const double &phi_mat, const int &nprops, const arma::vec &props);

        /**
         * @brief Copy constructor.
         * @param mc Material characteristics to copy
         */
        material_characteristics(const material_characteristics &mc);
        
        /**
         * @brief Virtual destructor.
         */
        virtual ~material_characteristics();

        /**
         * @brief Resize arrays to default size.
         */
		virtual void resize();
        
        /**
         * @brief Resize arrays for given number of properties.
         * @param nprops Number of material properties
         * @param init If true, initialize arrays to zero
         * @param value Initial value for properties
         */
        virtual void resize(const int &nprops, const bool &init = true, const double &value = 0.);
        
        /**
         * @brief Update all material characteristics.
         * @param number Material identification number
         * @param umat_name Name of the constitutive model
         * @param save Flag for saving results
         * @param psi_mat First Euler angle (rad)
         * @param theta_mat Second Euler angle (rad)
         * @param phi_mat Third Euler angle (rad)
         * @param nprops Number of properties
         * @param props Properties vector
         */
        virtual void update(const int &number, const std::string &umat_name, const int &save, const double &psi_mat, const double &theta_mat, const double &phi_mat, const int &nprops, const arma::vec &props);
        
        /**
         * @brief Get the number of material properties.
         * @return Number of properties (nprops)
         */
		virtual int dimprops () const {return nprops;}
    
        /**
         * @brief Assignment operator.
         * @param mc Material characteristics to assign
         * @return Reference to this object
         */
        virtual material_characteristics& operator = (const material_characteristics &mc);
		
        /**
         * @brief Stream output operator.
         * @param os Output stream
         * @param mc Material characteristics to output
         * @return Output stream
         */
        friend std::ostream& operator << (std::ostream& os, const material_characteristics &mc);
};


/** @} */ // end of phase group

} //namespace simcoon
