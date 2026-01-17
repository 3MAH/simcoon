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

///@file state_variables.hpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

namespace simcoon{

/**
 * @file state_variables.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */


/**
 * @brief Base class representing state variables of a material phase.
 * 
 * This class stores all mechanical state variables (strains, stresses, deformation gradients)
 * in both current and reference configurations. It supports finite strain formulations with
 * various stress measures (Cauchy, Kirchhoff, 2nd Piola-Kirchhoff) and strain measures
 * (Green-Lagrange, logarithmic).
 * 
 * @details The class maintains:
 * - Strain measures: Green-Lagrange \f$ \mathbf{E} \f$, logarithmic \f$ \boldsymbol{\varepsilon} \f$
 * - Stress measures: Cauchy \f$ \boldsymbol{\sigma} \f$, Kirchhoff \f$ \boldsymbol{\tau} \f$, 
 *   2nd Piola-Kirchhoff \f$ \mathbf{S} \f$
 * - Deformation: Deformation gradient \f$ \mathbf{F} \f$, stretch \f$ \mathbf{U} \f$, rotation \f$ \mathbf{R} \f$
 * - Internal state variables for constitutive laws
 */
class state_variables
{
	private:

	protected:

	public :

		arma::vec Etot; ///< Green-Lagrange strain tensor (Voigt notation)
		arma::vec DEtot; ///< Increment of Green-Lagrange strain
		arma::vec etot; ///< Logarithmic (Hencky) strain tensor (Voigt notation)
		arma::vec Detot; ///< Increment of logarithmic strain
		arma::vec PKII; ///< 2nd Piola-Kirchhoff stress tensor (Voigt notation)
		arma::vec PKII_start; ///< 2nd Piola-Kirchhoff stress at start of increment
		arma::vec tau; ///< Kirchhoff stress tensor (Voigt notation)
		arma::vec tau_start; ///< Kirchhoff stress at start of increment
		arma::vec sigma; ///< Cauchy stress tensor (Voigt notation)
		arma::vec sigma_start; ///< Cauchy stress at start of increment
        arma::mat F0; ///< Deformation gradient at start of increment
        arma::mat F1; ///< Deformation gradient at end of increment
        arma::mat U0; ///< Right stretch tensor at start of increment
        arma::mat U1; ///< Right stretch tensor at end of increment
        arma::mat R; ///< Rotation tensor
        arma::mat DR; ///< Increment of rotation tensor
        double T; ///< Current temperature
        double DT; ///< Temperature increment
    
        int nstatev; ///< Number of internal state variables
        arma::vec statev; ///< Internal state variables vector
        arma::vec statev_start; ///< Internal state variables at start of increment
    
        natural_basis nb; ///< Natural basis for covariant/contravariant operations
    
        /**
         * @brief Default constructor.
         */
		state_variables();
        
        /**
         * @brief Constructor allocating memory for state variables.
         * @param nstatev Number of internal state variables
         * @param init If true, initialize arrays to zero (default: true)
         * @param value Initial value for state variables (default: 0.0)
         */
		state_variables(const int &, const bool& = true, const double& = 0.);
        
        /**
         * @brief Full constructor with all parameters.
         * @param Etot Green-Lagrange strain
         * @param DEtot Increment of Green-Lagrange strain
         * @param etot Logarithmic strain
         * @param Detot Increment of logarithmic strain
         * @param PKII 2nd Piola-Kirchhoff stress
         * @param PKII_start 2nd Piola-Kirchhoff stress at start
         * @param tau Kirchhoff stress
         * @param tau_start Kirchhoff stress at start
         * @param sigma Cauchy stress
         * @param sigma_start Cauchy stress at start
         * @param F0 Deformation gradient at start
         * @param F1 Deformation gradient at end
         * @param U0 Stretch tensor at start
         * @param U1 Stretch tensor at end
         * @param R Rotation tensor
         * @param DR Rotation increment
         * @param T Temperature
         * @param DT Temperature increment
         * @param nstatev Number of state variables
         * @param statev State variables
         * @param statev_start State variables at start
         * @param nb Natural basis
         */
        state_variables(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const natural_basis &);
        
        /**
         * @brief Copy constructor.
         * @param sv State variables to copy
         */
		state_variables(const state_variables &);
        
        /**
         * @brief Virtual destructor.
         */
		virtual ~state_variables();
		
        /**
         * @brief Assignment operator.
         * @param sv State variables to assign
         * @return Reference to this object
         */
		virtual state_variables& operator = (const state_variables&);
		
        /**
         * @brief Copy field values from another state_variables object.
         * @param sv Source state variables
         * @return Reference to this object
         */
		virtual state_variables& copy_fields (const state_variables&);

        /**
         * @brief Resize arrays to default size.
         */
		virtual void resize();
        
        /**
         * @brief Resize arrays for given number of state variables.
         * @param nstatev Number of internal state variables
         * @param init If true, initialize arrays to zero
         * @param value Initial value for state variables
         */
		virtual void resize(const int &, const bool& = true, const double& = 0.);
        
        /**
         * @brief Update all state variables with new values.
         */
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const natural_basis &);
        
        /**
         * @brief Get the number of internal state variables.
         * @return Number of state variables (nstatev)
         */
		virtual int dimstatev () const {return nstatev;}
        
        /**
         * @brief Copy current values to start-of-increment values.
         */
        virtual void to_start();
        
        /**
         * @brief Set current values from start-of-increment values.
         * @param control Control flag for selective update
         */
        virtual void set_start(const int &);
    
        /**
         * @brief Compute 1st Piola-Kirchhoff stress tensor.
         * @return 1st Piola-Kirchhoff stress as 3x3 matrix
         */
		virtual arma::mat PKI_stress();
        
        /**
         * @brief Compute 1st Piola-Kirchhoff stress at start of increment.
         * @return 1st Piola-Kirchhoff stress as 3x3 matrix
         */
		virtual arma::mat PKI_stress_start();
        
        /**
         * @brief Compute Biot stress tensor.
         * @return Biot stress as 3x3 matrix
         */
		virtual arma::mat Biot_stress();
        
        /**
         * @brief Compute Biot stress at start of increment.
         * @return Biot stress as 3x3 matrix
         */
		virtual arma::mat Biot_stress_start();

        /**
         * @brief Rotate state variables from local to global frame.
         * @param sv Source state variables
         * @param psi First Euler angle (rad)
         * @param theta Second Euler angle (rad)
         * @param phi Third Euler angle (rad)
         * @return Reference to this object with rotated values
         */
        virtual state_variables& rotate_l2g(const state_variables&, const double&, const double&, const double&);
        
        /**
         * @brief Rotate state variables from global to local frame.
         * @param sv Source state variables
         * @param psi First Euler angle (rad)
         * @param theta Second Euler angle (rad)
         * @param phi Third Euler angle (rad)
         * @return Reference to this object with rotated values
         */
        virtual state_variables& rotate_g2l(const state_variables&, const double&, const double&, const double&);
    
        /**
         * @brief Stream output operator.
         * @param os Output stream
         * @param sv State variables to output
         * @return Output stream
         */
        friend std::ostream& operator << (std::ostream&, const state_variables&);
};


/** @} */ // end of phase group

} //namespace simcoon
