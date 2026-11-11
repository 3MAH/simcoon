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

///@file state_variables_M.hpp
///@brief State variables of a phase, in a defined coordinate system:
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include <simcoon/Simulation/Phase/state_variables.hpp>


namespace simcoon{

/**
 * @file state_variables_M.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */


//======================================
class state_variables_M : public state_variables
//======================================
{
	private:

	protected:

	public :
		
    
        arma::vec::fixed<4> Wm; ///< mechanical works [total, reversible, irreversible, dissipated]; crosses both ways
        arma::vec::fixed<4> Wm_start; ///< set by set_start on each copy, never crosses
    
		arma::mat::fixed<6,6> L; ///< elastic stiffness; plugin-kernel input, crosses g2l only (the solver zeroes the global one)
		arma::mat::fixed<6,6> Lt; ///< tangent stiffness; crosses both ways
		
		/** @brief Default constructor (zero strains and stresses, identity F, U, R, DR). */
		state_variables_M();
        /**
         * @brief Constructor with parameters; asserts the Voigt and matrix sizes.
         * @param Etot Green-Lagrange strain
         * @param DEtot increment of Green-Lagrange strain
         * @param etot logarithmic strain
         * @param Detot increment of logarithmic strain
         * @param PKII 2nd Piola-Kirchhoff stress
         * @param PKII_start 2nd Piola-Kirchhoff stress at start
         * @param tau Kirchhoff stress
         * @param tau_start Kirchhoff stress at start
         * @param sigma Cauchy stress
         * @param sigma_start Cauchy stress at start
         * @param F0 deformation gradient at start
         * @param F1 deformation gradient at end
         * @param U0 stretch tensor at start
         * @param U1 stretch tensor at end
         * @param R rotation tensor
         * @param DR rotation increment
         * @param T temperature
         * @param DT temperature increment
         * @param nstatev number of internal state variables
         * @param statev internal state variables
         * @param statev_start internal state variables at start
         * @param Wm mechanical works (size 4)
         * @param Wm_start mechanical works at start (size 4)
         * @param L elastic stiffness (6x6)
         * @param Lt tangent stiffness (6x6)
         */
        state_variables_M(const arma::vec &Etot, const arma::vec &DEtot, const arma::vec &etot, const arma::vec &Detot, const arma::vec &PKII, const arma::vec &PKII_start, const arma::vec &tau, const arma::vec &tau_start, const arma::vec &sigma, const arma::vec &sigma_start, const arma::mat &F0, const arma::mat &F1, const arma::mat &U0, const arma::mat &U1, const arma::mat &R, const arma::mat &DR, const double &T, const double &DT, const int &nstatev, const arma::vec &statev, const arma::vec &statev_start, const arma::vec &Wm, const arma::vec &Wm_start, const arma::mat &L, const arma::mat &Lt);
    
		/** @brief Copy constructor. */
		state_variables_M(const state_variables_M &);
		/** @brief Destructor. */
		virtual ~state_variables_M();
		
		/** @brief Copy assignment: state_variables::operator= (tangent_mode, ncp and statev included), then the mechanical members. */
		state_variables_M& operator = (const state_variables_M&);
		
        /**
         * @brief Copy the field values from another state_variables_M object.
         * The common fields of state_variables::copy_fields (U0, U1, R and DR included), then
         * L, Lt, Wm and Wm_start; tangent_mode, ncp and statev untouched.
         */
		virtual state_variables_M& copy_fields_M (const state_variables_M&);
		
        using state_variables::update;
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &); //Initialize with parameters
        virtual void to_start(); //rollback: Wm_start goes to Wm
        virtual void set_start(const int &); //accept: Wm goes to Wm_start
    
        using state_variables::rotate_l2g;
        /// state_variables::rotate_l2g plus Wm and the rotated tangent Lt (ownership on the declarations).
        virtual state_variables_M& rotate_l2g(const state_variables_M&, const frame_rotation&);
        using state_variables::rotate_g2l;
        /// state_variables::rotate_g2l plus Wm, L and Lt.
        virtual state_variables_M& rotate_g2l(const state_variables_M&, const frame_rotation&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables_M&);
};


/** @} */ // end of phase group

} //namespace simcoon
