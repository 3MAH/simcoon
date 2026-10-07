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
#include <simcoon/Simulation/Phase/state_variables.hpp>

namespace simcoon{

/**
 * @file state_variables_T.hpp
 * @brief Phase and state variable management.
 */

/** @addtogroup phase
 *  @{
 */


//======================================
    class state_variables_T : public state_variables
//======================================
{
	private:

	protected:

	public :
    
        arma::vec sigma_in; ///< inelastic stress; no in-tree writer (a plugin output in umat_plugin_api.hpp, unused), never crosses the frame
        arma::vec sigma_in_start; ///< never crosses
        arma::vec Wm; ///< mechanical works; crosses both ways
        arma::vec Wt; ///< thermal works; crosses both ways
        arma::vec Wm_start; ///< set by set_start on each copy, never crosses
        arma::vec Wt_start; ///< set by set_start on each copy, never crosses
		
        arma::mat dSdE; ///< mechanical tangent; kernel output, crosses l2g only
        arma::mat dSdEt; ///< read by nobody, never crosses
        arma::mat dSdT; ///< thermal stress tangent, a Voigt vector stored as a mat (1x6 by the constructors, 6x1 by the solver); kernel output, crosses l2g only
        double Q; ///< heat flux, set by the solver as -r; crosses g2l only
        double r; ///< heat source; kernel output, crosses l2g only
        double r_in; ///< never crosses
    
        /**
         * Heat source strain tangent \f$ \partial r / \partial \boldsymbol{\varepsilon} \f$, a
         * Voigt vector stored as a mat (1x6 by the constructors, 6x1 by the solver); kernel
         * output, crosses l2g only. Dual to the engineering strain, it rotates with the STRESS
         * operator \f$ \mathbf{Q}_S = \mathbf{Q}_E^{-T} \f$: \f$ \partial r/\partial \boldsymbol{\varepsilon}'
         * = \mathbf{Q}_E^{-T} \, \partial r/\partial \boldsymbol{\varepsilon} \f$ for
         * \f$ \boldsymbol{\varepsilon}' = \mathbf{Q}_E \boldsymbol{\varepsilon} \f$.
         */
        arma::mat drdE;
        arma::mat drdT; ///< heat source temperature tangent; kernel output, crosses l2g only

		state_variables_T(); 	//default constructor
    state_variables_T(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::vec &, const arma::vec &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &); //Constructor with parameters
		state_variables_T(const state_variables_T &);	//Copy constructor
		virtual ~state_variables_T();
		
		virtual state_variables_T& operator = (const state_variables_T&);
		
		virtual state_variables_T& copy_fields_T (const state_variables_T&);
		
		using state_variables::update;
		virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::vec &, const arma::vec &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &);
        virtual void to_start(); //rollback: Wm_start & Wt_start go to Wm & Wt, respectively
        virtual void set_start(const int &); //accept: Wm & Wt go to Wm_start & Wt_start, respectively
    
        using state_variables::rotate_l2g;
        /// state_variables::rotate_l2g plus the thermomechanical members (ownership on their declarations).
        virtual state_variables_T& rotate_l2g(const state_variables_T&, const frame_rotation&);
        using state_variables::rotate_g2l;
        /// state_variables::rotate_g2l plus the thermomechanical members.
        virtual state_variables_T& rotate_g2l(const state_variables_T&, const frame_rotation&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables_T&);
};


/** @} */ // end of phase group

} //namespace simcoon
