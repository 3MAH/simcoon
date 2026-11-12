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

///@file state_variables_T.hpp
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
    
        arma::vec::fixed<4> Wm; ///< mechanical works; crosses both ways
        arma::vec::fixed<3> Wt; ///< thermal works; crosses both ways
        arma::vec::fixed<4> Wm_start; ///< set by set_start on each copy, never crosses
        arma::vec::fixed<3> Wt_start; ///< set by set_start on each copy, never crosses
		
        arma::mat::fixed<6,6> dSdE; ///< mechanical tangent; kernel output, crosses l2g only
        arma::mat::fixed<6,6> dSdEt; ///< read by nobody, never crosses
        arma::vec::fixed<6> dSdT; ///< thermal stress tangent (Voigt); kernel output, crosses l2g only
        double Q; ///< heat flux (rate), set by the solver as -r/DTime (0 at DTime = 0); crosses g2l only
        /**
         * Heat of the increment per unit reference volume (MPa), storage included:
         * \f$ r = \Delta\gamma_{loc} - \theta_{n+1}\,\Delta\eta \f$ (heat_source_increment, tangent_assembly.hpp),
         * i.e. the time increment times the heat source. \f$ r = 0 \f$ is the adiabatic condition, also at
         * \f$ \Delta t = 0 \f$. Kernel output, crosses l2g only.
         */
        double r;
        double r_in; ///< never crosses
    
        /**
         * Strain derivative of the heat of the increment, \f$ \partial r / \partial \boldsymbol{\varepsilon} \f$ (Voigt);
         * kernel output, crosses l2g only. Dual to the engineering strain, it rotates with the STRESS
         * operator \f$ \mathbf{Q}_S = \mathbf{Q}_E^{-T} \f$: \f$ \partial r/\partial \boldsymbol{\varepsilon}'
         * = \mathbf{Q}_E^{-T} \, \partial r/\partial \boldsymbol{\varepsilon} \f$ for
         * \f$ \boldsymbol{\varepsilon}' = \mathbf{Q}_E \boldsymbol{\varepsilon} \f$.
         */
        arma::vec::fixed<6> drdE;
        double drdT; ///< temperature derivative of the heat of the increment r; kernel output, crosses l2g only

		/** @brief Default constructor (zero strains and stresses, identity F, U, R, DR). */
		state_variables_T();
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
         * @param Q heat flux (rate)
         * @param r heat of the increment
         * @param r_in stored as is (read by nobody in-tree)
         * @param Wm mechanical works (size 4)
         * @param Wt thermal works (size 3)
         * @param Wm_start mechanical works at start (size 4)
         * @param Wt_start thermal works at start (size 3)
         * @param dSdE mechanical tangent (6x6)
         * @param dSdEt stored as is (6x6, read by nobody)
         * @param dSdT thermal stress tangent (6)
         * @param drdE strain derivative of r (6)
         * @param drdT temperature derivative of r
         */
        state_variables_T(const arma::vec &Etot, const arma::vec &DEtot, const arma::vec &etot, const arma::vec &Detot, const arma::vec &PKII, const arma::vec &PKII_start, const arma::vec &tau, const arma::vec &tau_start, const arma::vec &sigma, const arma::vec &sigma_start, const arma::mat &F0, const arma::mat &F1, const arma::mat &U0, const arma::mat &U1, const arma::mat &R, const arma::mat &DR, const double &T, const double &DT, const int &nstatev, const arma::vec &statev, const arma::vec &statev_start, const double &Q, const double &r, const double &r_in, const arma::vec &Wm, const arma::vec &Wt, const arma::vec &Wm_start, const arma::vec &Wt_start, const arma::mat &dSdE, const arma::mat &dSdEt, const arma::vec &dSdT, const arma::vec &drdE, const double &drdT);
		/** @brief Copy constructor. */
		state_variables_T(const state_variables_T &);
		/** @brief Destructor. */
		virtual ~state_variables_T();
		
		/** @brief Copy assignment: state_variables::operator= (tangent_mode, ncp and statev included), then the thermomechanical members. */
		state_variables_T& operator = (const state_variables_T&);
		
        /**
         * @brief Copy the field values from another state_variables_T object.
         * The common fields of state_variables::copy_fields (U0, U1, R and DR included), then
         * the operators dSdE, dSdEt, dSdT, drdE, drdT, Wm, Wt, Wm_start,
         * Wt_start, Q, r and r_in; tangent_mode, ncp and statev untouched.
         */
		virtual state_variables_T& copy_fields_T (const state_variables_T&);
		
		using state_variables::update;
		virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const double &, const double &, const double &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::vec &, const arma::vec &, const double &);
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
