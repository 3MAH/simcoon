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

///@file step_thermomeca.hpp
///@brief object that defines a thermomechanical step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "step.hpp"
#include "../Phase/state_variables_T.hpp"

namespace simcoon{

/**
 * @file step_thermomeca.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


//======================================
class step_thermomeca : public step
//======================================
{
private:
    
protected:
    
	public :

    arma::Col<int> cBC_meca; //True is for stress (flux), false if for strain (state)
    arma::vec BC_meca;
    arma::mat mecas;
    arma::mat BC_mecas;
    arma::mat BC_w;
    arma::mat BC_R;    
    double BC_T;
    int cBC_T;         //True (1) is for a heat flux entering in a material point, 0 is for fixed temperature
    arma::vec Ts;
    arma::vec BC_Ts;
    
    /** @brief Default constructor: BC_T = 0, cBC_T = 0. */
    step_thermomeca();
    /**
     * @brief Constructor that allocates cBC_meca and BC_meca (size 9 for control_type 5/6, else 6), BC_w = 0, BC_R = I.
     * @param control_type control type of the step
     */
    step_thermomeca(const unsigned int &control_type);
    /**
     * @brief Constructor with parameters.
     * @param number step identification number
     * @param Dn_init initial time increment fraction
     * @param Dn_mini minimum time increment fraction
     * @param Dn_inc maximum time increment fraction
     * @param mode loading mode
     * @param control_type control type
     * @param cBC_meca type of each mechanical boundary condition (1 stress, 0 strain)
     * @param BC_meca target values of the mechanical boundary conditions
     * @param mecas mechanical increments, one row per increment
     * @param BC_mecas cumulated mechanical boundary values, one row per increment
     * @param BC_T target temperature or heat flux (see cBC_T)
     * @param cBC_T type of the thermal boundary condition (1 heat flux, 0 fixed temperature)
     * @param Ts temperature/heat increments, one per increment
     * @param BC_Ts cumulated thermal boundary values, one per increment
     * @param BC_w spin of the step (3x3)
     * @param BC_R rotation of the step (3x3)
     */
    step_thermomeca(const int &number, const double &Dn_init, const double &Dn_mini, const double &Dn_inc, const int &mode, const unsigned int &control_type, const arma::Col<int>& cBC_meca, const arma::vec& BC_meca, const arma::mat& mecas, const arma::mat& BC_mecas, const double& BC_T, const int& cBC_T, const arma::vec& Ts, const arma::vec& BC_Ts, const arma::mat& BC_w, const arma::mat& BC_R);
    /** @brief Copy constructor. */
    step_thermomeca(const step_thermomeca&) = default;
    /** @brief Destructor. */
    virtual ~step_thermomeca();
    
    using step::generate;
    virtual void generate(const double&, const arma::vec&, const arma::vec&, const double&);
    /**
     * @brief Generate the increment tables for a fully kinematic thermomechanical step
     * (control types 5/6); same geodesic-F conventions as step_meca::generate_kin.
     * @param[in] mTime current (absolute) time at the start of the step
     * @param[in] mF current 3x3 deformation gradient \f$ \mathbf{F} \f$ (the step anchor)
     * @param[in] mT current temperature
     */
    virtual void generate_kin(const double&, const arma::mat&m, const double &);    
    virtual void assess_inc(const double &, double &, const double &, phase_characteristics &, double &, const double &, const arma::mat &, const int &);
    
    /** @brief Copy assignment. */
    virtual step_thermomeca& operator = (const step_thermomeca&) = default;
        
    friend  std::ostream& operator << (std::ostream&, const step_thermomeca&);
};


/** @} */ // end of solver group

} //namespace simcoon
