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

///@file step_meca.hpp
///@brief object that defines a mechanical step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "step.hpp"
#include "../Phase/state_variables_M.hpp"

namespace simcoon{

/**
 * @file step_meca.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


//======================================
class step_meca : public step
//======================================
{
private:
    
protected:
    
	public :

    arma::Col<int> cBC_meca; //1 is for stress (flux), 0 if for strain (state)
    arma::vec BC_meca;
    arma::mat mecas;
    arma::mat BC_mecas;
    arma::mat BC_w;
    arma::mat BC_R;
    double BC_T;
    int cBC_T;
    arma::vec Ts;
    arma::vec BC_Ts;
    
    step_meca(); 	//default constructor
    step_meca(const unsigned int &size); 	//constructor that allocates BC_meca and cBC_meca
    step_meca(const int &number, const double &Dn_init, const double &Dn_mini, const double &Dn_inc, const int &mode, const unsigned int &control_type, const arma::Col<int>& cBC_meca, const arma::vec& BC_meca, const arma::mat& mecas, const arma::mat& BC_mecas, const double& BC_T, const int& cBC_T, const arma::vec& Ts, const arma::vec& BC_Ts, const arma::mat& BC_w, const arma::mat& BC_R); //Constructor with parameters

    step_meca(const step_meca &s); 	//Copy constructor
    virtual ~step_meca();
    
    using step::generate;
    virtual void generate(const double&, const arma::vec &, const arma::vec&, const double&);
    virtual void generate_kin(const double&, const arma::mat &, const double &);
    virtual void assess_inc(const double &, double &, const double &, phase_characteristics &, double &, const double &, const arma::mat &, const int &);
    
    virtual step_meca& operator = (const step_meca &s);

    friend  std::ostream& operator << (std::ostream& os, const step_meca &s);
};


/** @} */ // end of solver group

} //namespace simcoon
