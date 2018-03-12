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
#include <simcoon/simulation/Phase/state_variables.hpp>

namespace simcoon{

//======================================
class state_variables_M : public state_variables
//======================================
{
	private:

	protected:

	public :
		
        arma::vec sigma_in;
        arma::vec sigma_in_start;
    
        arma::vec Wm;
        arma::vec Wm_start;
    
		arma::mat L;
		arma::mat Lt;
		
		state_variables_M(); 	//default constructor
        state_variables_M(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::vec &, const arma::vec &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &); //Constructor with parameters
    
		state_variables_M(const state_variables_M &);	//Copy constructor
		virtual ~state_variables_M();
		
		virtual state_variables_M& operator = (const state_variables_M&);
		
		virtual state_variables_M& copy_fields_M (const state_variables_M&);
		
        using state_variables::update;
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::vec &, const arma::vec &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &); //Initialize with parameters
        virtual void to_start(); //Wm goes to Wm_start
        virtual void set_start(); //Wm_start goes to Wm
    
        using state_variables::rotate_l2g;
        virtual state_variables_M& rotate_l2g(const state_variables_M&, const double&, const double&, const double&);
        using state_variables::rotate_g2l;
        virtual state_variables_M& rotate_g2l(const state_variables_M&, const double&, const double&, const double&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables_M&);
};

} //namespace simcoon
