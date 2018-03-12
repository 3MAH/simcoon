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

namespace simcoon{

//======================================
class state_variables
//======================================
{
	private:

	protected:

	public :
		
		arma::vec Etot;
		arma::vec DEtot;
		arma::vec sigma;
		arma::vec sigma_start;
        arma::mat F0;
        arma::mat F1;
        double T;
        double DT;
    
        int nstatev;
        arma::vec statev;
        arma::vec statev_start;
    
		state_variables(); 	//default constructor
		state_variables(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for statev
    state_variables(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &,const double &, const double &, const int &, const arma::vec &, const arma::vec &); //Constructor with parameters
		state_variables(const state_variables &);	//Copy constructor
		virtual ~state_variables();
		
		virtual state_variables& operator = (const state_variables&);
		
		virtual state_variables& copy_fields (const state_variables&);

		virtual void resize();	//constructor - allocates memory for statev
		virtual void resize(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for statev
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &); //Initialize with parameters
		virtual int dimstatev () const {return nstatev;}       // returns the number of statev, nstatev    
        virtual void to_start(); //sigma goes to sigma_start
        virtual void set_start(); //sigma_start goes to sigma
    
        virtual state_variables& rotate_l2g(const state_variables&, const double&, const double&, const double&);
        virtual state_variables& rotate_g2l(const state_variables&, const double&, const double&, const double&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables&);
};

} //namespace simcoon
