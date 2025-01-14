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

//======================================
class state_variables
//======================================
{
	private:

	protected:

	public :

		arma::vec Etot; // Green-Lagrange strain
		arma::vec DEtot; // increment of Green-Lagrange strain
		arma::vec etot; // logarithmic strain
		arma::vec Detot; // increment of logarithmic strain
		arma::vec PKII; //2nd Piola Kirchoff stress
		arma::vec PKII_start; //2nd Piola Kirchoff stress
		arma::vec tau; //Kirchoff stress
		arma::vec tau_start; //Kirchoff stress
		arma::vec sigma; //Cauchy stress
		arma::vec sigma_start; //Cauchy stress
        arma::mat F0;
        arma::mat F1;
        arma::mat U0;
        arma::mat U1;
        arma::mat R;
        arma::mat DR;
        double T;
        double DT;
    
        int nstatev;
        arma::vec statev;
        arma::vec statev_start;
    
        natural_basis nb;
    
		state_variables(); 	//default constructor
		state_variables(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for statev
        state_variables(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const natural_basis &); //Constructor with parameters
		state_variables(const state_variables &);	//Copy constructor
		virtual ~state_variables();
		
		virtual state_variables& operator = (const state_variables&);
		
		virtual state_variables& copy_fields (const state_variables&);

		virtual void resize();	//constructor - allocates memory for statev
		virtual void resize(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for statev
        virtual void update(const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const arma::mat &, const double &, const double &, const int &, const arma::vec &, const arma::vec &, const natural_basis &); //Initialize with parameters
		virtual int dimstatev () const {return nstatev;}       // returns the number of statev, nstatev    
        virtual void to_start(); //sigma goes to sigma_start
        virtual void set_start(const int &); //sigma_start goes to sigma
    
		virtual arma::mat PKI_stress();
		virtual arma::mat PKI_stress_start();
		virtual arma::mat Biot_stress();				
		virtual arma::mat Biot_stress_start();

        virtual state_variables& rotate_l2g(const state_variables&, const double&, const double&, const double&);
        virtual state_variables& rotate_g2l(const state_variables&, const double&, const double&, const double&);
    
        friend std::ostream& operator << (std::ostream&, const state_variables&);
};

} //namespace simcoon
