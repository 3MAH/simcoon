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

///@file phase_multi.hpp
///@brief Characteristics of a phase for multiphase simulations, the parent class of:
// - ellipsoid_multi
// - layer_multi
// - cylinder_multi
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>

namespace simcoon{

//======================================
class phase_multi
//======================================
{
	private:

	protected:

	public :
		
		arma::mat A;	//Concentration tensor (strain)
		arma::mat A_start;	//Concentration tensor (strain)
		arma::mat B;	//Concentration tensor (stress)
		arma::mat B_start;	//Concentration tensor (stress)

        arma::vec A_in;	//Inelastic concentration tensor (strain vector)
        
		phase_multi(); 	//default constructor
        phase_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&); //Constructor with parameters
		phase_multi(const phase_multi&);	//Copy constructor
        ~phase_multi();
		
        void to_start();
        void set_start();
    
		virtual phase_multi& operator = (const phase_multi&);
		
        friend std::ostream& operator << (std::ostream&, const phase_multi&);
};

} //namespace simcoon
