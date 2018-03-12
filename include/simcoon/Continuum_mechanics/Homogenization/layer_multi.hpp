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

///@file layer_multi.hpp
///@brief Characteristics of a layer, similar to a phase in this version for homogenization purposes
///@version 1.0
//
#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

//======================================
class layer_multi : public phase_multi
//======================================
{
	private:

	protected:
    
	public :

/*        static arma::mat sumDnn;
        static arma::mat sumcDsig;
        static double m;*/
        //Those static variables are redefined in fonction to help multithread in case of n-order micromechanics
        
        //Part of the tangent modulus (usefull derivative)
        arma::mat Dnn;
        arma::mat Dnt;
        //Derivatives of the gradient / x1
        arma::mat dXn;
        arma::mat dXt;
    
        //local and global increments of the concentration strain / x1
        arma::vec sigma_hat;
        arma::vec dzdx1;
    
        layer_multi(); 	//default constructor
        layer_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::vec&); //Constructor with parameters
        layer_multi(const layer_multi&);	//Copy constructor
            ~layer_multi();
    
        virtual layer_multi& operator = (const layer_multi&);
        
        friend std::ostream& operator << (std::ostream&, const layer_multi&);
    
};

} //namespace simcoon
