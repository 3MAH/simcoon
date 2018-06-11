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

//======================================
class step_meca : public step
//======================================
{
private:
    
protected:
    
	public :

    arma::Col<int> cBC_meca; //True is for stress (flux), false if for strain (state)
    arma::vec BC_meca;
    arma::mat mecas;
    double BC_T;
    int cBC_T;
    arma::vec Ts;
    
    step_meca(); 	//default constructor
    step_meca(const unsigned int &); 	//constructor that allocates BC_meca and cBC_meca
    step_meca(const int &, const double &, const double &, const double &, const int &, const unsigned int &, const arma::Col<int>&, const arma::vec&, const arma::mat&, const double&, const int&, const arma::vec&); //Constructor with parameters
    
    step_meca(const step_meca&);	//Copy constructor
    virtual ~step_meca();
    
    using step::generate;
    virtual void generate(const double&, const arma::vec&, const arma::vec&, const double&);
    virtual void assess_inc(const double &, double &, const double &, phase_characteristics &, double &, const double &);    
    
    virtual step_meca& operator = (const step_meca&);
        
    friend  std::ostream& operator << (std::ostream&, const step_meca&);
};

} //namespace simcoon
