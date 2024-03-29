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
    
    step_thermomeca(); 	//default constructor
    step_thermomeca(const unsigned int &); 	//constructor that allocates BC_meca and cBC_meca
    step_thermomeca(const int &, const double &, const double &, const double &, const int &, const unsigned int &, const arma::Col<int>&, const arma::vec&, const arma::mat&, const arma::mat&, const double&, const int&, const arma::vec&, const arma::vec&, const arma::mat&, const arma::mat&); //Constructor with parameters
    step_thermomeca(const step_thermomeca&);	//Copy constructor
    virtual ~step_thermomeca();
    
    using step::generate;
    virtual void generate(const double&, const arma::vec&, const arma::vec&, const double&);
    virtual void generate_kin(const double&, const arma::mat&m, const double &);    
    virtual void assess_inc(const double &, double &, const double &, phase_characteristics &, double &, const double &, const arma::mat &, const int &);
    
    virtual step_thermomeca& operator = (const step_thermomeca&);
        
    friend  std::ostream& operator << (std::ostream&, const step_thermomeca&);
};

} //namespace simcoon
