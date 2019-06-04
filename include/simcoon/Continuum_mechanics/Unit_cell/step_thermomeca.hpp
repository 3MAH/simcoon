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

///@file material_characteristics.hpp
///@brief Characteristics of a material, the parent class of:
// - phase_characteristics
// - ellipsoid_characteristics
// - layer_characteristics
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>

namespace simcoon{
    
//======================================
class aba_step_thermomeca : public step_thermomeca
//======================================
{
private:
    
protected:
    
    public :
    
    std::string name;     //name of the step
    bool nlgeom;
    int type; //Automatic or fixed
    double max_temp;
        
    aba_step_thermomeca(); 	//default constructor
    aba_step_thermomeca(const std::string &, const bool &, const int &, const double &, const int &, const double &, const double &, const double &, const int &, const unsigned int &, const arma::Col<int>&, const arma::vec&, const arma::mat&, const double&, const int&, const arma::vec&,  const arma::mat&, const arma::mat&); //Constructor with parameters
    
    aba_step_thermomeca(const aba_step_thermomeca&);	//Copy constructor
    virtual ~aba_step_thermomeca();

    using step::generate;
    virtual void update(const step_thermomeca&, const std::string &, const bool &, const int &, const int &);
    
    virtual aba_step_thermomeca& operator = (const aba_step_thermomeca&);
    
    virtual void write(const std::string &, const std::string &inputfile, const bool & = true);
    
    friend  std::ostream& operator << (std::ostream&, const aba_step_thermomeca&);
};
    
} //namespace simcoon
