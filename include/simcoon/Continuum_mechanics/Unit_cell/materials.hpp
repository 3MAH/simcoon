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
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>

namespace simcoon{
    
//======================================
class aba_material : public material_characteristics
//======================================
{
private:
    
protected:
    
    public :
    
    int id;     //id key of the considered phase
    int nstatev; //Taken from the state_variables object
    double density;
    double conductivity;
    
    aba_material(); 	//default constructor
    aba_material(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for props
    
    aba_material(const int &, const int &, const std::string &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const int &, const arma::vec &); //constructor with parameters (2nd : id, 9th : nstatev)
    
    aba_material(const aba_material&);	//Copy constructor
    virtual ~aba_material();
    
    using material_characteristics::update;
    virtual void update(const int &, const int &, const std::string &, const int &, const double &, const double &, const double &, const double &, const double &, const int &, const int &, const arma::vec &);
    virtual void update(const material_characteristics&, const double &, const double &, const state_variables&, const int &);
    
    virtual int dimstatev () const {return nstatev;}       // returns the number of props, nprops
    
    virtual aba_material& operator = (const aba_material&);
    
    virtual void write(const unsigned int &, const std::string &, const std::string &inputfile, const bool & = true);
    
    friend std::ostream& operator << (std::ostream&, const aba_material&);
};
    
} //namespace simcoon
