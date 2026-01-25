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

///@file cylinder_multi.hpp
///@brief characteristics of an cylinder for homogenization purposes
///@version 1.0

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>

namespace simcoon{

/**
 * @file cylinder_multi.hpp
 * @brief Cylindrical inclusion homogenization functions.
 */

/** @addtogroup homogenization
 *  @{
 */


//======================================
class cylinder_multi : public phase_multi
//======================================
{
private:
    
protected:
    
	public :

    arma::mat T_loc;
    arma::mat T;
    
    arma::mat A_loc;
    arma::mat B_loc;
    
    cylinder_multi(); //default constructor
    cylinder_multi(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::vec&, const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&); //Constructor with parameters
    cylinder_multi(const cylinder_multi&);	//Copy constructor
    ~cylinder_multi();
    
//    virtual void l2g_T();
//    virtual void g2l_T();
    
	virtual cylinder_multi& operator = (const cylinder_multi&);
    
    friend std::ostream& operator << (std::ostream&, const cylinder_multi&);
    
};


/** @} */ // end of homogenization group

} //namespace simcoon
