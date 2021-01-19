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

///@file output.hpp
///@brief object that defines the output
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class solver_output
//======================================
{
private:
    
protected:
    
public :

    //output values
//    int o_nb_meca;
//    arma::Col<int> o_meca;
    int o_nb_strain;
    arma::Col<int> o_strain;
    int o_nb_stress;
    arma::Col<int> o_stress;
    int o_nb_T;
    
    int o_strain_type; //0 for Green-Lagrange, 1 for Biot, 2 for logarithmic
    int o_stress_type;  //0 for PKII, 1 for Nominal, 2 for PKI, 3 for Kirchoff, 4 for Cauchy
    int o_rotation_type;  //0 for nothing, 1 for rotation, 2 for spin, 3 for both
    
	int o_nw_statev;
	arma::Col<int> o_wanted_statev;
	arma::Col<int> o_range_statev;
    
    arma::Col<int> o_type;
    arma::Col<int> o_nfreq;
    arma::vec o_tfreq;
    
    solver_output(); 	//default constructor
    solver_output(const int&);	//Constructor with parameters
    solver_output(const solver_output &);	//Copy constructor
    ~solver_output();
    
    virtual solver_output& operator = (const solver_output&);
    
    friend  std::ostream& operator << (std::ostream&, const solver_output&);
};

} //namespace simcoon
