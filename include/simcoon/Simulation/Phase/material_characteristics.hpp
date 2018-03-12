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

namespace simcoon{

//======================================
class material_characteristics
//======================================
{
	private:

	protected:

	public :

		int number;
        std::string umat_name;
        int save;   //If the restults of this material being saved or not
        double psi_mat;
        double theta_mat;
        double phi_mat;
        
		int nprops;
		arma::vec props;
    
		material_characteristics(); 	//default constructor
		material_characteristics(const int &, const bool& = true, const double& = 0.);	//constructor - allocates memory for props
    
        material_characteristics(const int &, const std::string &, const int &, const double &, const double &, const double &, const int &, const arma::vec &);

		material_characteristics(const material_characteristics&);	//Copy constructor
        virtual ~material_characteristics();

		virtual void resize();
		virtual void resize(const int &, const bool & = true, const double & = 0.);
		virtual void update(const int &, const std::string &, const int &, const double &, const double &, const double &, const int &, const arma::vec &);
		virtual int dimprops () const {return nprops;}       // returns the number of props, nprops
    
		virtual material_characteristics& operator = (const material_characteristics&);
		
        friend std::ostream& operator << (std::ostream&, const material_characteristics&);
};

} //namespace simcoon
