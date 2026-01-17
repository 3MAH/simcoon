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

///@file ODF.hpp
///@brief Characteristics of an ODF
///@version 0.9

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Material/peak.hpp>

namespace simcoon{

/**
 * @file PDF.hpp
 * @brief Probability Distribution Function (PDF) management.
 */

/** @addtogroup material
 *  @{
 */


//======================================
class PDF
//======================================
{
	private:

	protected:

	public :

        int Nphases; //Number of phases for the ODF discretization
        int Parameter;  //Index of the material parameter to change
        int n_densities; //number of angle to evaluate the densities (mostly for plotting purposes)
    
        std::vector<peak> peaks; //peaks
        double norm;
    
        arma::vec densities; //sum of densities
        arma::vec limits; //minimal and maximal angles of orientation
    
		PDF(); 	//default constructor
		PDF(const int &, const double &, const double &); //Partial constructor
        PDF(const int &, const int &, const int &, const std::vector<peak> &, const double &, const arma::vec &, const arma::vec &); //Full constructor

		PDF(const PDF&);	//Copy constructor
        virtual ~PDF();

        void construct(const int&);
        double density(const double &);
    
		virtual PDF& operator = (const PDF&);
    
        friend std::ostream& operator << (std::ostream&, const PDF&);
};


/** @} */ // end of material group

} //namespace simcoon
