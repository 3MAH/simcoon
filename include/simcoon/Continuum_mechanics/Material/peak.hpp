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

///@file peak.hpp
///@brief Characteristics of a peak
///@version 0.9

#pragma once

#include <iostream>
#include <string>
#include <armadillo>

namespace simcoon{

/**
 * @file peak.hpp
 * @brief Peak characteristics for ODF/PDF.
 */

/** @addtogroup material
 *  @{
 */


//======================================
class peak
//======================================
{
	private:

	protected:

	public :

        int number;
        int method;
		double mean;
		double s_dev;
		double width;
        double ampl;
        arma::vec params; //Potential additional parameters
        arma::vec densities; //density
    
		peak(); 	//default constructor
    
        peak(const int &, const int &, const double &, const double &, const double &, const double &, const arma::vec &, const arma::vec &);

		peak(const peak&);	//Copy constructor
        virtual ~peak();

        virtual double get_density_ODF(const double &);

        virtual double get_density_PDF(const double &);
    
		virtual peak& operator = (const peak&);
    
        friend std::ostream& operator << (std::ostream&, const peak&);
};


/** @} */ // end of material group

} //namespace simcoon
