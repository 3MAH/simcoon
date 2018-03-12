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

///@file parameters.hpp
///@brief Handle of input parameters
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class constants
//======================================
{
	private:

	protected:

	public :
    
		int number;     //number of the constant
        double value;       //Value of the constant
    
        arma::vec input_values;   //values of the constant for each input file considered (test)
    
        std::string key;    //A unique key utilized to replace the constants in file(s)
        int ninput_files;
        std::vector<std::string> input_files;   //vector of files impacted (automaticaly filed for some parameter types)
    
		constants(); 	//default constructor
		constants(const int&, const int&);	//constructor - number, min and max values
        constants(const int&, const double&, const arma::vec&, const std::string&, const int &, const std::vector<std::string>&); //Constructor with parameters
		constants(const constants &);	//Copy constructor
		~constants();
		
		int dimfiles () const {return ninput_files;}  // returns the number of files associated to this parameter
    
        void update(const int&); //Update value based on the number in the vec input_value
        void resize(const int&, const int&);
				
		virtual constants& operator = (const constants&);
		
        friend  std::ostream& operator << (std::ostream&, const constants&);

};

} //namespace simcoon
