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

///@file opti_data.hpp
///@brief Handle of data from optimization
///@author Chemisky & Despringre
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>

namespace simcoon{

//======================================
class opti_data
//======================================
{
	private:

	protected:

	public :
        std::string name;
		int number;
		int ndata;
		int ninfo;
		int ncolumns;
        arma::Col<int> c_data;
        arma::mat data;
        int skiplines;
		
		opti_data(); 	//default constructor
		opti_data(int, int);	//constructor - allocates memory for statev
        opti_data(std::string, int, int, int, int, int); //Constructor with parameters
		opti_data(const opti_data &);	//Copy constructor
		~opti_data();
		
		int dimdata () const {return ndata;}       // returns the number of data points
		int diminfo () const {return ninfo;}       // returns the number of informations at each datapoint

		void constructc_data();
		void constructdata();

        void import(std::string, int=0);
				
		virtual opti_data& operator = (const opti_data&);
		
        friend  std::ostream& operator << (std::ostream&, const opti_data&);
};

} //namespace simcoon
