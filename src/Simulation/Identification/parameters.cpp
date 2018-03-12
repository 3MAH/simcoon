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

///@file parameters.cpp
///@brief Handle of input parameters
///@version 0.9

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Identification/parameters.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for parameters===================================

//=====Public methods for parameters============================================

//@brief default constructor
//-------------------------------------------------------------
parameters::parameters()
//-------------------------------------------------------------
{
    
    number = 0;
    value = 0.;  //Initial Value of the parameter
    min_value = 0.;   //Minimum value of the parameter
    max_value = 0.;   //Maximum value of the parameter
    ninput_files = 0;
}

/*!
 \brief Constructor
 \param mnumber : number of the parameter
 \param mmin_value : Minimal value of the parameter
 \param mmax_value : Maximal value of the parameter
 */

//-------------------------------------------------------------
parameters::parameters(const int &mnumber, const double &mmin_value, const double &mmax_value)
//-------------------------------------------------------------
{
    number = mnumber;
    min_value = mmin_value;
    max_value = mmax_value;
    value = (min_value+max_value)/2.;
    ninput_files = 0;
}

/*!
 \brief Constructor with parameters
 \param mnumber : number of the parameter
 \param mmin_value : Minimal value of the parameter
 \param mmax_value : Maximal value of the parameter
 */

//-------------------------------------------------------------
parameters::parameters(const int &mnumber, const double &mmin_value, const double &mmax_value, const string &mkey, const int &mninput_files, const std::vector<std::string> &minput_files)
//-------------------------------------------------------------
{
    number = mnumber;
    min_value = mmin_value;
    max_value = mmax_value;
    value = (min_value+max_value)/2.;

    key = mkey;
    ninput_files = mninput_files;
    input_files = minput_files;
}

/*!
 \brief Copy constructor
 \param ed opti_data object to duplicate
 */

//------------------------------------------------------
parameters::parameters(const parameters& pa)
//------------------------------------------------------
{
    number=pa.number;
    min_value = pa.min_value;
    max_value = pa.max_value;
    value = pa.value;
    
    key = pa.key;
    ninput_files = pa.ninput_files;
    input_files = pa.input_files;
}

/*!
 \brief destructor
 */

parameters::~parameters() {}

//-------------------------------------------------------------
void parameters::update(const double &p)
//-------------------------------------------------------------
{
    value = p;
}
    
    
//-------------------------------------------------------------
void parameters::resize(const int &n)
//-------------------------------------------------------------
{
    
    ninput_files = n;
    input_files.resize(n);
}
    
//----------------------------------------------------------------------
parameters& parameters::operator = (const parameters& pa)
//----------------------------------------------------------------------
{
    number=pa.number;
    min_value = pa.min_value;
    max_value = pa.max_value;
    value = pa.value;
    
    key = pa.key;
    ninput_files = pa.ninput_files;
    input_files = pa.input_files;
    
    return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const parameters& pa)
//--------------------------------------------------------------------------
{
    
    s << "Display info on the parameter data\n";
    s << "Number of the parameter: " << pa.number << "\n";
    s << "Bounds (Min and Max) values: " << pa.min_value << "\t" << pa.max_value << "\n";
    s << "Number of files impacted and list of files: " << pa.input_files.size() << "\n";
    
    for (vector<string>::const_iterator iter = pa.input_files.begin(); iter != pa.input_files.end(); iter++) {
        cout << *iter << "\n";
    }
    
    return s;
}
    
} //namespace simcoon