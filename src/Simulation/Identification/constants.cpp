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

///@file constants.cpp
///@brief Handle of input constants
///@version 0.9

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Identification/constants.hpp>

using namespace std;
using namespace arma;

namespace simcoon{
    
//=====Private methods for constants===================================

//=====Public methods for constants============================================

//@brief default constructor
//-------------------------------------------------------------
constants::constants()
//-------------------------------------------------------------
{
    
    number = 0;
    value = 0.;  //Initial Value of the parameter
    ninput_files = 0;
}

/*!
 \brief Constructor
 \param mnumber : number of the constant
 \param nfiles : Number of files where the constant is present
 */

//-------------------------------------------------------------
constants::constants(const int &mnumber, const int& nfiles)
//-------------------------------------------------------------
{
    assert(nfiles > 0);
    
    number = mnumber;
    input_values.resize(nfiles);
    value = 0.;
    ninput_files = 0;
}

/*!
 \brief Constructor with parameters
 \param mnumber : number of the constant
 */

//-------------------------------------------------------------
constants::constants(const int &mnumber, const double &mvalue, const vec &minput_values, const string &mkey, const int &mninput_files, const std::vector<std::string> &minput_files)
//-------------------------------------------------------------
{
    number = mnumber;
    value = mvalue;
    input_values = minput_values;
    key = mkey;
    ninput_files = mninput_files;
    input_files = minput_files;
}

/*!
 \brief Copy constructor
 \param ed opti_data object to duplicate
 */

//------------------------------------------------------
constants::constants(const constants& co)
//------------------------------------------------------
{
    number=co.number;
    value = co.value;
    input_values = co.input_values;
    key = co.key;
    ninput_files = co.ninput_files;
    input_files = co.input_files;
}

/*!
 \brief destructor
 */

constants::~constants() {}

//-------------------------------------------------------------
void constants::update(const int &file)
//-------------------------------------------------------------
{
    value = input_values(file);
}
    
//-------------------------------------------------------------
void constants::resize(const int &n, const int &nfiles)
//-------------------------------------------------------------
{
    ninput_files = n;
    input_files.resize(n);
    input_values.resize(nfiles);
}
    
//----------------------------------------------------------------------
constants& constants::operator = (const constants& co)
//----------------------------------------------------------------------
{
    number=co.number;
    value = co.value;
    input_values = co.input_values;
    key = co.key;
    ninput_files = co.ninput_files;
    input_files = co.input_files;
    
    return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const constants& co)
//--------------------------------------------------------------------------
{
    
    s << "Display info on the parameter data\n";
    s << "Number of the parameter: " << co.number << "\n";
    s << "Number of files impacted and list of files: " << co.input_files.size() << "\n";
    
    for (vector<string>::const_iterator iter = co.input_files.begin(); iter != co.input_files.end(); iter++) {
        cout << *iter << "\n";
    }
    
    return s;
}
    
} //namespace simcoon