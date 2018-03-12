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

///@file opti_data.cpp
///@brief object that stores exp/num data
///@author Chemisky & Despringre
///@version 1.0
///@date 05/28/2014

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/Simulation/Identification/opti_data.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for phases_characteristics===================================

//=====Public methods for phases_characteristics============================================

//@brief default constructor
//-------------------------------------------------------------
opti_data::opti_data()
//-------------------------------------------------------------
{
	name="undefined";
	number=0;
	ninfo=0;
	ndata=0;
	ncolumns=0;
    skiplines = 0;
}

/*!
  \brief Constructor
  \param n : number of data points
  \param m : number of information in each point 
*/

//-------------------------------------------------------------
opti_data::opti_data(int n, int m)
//-------------------------------------------------------------
{
	assert(n>0);
	assert(m>0);

	name="undefined";
	number = 0;
	ndata = n;
	ninfo = m;
	ncolumns=0;
    skiplines = 0;
	
	c_data.zeros(m);
	data = zeros(n,m);
}

/*!
  \brief Constructor with parameters
  \param mname : name of the experimental data file (with the extension)
  \param mnumber : number of the experimental file
  \param mndata : number of data points
*/

//-------------------------------------------------------------
opti_data::opti_data(string mname, int mnumber, int mninfo, int mndata, int mncolumns, int mskiplines)
//-------------------------------------------------------------
{	
	assert(mndata>0);
	assert(mninfo>0);	

	name=mname;
	number = mnumber;
	ndata = mndata;
	ninfo = mninfo;
	ncolumns = mncolumns;
    skiplines = mskiplines;
		
	c_data.zeros(mninfo);
	data = zeros(mndata,mninfo);
}

/*!
  \brief Copy constructor
  \param ed opti_data object to duplicate
*/

//------------------------------------------------------
opti_data::opti_data(const opti_data& ed)
//------------------------------------------------------
{
	assert(ed.ndata>0);
	assert(ed.ninfo>0);	
  
	name=ed.name;
	number = ed.number;
	ndata = ed.ndata;
	ninfo = ed.ninfo;
	ncolumns = ed.ncolumns;
    skiplines = ed.skiplines;

	c_data = ed.c_data;
	data = ed.data;
}

/*!
  \brief destructor
*/

opti_data::~opti_data() {}

//-------------------------------------------------------------
void opti_data::constructc_data()
//-------------------------------------------------------------
{
	assert(ninfo>0);
	
	c_data.zeros(ninfo);
}

//-------------------------------------------------------------
void opti_data::constructdata()
//-------------------------------------------------------------
{
	assert(ninfo>0);
	assert(ndata>0);
	
	data = zeros(ndata,ninfo);
}

//-------------------------------------------------------------
void opti_data::import(string folder, int nexp)
//-------------------------------------------------------------
{
    assert(ninfo>0);
    assert(ncolumns>0);
    if (nexp > 0) {
        ndata = nexp;
        constructdata();
    }
    ndata = 0;
    ifstream ifdata;
    string buffer;
    
    string temp_string;
    string path = folder + "/" + name;
    
    ifdata.open(path, ios::in);
    while (!ifdata.eof())
    {
        getline (ifdata,buffer);
        if (buffer != "") {
            ndata++;
        }
    }
    ndata -= skiplines;
    
    ifdata.close();
    
    ifdata.open(path, ios::in);
    
    for (int i=0; i<skiplines; i++) {
        getline (ifdata,buffer);
    }
    
    if ((nexp == 0)||(ndata < nexp)) {
        constructdata();
        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < ncolumns; j++) {
                ifdata >> temp_string;
                for(int k=0; k<ninfo; k++) {
                    if((j)==c_data(k)) {
                        data(i, k)=stod(temp_string);
                    }
                }				
            }
        }
    }
    else if (ndata > nexp) {
        for (int i = 0; i < nexp; i++) {
            for (int j = 0; j < ncolumns; j++) {
                ifdata >> temp_string;
                for(int k=0; k<ninfo; k++) {
                    if((j)==c_data(k)) {
                        data(i, k)=stod(temp_string);
                    }
                }
            }
        }
    }
    
    ifdata.close();
    ifdata.clear();
}

/*!
  \brief Standard operator = for opti_data
*/

//----------------------------------------------------------------------
opti_data& opti_data::operator = (const opti_data& ed)
//----------------------------------------------------------------------
{
	assert(ed.ndata>0);
	assert(ed.ninfo>0);	
  
	name=ed.name;
	number = ed.number;
	ndata = ed.ndata;
	ninfo = ed.ninfo;
	ncolumns = ed.ncolumns;	

	c_data = ed.c_data;
	data = ed.data;

	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const opti_data& ed)
//--------------------------------------------------------------------------
{
//	assert(ed.ndata>0);
//	assert(ed.ninfo>0);	

	s << "Display info on the experimental data\n";
	s << "Name of the experimental data file: " << ed.name << "\n";
	s << "Number of the experimental data file: " << ed.number << "\n";
	s << "Number of data points: " << ed.ndata << "\n";
	s << "Number of informations per data point: " << ed.ninfo << "\n";
	s << "Number of columns in the file: " << ed.ncolumns << "\n";
	s << "Number of lines skipped at the beginning of the file: " << ed.skiplines << "\n";
    
/*	for (int i=1; i<=ed.ndata; i++) {
	  
		s << i;
		for (int j=1; j<=ed.ninfo; j++) {
		      s << "\t" << ed.E(i,j);
		}
		s << "\n";	
	}

	s << "\n\n";*/

	return s;
}
    
} //namespace simcoon
