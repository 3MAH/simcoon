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

///@file read.cpp
///@brief To read ODF or PDF from Npeak.dat
///@version 1.0

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_Mechanics/Material/ODF.hpp>
#include <simcoon/Continuum_Mechanics/Material/PDF.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void read_peak(ODF &odf_rve, const string &path_data, const string &inputfile) {
    
    unsigned int npeaks = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramodf;
    
    paramodf.open(path_inputfile, ios::in);
    if(paramodf) {
        while (!paramodf.eof())
        {
            getline (paramodf,buffer);
            if (buffer != "") {
                npeaks++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the peak characteristics in the folder: " << path_data << endl;
        return;
    }
    paramodf.close();
    npeaks--;
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    odf_rve.construct(npeaks);
    
    int nparams = 0;
    
    paramodf.open(path_inputfile, ios::in);
    paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nparams;
        
        odf_rve.peaks[i].params = zeros(nparams);
        for(unsigned int j=0; j<odf_rve.peaks[i].params.n_elem; j++) {
            paramodf >> buffer;
        }
    }
    paramodf.close();
    
    paramodf.open(path_inputfile, ios::in);
    paramodf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        paramodf >> odf_rve.peaks[i].number >> odf_rve.peaks[i].method >> odf_rve.peaks[i].mean >> odf_rve.peaks[i].s_dev >> odf_rve.peaks[i].width >> odf_rve.peaks[i].ampl >> buffer;
        
        for(unsigned int j=0; j<odf_rve.peaks[i].params.n_elem; j++) {
            paramodf >> odf_rve.peaks[i].params(j);
        }
        
        if(odf_rve.radian == false) {
            odf_rve.peaks[i].mean *=(pi/180.);
            odf_rve.peaks[i].s_dev *=(pi/180.);
            odf_rve.peaks[i].width *=(pi/180.);
        }
    }
    
    paramodf.close();
}



//Overload of the read_peak function for PDF
void read_peak(PDF &pdf_rve, const string &path_data, const string &inputfile) {
    
    unsigned int npeaks = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream parampdf;
    
    parampdf.open(path_inputfile, ios::in);
    if(parampdf) {
        while (!parampdf.eof())
        {
            getline (parampdf,buffer);
            if (buffer != "") {
                npeaks++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the peak characteristics in the folder: " << path_data << endl;
        return;
    }
    parampdf.close();
    npeaks--;
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    pdf_rve.construct(npeaks);
    
    int nparams = 0;
    
    parampdf.open(path_inputfile, ios::in);
    parampdf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        parampdf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nparams;
        
        pdf_rve.peaks[i].params = zeros(nparams);
        for(unsigned int j=0; j<pdf_rve.peaks[i].params.n_elem; j++) {
            parampdf >> buffer;
        }
    }
    parampdf.close();
    
    parampdf.open(path_inputfile, ios::in);
    parampdf >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<npeaks; i++) {
        
        parampdf >> pdf_rve.peaks[i].number >> pdf_rve.peaks[i].method >> pdf_rve.peaks[i].mean >> pdf_rve.peaks[i].s_dev >> pdf_rve.peaks[i].width >> pdf_rve.peaks[i].ampl >> buffer;
        
        for(unsigned int j=0; j<pdf_rve.peaks[i].params.n_elem; j++) {
            parampdf >> pdf_rve.peaks[i].params(j);
        }
        
    }
    
    parampdf.close();
}

} //namespace simcoon
