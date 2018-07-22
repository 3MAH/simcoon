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
///@brief To read from NsectionsX.dat
///@version 1.0

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/read.hpp>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

using namespace std;
using namespace arma;

namespace simcoon{

std::vector<std::string> split(const std::string &s, const char &delim)
{
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}
    
/*void unit_cell_essentials(unsigned int &loading_type, unsigned int &BC_type, int &max_temp, const string &path, const string &filename) {
    
    string pathfile = path + "/" + filename;
    ifstream unit_cell_essentials;
    string buffer;
    
    unit_cell_essentials.open(pathfile, ios::in);
    if(!unit_cell_essentials) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
    
    ///Get the control values for the solver
    unit_cell_essentials >> buffer >> loading_type >> buffer >> BC_type >> buffer >> max_temp;
    unit_cell_essentials.close();
}*/
    
void read_sections(std::vector<section_characteristics> &sections, const unsigned int &loading_type, const string &path_data, const string &inputfile) {
    
    unsigned int nsections = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nsections++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the phase characteristics in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nsections--;
    
    //Resize the vector
    sections.resize(nsections);
    int nprops = 0;
    int nstatev = 0;
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(unsigned int i=0; i<nsections; i++) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nprops >> nstatev;
        

        sections[i].abamat.resize(nprops, nstatev);
        sections[i].abamat.update(0, 0, "umat", 1, 0., 0., 0., 0., 0., nprops, nstatev, zeros(nprops));
        
        for(int j=0; j<sections[i].abamat.nprops; j++) {
            paramphases >> buffer;
        }
    }
    paramphases.close();
    
    paramphases.open(path_inputfile, ios::in);
    if (loading_type == 1) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    }
    else if((loading_type == 2)||(loading_type == 3)){
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    }
    else
        cout << "error in Continuum_mechanics/Unit_cell/read : loading type should take values in the range (1,3)" << endl;

    
    for(unsigned int i=0; i<nsections; i++) {
        
        if (loading_type == 1) {
            paramphases >> sections[i].abamat.number >> sections[i].elset_name >> sections[i].abamat.umat_name >> sections[i].abamat.psi_mat >> sections[i].abamat.theta_mat >> sections[i].abamat.phi_mat >> buffer >> buffer;
        }
        else if((loading_type == 2)||(loading_type == 3)){
            paramphases >> sections[i].abamat.number >> sections[i].elset_name >> sections[i].abamat.umat_name >> sections[i].abamat.psi_mat >> sections[i].abamat.theta_mat >> sections[i].abamat.phi_mat >> sections[i].abamat.conductivity >> sections[i].abamat.density >> buffer >> buffer;
        }
        
        sections[i].abamat.id = sections[i].abamat.number;
        
        for(int j=0; j<sections[i].abamat.nprops; j++) {
            paramphases >> sections[i].abamat.props(j);
        }
        
        sections[i].abamat.psi_mat*=(sim_pi/180.);
        sections[i].abamat.theta_mat*=(sim_pi/180.);
        sections[i].abamat.phi_mat*=(sim_pi/180.);
    }
    
    paramphases.close();
}

void read_mesh(std::vector<Node> &nodes, const string &path_data, const string &inputfile) {
    
    unsigned int nnodes = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nnodes++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the node list in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nnodes--;
    
    //Resize the vector
    nodes.resize(nnodes);
    paramphases.open(path_inputfile, ios::in);
    vec coords = zeros(3);

    getline(paramphases,buffer);
    vector<string> strings;
    if(paramphases) {
        for(unsigned int i=0; i<nnodes; i++) {
            getline(paramphases,buffer);
            if (buffer != "") {
                strings = split(buffer, ',');
                nodes[i].number = stoi(strings[0]);
                coords(0) = stod(strings[1]);
                coords(1) = stod(strings[2]);
                coords(2) = stod(strings[3]);
                Point p(coords(0),coords(1),coords(2));
                nodes[i].coords = p;
            }
        }
    }
    paramphases.close();
}
    

} //namespace simcoon
