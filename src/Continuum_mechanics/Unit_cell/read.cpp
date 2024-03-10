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
#include <simcoon/Continuum_mechanics/Unit_cell/element.hpp>
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
    
void unit_cell_essentials(unsigned int &simul_type, unsigned int &loading_type, unsigned int &BC_type, unsigned int &max_temp, const string &path, const string &filename) {
    
    string pathfile = path + "/" + filename;
    ifstream unit_cell_essentials;
    string buffer;
    
    unit_cell_essentials.open(pathfile, ios::in);
    if(!unit_cell_essentials) {
        cout << "Error: cannot open : " << filename << " in :" << path << endl;
        return;
    }
    ///Get the control values for the solver
    unit_cell_essentials >> buffer >> simul_type >> buffer >> loading_type >> buffer >> BC_type >> buffer >> max_temp;
    unit_cell_essentials.close();
}
    
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

void read_abaqus_full_mesh(std::vector<section_characteristics> &sc, std::vector<Node> &nodes_full, std::vector<Element> & elements_full, std::vector<std::vector<std::string> > &sc_mesh, const unsigned int &loading_type, const std::vector<std::string> &elements_type, const string &path_data, const std::string &inputfile_sections, const std::string &inputfile_sc_def, const std::string &inputfile_mesh) {
   
    std::string buffer;
    //Read the sections first
    read_sections(sc, loading_type, path_data, inputfile_sections);
    
    //Read the sc_vector next
    std::string path_sc_mesh = path_data + "/" + inputfile_sc_def;
    std::ifstream paramphases;
    vector<string> strings;
    paramphases.open(path_sc_mesh, ios::in);
    if(paramphases) {
        while(getline(paramphases, buffer, '\n')) {
            if (buffer[buffer.size() - 1] == '\r')
            buffer.resize(buffer.size() - 1);
            
            std::vector<std::string> object;
            if (buffer != "") {
                strings = split(buffer, ',');
                for (auto s:strings) {
                    object.push_back(s);
                }
            }
            sc_mesh.push_back(object);
        }
    }
        
/*    for (auto object_list:sc_mesh) {
        for (auto s:object_list) {
            cout << s << "\t";
        }
        cout << endl;
    }
*/
    
    paramphases.close();
    assert(sc_mesh.size() == sc.size());
    
    //open the reading file
    std::string path_mesh = path_data + "/" + inputfile_mesh;
    std::ifstream para_mesh;
    para_mesh.open(path_mesh, ios::in);
    int Node_temp = 0;
    int Element_temp = 0;
    bool badInput = false;

    std::vector<unsigned int> NSet_number;
    std::vector<unsigned int> ElSet_number;
    std::string el_type_found;
    std::string nset_name_found;
    std::string elset_name_found;
    std::vector<std::string> search;
    std::map<std::string, bool> NSet_is_generate;
    std::map<std::string, bool> ElSet_is_generate;

    bool Nodefound = false;
    bool Elementfound = false;
    bool NSetfound = false;
    bool ElSetfound = false;
    bool switch_on = false;
    bool is_generate = false;
    
    //Node preparation
    vec coords = zeros(3);
    nodes_full.clear();
    
    if (para_mesh) {
        while(getline(para_mesh, buffer, '\n')) {
            if (buffer[buffer.size() - 1] == '\r') {
                buffer.resize(buffer.size() - 1);
                std::transform(buffer.begin(), buffer.end(), buffer.begin(),[](unsigned char c){ return std::tolower(c); });
            }
            
            switch_on = false;
            badInput = false;
            is_generate = false;

            if (buffer.find("*", 0) != string::npos) {
                Nodefound = false;
                Elementfound = false;
                NSetfound = false;
                ElSetfound = false;
            }
            
            search.push_back("*Node");
            search.push_back("*node");
            search.push_back("*NODE");
            for(auto s:search) {
                if (buffer.find(s, 0) != string::npos) {
                    Nodefound = true;
                    Elementfound = false;
                    NSetfound = false;
                    ElSetfound = false;

                    switch_on = true;
                }
            }
            search.clear();
            
            for (auto el_type:elements_type) {
                search.clear();
                search.push_back("*element, type=" + el_type);
                search.push_back("*Element, type=" + el_type);
                search.push_back("*Element, Type=" + el_type);
                search.push_back("*ELEMENT, type=" + el_type);
                for(auto s:search) {
                    if (buffer.find(s, 0) != string::npos) {
                        Nodefound = false;
                        Elementfound = true;
                        NSetfound = false;
                        ElSetfound = false;
                        
                        switch_on = true;
                        el_type_found = el_type;
                    }
                }
            }
            search.clear();

            search.push_back("*NSet, NSet=");
            search.push_back("*NSet,NSet=");
            search.push_back("*NSET, NSET=");
            search.push_back("*NSET,NSET=");
            search.push_back("*Nset, nset=");
            for(auto s:search) {
                if (buffer.find(s, 0) != string::npos) {
                    Nodefound = false;
                    Elementfound = false;
                    NSetfound = true;
                    ElSetfound = false;
                    
                    switch_on = true;
                    strings = split(buffer, ',');
                    strings = split(strings[1], '=');
                    nset_name_found = strings[1];
//                cout << "string::npos : " << string::npos << endl;
                }
            }
            search.clear();
            
            search.push_back("*ElSet, ElSet=");
            search.push_back("*ElSet,ElSet=");
            search.push_back("*ELSET, ELSET=");
            search.push_back("*ELSET,ELSET=");
            search.push_back("*Elset, elset=");
            for(auto s:search) {
                if (buffer.find(s, 0) != string::npos) {
                    Nodefound = false;
                    Elementfound = false;
                    NSetfound = false;
                    ElSetfound = true;
                    
                    switch_on = true;
                    strings = split(buffer, ',');
                    strings = split(strings[1], '=');
                    elset_name_found = strings[1];
                }
            }
            search.clear();

            search.push_back("generate");
            search.push_back("Generate");
            search.push_back("GENERATE");
            for(auto s:search) {
                if (buffer.find(s, 0) != string::npos) {
                    is_generate = true;
                }
            }
            search.clear();

            search.push_back("output");
            search.push_back("Output");
            search.push_back("OUTPUT");
            for(auto s:search) {
                if (buffer.find(s, 0) != string::npos) {
                    Nodefound = false;
                    Elementfound = false;
                    NSetfound = false;
                    ElSetfound = false;
                }
            }
            search.clear();

            if (buffer == "") {
                Nodefound = false;
                Elementfound = false;
                NSetfound = false;
                ElSetfound = false;
                switch_on = true;
            }
            search.clear();
                        
            if((Nodefound)&&(switch_on==false)) {
                strings = split(buffer, ',');
                coords(0) = stod(strings[1]);
                coords(1) = stod(strings[2]);
                coords(2) = stod(strings[3]);
                Node temp_node(stoi(strings[0]), Point(coords(0),coords(1),coords(2)));
                nodes_full.push_back(temp_node);
//                cout << "pushed node : " << temp_node << endl;
            }
            if((Elementfound)&&(switch_on==false)) {
                strings = split(buffer, ',');
                Element temp_el;
                temp_el.number = stoi(strings[0]);
                temp_el.type = el_type_found;
                std::vector<Node> temp_el_node_list;
                unsigned int number_nodes_per_el = strings.size()-1;
                for(unsigned int j=0; j<number_nodes_per_el; j++) {
                    auto it_node_inlist = std::find(nodes_full.begin(), nodes_full.end(), Node(stoi(strings[j+1]),Point(0.,0.,0.)));
                    Node node_inlist = *it_node_inlist;
                    temp_el_node_list.push_back(node_inlist);
                }
                temp_el.nodes = temp_el_node_list;
                elements_full.push_back(temp_el);
//                cout << "pushed element : " << temp_el << endl;
            }
            if(NSetfound) {
                if (switch_on==true) {
                    for(unsigned int i=0; i<sc_mesh.size(); i++) {
                        for(unsigned int j=0; j<sc_mesh[i].size(); j++) {
                            if (sc_mesh[i][j] == nset_name_found) {
                                NSet_is_generate.insert ( std::pair<std::string,bool>(nset_name_found,is_generate) );
                            }
                        }
                    }
                }
                else {
    //cout << "found nset : " << NSetfound << endl;
                    strings = split(buffer, ',');
                    for(unsigned int i=0; i<sc_mesh.size(); i++) {
                        for(unsigned int j=0; j<sc_mesh[i].size(); j++) {
                            if (sc_mesh[i][j] == nset_name_found) {
                                if (NSet_is_generate[nset_name_found] == false) {
                                    unsigned int number_nodes_per_line = strings.size();
                                    for(unsigned int k=0; k<number_nodes_per_line; k++) {
                                        try {
                                            Node_temp = stoi(strings[k]);
                                        }
                                        catch (...) {
                                            badInput = true;
                                        }
                                        if (badInput == false) {
                                            auto it_node_inlist = std::find(nodes_full.begin(), nodes_full.end(), Node(Node_temp,Point(0.,0.,0.)));
                                            sc[i].nodes.push_back(*it_node_inlist);
                                        }
                                    }
                                }
                                else {
                                    int begin = stoi(strings[0]);
                                    int end = stoi(strings[1]);
                                    int pas = stoi(strings[2]);
                                    for(unsigned int k=begin; k<=end; k+=pas) {
                                        auto it_node_inlist = std::find(nodes_full.begin(), nodes_full.end(), Node(k,Point(0.,0.,0.)));
                                        sc[i].nodes.push_back(*it_node_inlist);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if((ElSetfound)&&(switch_on==false)) {
                if (switch_on==false) {
                    for(unsigned int i=0; i<sc_mesh.size(); i++) {
                        for(unsigned int j=0; j<sc_mesh[i].size(); j++) {
                            if (sc_mesh[i][j] == elset_name_found) {
                                ElSet_is_generate.insert ( std::pair<std::string,bool>(elset_name_found,is_generate) );
                            }
                        }
                    }
                }
                else {
                    strings = split(buffer, ',');
                    for(unsigned int i=0; i<sc_mesh.size(); i++) {
                        for(unsigned int j=0; j<sc_mesh[i].size(); j++) {
                            if (sc_mesh[i][j] == elset_name_found) {
                                if (is_generate == false) {
                                    unsigned int number_elements_per_line = strings.size();
                                    for(unsigned int k=0; k<number_elements_per_line; k++) {
                                        try {
                                            Element_temp = stoi(strings[k]);
                                        }
                                        catch (...) {
                                            badInput = true;
                                        }
                                        if (badInput == false) {
                                            std::vector<Node> null;
                                            string el_type_temp;
                                            auto it_element_inlist = std::find(elements_full.begin(), elements_full.end(), Element(el_type_temp,Element_temp,null));
                                            sc[i].elements.push_back(*it_element_inlist);
                                        }
                                    }
                                }
                                else {
                                    int begin = stoi(strings[0]);
                                    int end = stoi(strings[1]);
                                    int pas = stoi(strings[2]);
                                    for(unsigned int k=begin; k<=end; i+=pas) {
                                        std::vector<Node> null;
                                        string el_type_temp;
                                        auto it_element_inlist = std::find(elements_full.begin(), elements_full.end(), Element(el_type_temp,k,null));
                                        sc[i].elements.push_back(*it_element_inlist);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        cout << "Error: The file that describe the mesh (mesh.inp) has not been found" << endl;
    }
    para_mesh.close();
    
}

void read_salome_full_mesh(std::vector<section_characteristics> &sc, std::vector<Col<int> > &sc_mesh, const unsigned int &loading_type, const string &path_data, const string &inputfile_sections, const string &inputfile_sc_def, const string &inputfile_mesh) {
   
    std::string buffer;
    std::ifstream para_mesh;
    
    //Read the sections first
    read_sections(sc, loading_type, path_data, inputfile_sections);
    
    //Read the sc_vector next
    std::string path_sc_mesh = path_data + "/" + inputfile_sc_def;
    std::ifstream paramphases;
    vector<string> strings;
    paramphases.open(path_sc_mesh, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                strings = split(buffer, ',');
                Col<int> temp(strings.size());
                for (unsigned int i=0; i<strings.size(); i++) {
                    temp(i) = stoi(strings[i]);
                }
            sc_mesh.push_back(temp);
            }
        }
    }
    assert(sc_mesh.size() == sc.size());
    paramphases.close();

    //Read the vector of nodes
    std::vector<Node> nodes_temp;
    std::vector<Element> elements_temp;
    
    for(unsigned int i=0; i<sc_mesh.size(); i++) {
        for(unsigned int j=0; j<sc_mesh[i].n_elem; j++) {
            std::string inputfile_mesh_ext = inputfile_mesh.substr(inputfile_mesh.length()-4,inputfile_mesh.length());
            std::string inputfile_mesh_name_root = inputfile_mesh.substr(0,inputfile_mesh.length()-4); //to remove the extension
            std::string inputfile_ij = inputfile_mesh_name_root + '_' + to_string(sc_mesh[i](j)) + inputfile_mesh_ext;

            if (nodes_temp.size()>0) {
                nodes_temp.clear();
            }
            if (elements_temp.size()>0) {
                elements_temp.clear();
            }
            read_salome_mesh(nodes_temp, elements_temp, path_data, inputfile_ij);

            for(auto n: nodes_temp) {
                sc[i].nodes.push_back(n);
            }
            for(auto e: elements_temp) {
                sc[i].elements.push_back(e);
            }
        }
    }
}

void read_salome_mesh(std::vector<Node> &nodes, std::vector<Element> &elements, const string &path_data, const string &inputfile) {

    unsigned int number_nodes = 0;
    unsigned int number_elements = 0;

    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;

    std::ifstream para_mesh;
    para_mesh.open(path_inputfile, ios::in);
    vector<string> strings;
    if(para_mesh) {
        getline (para_mesh,buffer);
        if (buffer != "") {
            strings = split(buffer, ' ');
            number_nodes = stoi(strings[0]);
            number_elements = stoi(strings[1]);
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the node list in the folder :" << path_data << endl;
        return;
    }
    para_mesh.close();

    //open the reading file
    para_mesh.open(path_inputfile, ios::in);
    vec coords = zeros(3);
    if(para_mesh) {
        getline (para_mesh,buffer);
        for(unsigned int i=0; i<number_nodes; i++) {
            getline(para_mesh,buffer);
            if (buffer != "") {
                strings = split(buffer, ' ');
                coords(0) = stod(strings[1]);
                coords(1) = stod(strings[2]);
                coords(2) = stod(strings[3]);
                Node temp_node(stoi(strings[0]), Point(coords(0),coords(1),coords(2)));
                nodes.push_back(temp_node);
            }
        }
        for(unsigned int i=0; i<number_elements; i++) {
            getline(para_mesh,buffer);
            if (buffer != "") {
                strings = split(buffer, ' ');
                //edges
/*                if(strings[1].front() == '1') {
                   elements[i].number = stoi(strings[0]);
                   elements[i].type = strings[1];
                   unsigned int number_nodes_per_el = stoi(strings[1].substr(1,strings[1].length()));
                   for(unsigned int j=0; j<number_nodes_per_el; j++) {
                      auto it_node_inlist = std::find(nodes.begin(), nodes.end(), Node(stoi(strings[j+2]),Point(0.,0.,0.)));
                      Node node_inlist = *it_node_inlist;
                      elements[i].nodes.push_back(node_inlist);
                   }
                
                //faces
                if(strings[1].front() == '2') {
                   elements[i].number = stoi(strings[0]);
                   elements[i].type = strings[1];
                   unsigned int number_nodes_per_el = stoi(strings[1].substr(1,strings[1].length()));
                   for(unsigned int j=0; j<number_nodes_per_el; j++) {
                      auto it_node_inlist = std::find(nodes.begin(), nodes.end(), Node(stoi(strings[j+2]),Point(0.,0.,0.)));
                      Node node_inlist = *it_node_inlist;
                      elements[i].nodes.push_back(node_inlist);
                   }*/
                
                //volumes
                if(strings[1].front() == '3') {
                   Element temp_el;
                   temp_el.number = stoi(strings[0]);
                   temp_el.type = strings[1];
                   unsigned int number_nodes_per_el = stoi(strings[1].substr(1,strings[1].length()));
                   for(unsigned int j=0; j<number_nodes_per_el; j++) {
                      auto it_node_inlist = std::find(nodes.begin(), nodes.end(), Node(stoi(strings[j+2]),Point(0.,0.,0.)));
                      Node node_inlist = *it_node_inlist;
                      temp_el.nodes.push_back(node_inlist);
                   }
                   elements.push_back(temp_el);
                }
                
            }
        }
    //end of reading
    }
    para_mesh.close();
}

void read_nodes_file(std::vector<Node> &nodes, const string &path_data, const string &inputfile) {
    
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

void read_results(vec &E, vec &S, const string &path_data, const string &resultslfile) {
  string buffer;
  ifstream res;
  string path_res = path_data + "/" + resultslfile;
  res.open(path_res, ios::in);
  E = zeros(6);
  S = zeros(6);
  if (res) {
    res >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    for (unsigned int i=0;i<6;i++) {
      res >> E(i);
    }
    for (unsigned int i=0;i<105;i++) {
      res >> buffer;
    }
    for (unsigned int i=0;i<6;i++) {
      res >> S(i);
    }
  }
  else {
    cout << "Error: cannot open the file " << resultslfile << " in the folder :" << path_data << endl;
  }
  res.close();
}
    
void read_subphases_results(vec &E, vec &S, const unsigned int &ph, const string &path_data, const string &resultslfile) {
  string buffer;
  ifstream res;
  string path_res = path_data + "/" + resultslfile;
  res.open(path_res, ios::in);
  E = zeros(6);
  S = zeros(6);
  if (res) {
    for (unsigned int i=0;i<17+ph*8;i++) {
      res >> buffer;
    }
    for (unsigned int i=0;i<6;i++) {
      res >> E(i);
    }
    for (unsigned int i=0;i<105;i++) {
      res >> buffer;
    }
    for (unsigned int i=0;i<6;i++) {
      res >> S(i);
    }
  }
  else {
    cout << "Error: cannot open the file " << resultslfile << " in the folder :" << path_data << endl;
  }
  res.close();
}

} //namespace simcoon
