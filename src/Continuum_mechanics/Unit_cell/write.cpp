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

///@file write.cpp
///@brief To write NphasesX.dat and NlayerX.dat files
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <boost/filesystem.hpp>
#include <CGAL/Simple_cartesian.h>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/node.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/section_characteristics.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/materials.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_meca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/step_thermomeca.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_mesh.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/cubic_equation.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/interpolate.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/geom_functions.hpp>
#include <simcoon/Continuum_mechanics/Unit_cell/write.hpp>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;

using namespace std;
using namespace arma;

namespace simcoon{

void update_sections(section_characteristics &section_rve, const int &nsections, const int &nscale, const int &geom_type, const string &path_data)
{
    
    string umat_name_macro = "ABAPC";
    phase_characteristics rve;
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_rve = {(double)nsections,(double)nscale};
    
    rve.sptr_matprops->update(0, umat_name_macro, 1, psi_rve, theta_rve, phi_rve, props_rve.n_elem, props_rve);
    rve.construct(geom_type,1); //The rve is supposed to be mechanical only here
    string inputfile; //file # that stores the microstructure properties
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            inputfile = "Nphases" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_phase(rve, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
        case 2: {
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 3: {
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_cylinder(rve, path_data, inputfile);
            break;
        }
    }
    
    int id = 99999;
    section_rve.update(umat_name_macro, id, rve);
}

void write_section(section_characteristics &section_rve, const unsigned int &loading_type, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    param_mats.close();
    
    section_rve.abamat.write(loading_type,path_data,outputfile,true);
    
    param_mats.open(filename, ios::app);
    vec x = {1.0,0.,0.};
    vec y = {0.,1.0,0.};
    
    mat R = fillR(section_rve.abamat.psi_mat, section_rve.abamat.theta_mat, section_rve.abamat.phi_mat,true,"xzx");
    vec u = R*x;
    vec v = R*y;
    
    param_mats << "*Orientation, name=Ori-" << section_rve.abamat.number << endl;
    param_mats << u(0) << ", " << u(1) << ", " << u(2) << ", " << v(0) << ", " << v(1) << ", " << v(2) << endl;
    param_mats << "1, 0." << endl;
    param_mats << "*Solid Section, ElSet=" << section_rve.elset_name << ", orientation=Ori-" << section_rve.abamat.number << ", Material=" << section_rve.abamat.umat_name << "-" << section_rve.abamat.number << endl;
    param_mats.close();
}
    
    
void write_sections(section_characteristics &section_rve, const unsigned int &loading_type, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    param_mats.close();
    
    for(unsigned int i=0; i<section_rve.sub_sections.size(); i++) {
        section_rve.sub_sections[i].abamat.write(loading_type,path_data,outputfile,true);
    }
    
    param_mats.open(filename, ios::app);
    vec x = {1.0,0.,0.};
    vec y = {0.,1.0,0.};
    
    for(unsigned int i=0; i<section_rve.sub_sections.size(); i++) {
        
        mat R = fillR(section_rve.sub_sections[i].abamat.psi_mat, section_rve.sub_sections[i].abamat.theta_mat, section_rve.sub_sections[i].abamat.phi_mat,true,"zxz");
        vec u = R*x;
        vec v = R*y;
        
        param_mats << "*Orientation, name=Ori-" << section_rve.sub_sections[i].abamat.number << endl;
        param_mats << u(0) << ", " << u(1) << ", " << u(2) << ", " << v(0) << ", " << v(1) << ", " << v(2) << endl;
        param_mats << "1, 0." << endl;
        param_mats << "*Solid Section, ElSet=" << section_rve.sub_sections[i].elset_name << ", orientation=Ori-" << section_rve.sub_sections[i].abamat.number << ", Material=" << section_rve.sub_sections[i].abamat.umat_name << "-" << section_rve.sub_sections[i].abamat.number << endl;
    }
    param_mats.close();
}

void write_sections(std::vector<section_characteristics> &sections, const unsigned int &loading_type, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    param_mats.close();
    
    for(auto s : sections) {
        s.abamat.write(loading_type,path_data,outputfile,true);
    }
    
    param_mats.open(filename, ios::app);
    vec x = {1.0,0.,0.};
    vec y = {0.,1.0,0.};
    
    for(auto s : sections) {
        
        mat R = fillR(s.abamat.psi_mat, s.abamat.theta_mat, s.abamat.phi_mat,true,"zxz");
        vec u = R*x;
        vec v = R*y;

        param_mats << "*Orientation, name=Ori-" << s.abamat.number << endl;
        param_mats << u(0) << ", " << u(1) << ", " << u(2) << ", " << v(0) << ", " << v(1) << ", " << v(2) << endl;
        param_mats << "1, 0." << endl;
        param_mats << "*Solid Section, ElSet=" << s.elset_name << ", orientation=Ori-" << s.abamat.number << ", Material=" << s.abamat.umat_name << "-" << s.abamat.number << endl;
    }
    param_mats.close();
}
    
void update_materials(std::vector<aba_material> &aba_mats, const int &nphases, const int &nscale, const int &geom_type, const string &path_data) {

    aba_mats.clear();
    phase_characteristics rve;
    string umat_name_macro = "ABAPC";
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    vec props_rve = {(double)nphases,(double)nscale};
    
    rve.sptr_matprops->update(0, umat_name_macro, 1, psi_rve, theta_rve, phi_rve, props_rve.n_elem, props_rve);
    rve.construct(geom_type,1); //The rve is supposed to be mechanical only here
    string inputfile; //file # that stores the microstructure properties
    
    switch (geom_type) {
            
        case 0 : {
            //Definition from Nphases.dat
            inputfile = "Nphases" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_phase(rve, path_data, inputfile);
            break;
        }
        case 1: {
            //Definition from Nlayers.dat
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
        case 2: {
            //Definition from Nellipsoids.dat
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 3: {
            //Definition from Ncylinders.dat
            inputfile = "Ncylinders" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_cylinder(rve, path_data, inputfile);
            break;
        }
    }
    
    int id_mat = 0;
    
    for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
        aba_material temp;
        id_mat = 100000 + rve.sptr_matprops->props(1)*1000 + rve.sub_phases[i].sptr_matprops->number;
        temp.update(*rve.sub_phases[i].sptr_matprops, 0., 0., *rve.sub_phases[i].sptr_sv_global, id_mat);
        aba_mats.push_back(temp);
    }
    
}
    
void write_materials(std::vector<aba_material> &aba_mats, const unsigned int &loading_type, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);
    param_mats << "** ==================\n";
    param_mats << "** MATERIALS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
        
    for(unsigned int i=0; i<aba_mats.size(); i++) {
        aba_mats[i].write(loading_type,path_data,outputfile,true);
    }
    param_mats.close();
}
    
void update_steps(std::vector<std::shared_ptr<step> > &aba_steps, const std::vector<block> &blocks, const bool &nlgeom, const int &loading_type, const int &max_temp) {

    string name_step_ini = "Step";
    string name_step = name_step_ini;
    unsigned int type = 0;
    
    for (unsigned int i=0; i<blocks.size(); i++) {
        for (unsigned int j=0; j<blocks[i].ncycle; j++) {
            for (unsigned int k=0; k<blocks[i].nstep; k++) {
                
                switch (loading_type) {
                    case 1: {
                        
                        name_step = name_step_ini + to_string(i+1) + to_string(j+1) + to_string(k+1);
                        shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[k]);
                        if ((sptr_meca->Dn_init == 1.)&&(sptr_meca->Dn_mini == 1.)) {
                            type = 1;
                        }
                        else {
                            type = 0;
                        }
                        shared_ptr<aba_step_meca> sptr_meca_aba = make_shared<aba_step_meca>();
                        sptr_meca_aba->update(*sptr_meca, name_step, nlgeom, type);
                        aba_steps.push_back(sptr_meca_aba);
                        break;
                    }
                    case 2: {
                        
                        name_step = name_step_ini + to_string(i+1) + to_string(j+1) + to_string(k+1);
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[k]);
                        
                        if ((sptr_thermomeca->Dn_init == 1.)&&(sptr_thermomeca->Dn_mini == 1.)) {
                            type = 1;
                        }
                        else {
                            type = 0;
                        }
                        shared_ptr<aba_step_thermomeca> sptr_thermomeca_aba = make_shared<aba_step_thermomeca>();
                        sptr_thermomeca_aba->update(*sptr_thermomeca, name_step, nlgeom, type, max_temp);
                        aba_steps.push_back(sptr_thermomeca_aba);
                        break;
                    }
                    default: {
                        cout << "Please provide a consistent loading type for the block " << i << "\n";
                        break;
                    }
                }
            }
        }
    }
    
}
    
void write_steps(const std::vector<std::shared_ptr<step> > &aba_steps, const int &loading_type, const double &temp_init, const string &path_data, const string &outputfile) {
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream param_mats;
    
    param_mats.open(filename, ios::out);

    param_mats << "*Boundary\n";
    param_mats << "CentreNode, 1, 1\n";
    param_mats << "CentreNode, 2, 2\n";
    param_mats << "CentreNode, 3, 3\n";
    param_mats << "**\n";
    
    param_mats << "*Initial Conditions, type=TEMPERATURE\n";
    param_mats << "AllNodes, " << temp_init << "\n";
    param_mats << "**\n";
    
    param_mats << "** ==================\n";
    param_mats << "** STEPS\n";
    param_mats << "** ==================\n";
    param_mats << "**" << endl;
    
    switch (loading_type) {
        case 1: {
            
            for(unsigned int i=0; i<aba_steps.size(); i++) {
                shared_ptr<aba_step_meca> sptr_meca_aba = std::dynamic_pointer_cast<aba_step_meca>(aba_steps[i]);
                sptr_meca_aba->write(path_data,outputfile,true);
            }
            break;
        }
        case 2: {
            
            for(unsigned int i=0; i<aba_steps.size(); i++) {
                shared_ptr<aba_step_thermomeca> sptr_thermomeca_aba = std::dynamic_pointer_cast<aba_step_thermomeca>(aba_steps[i]);
                sptr_thermomeca_aba->write(path_data,outputfile,true);
            }
            break;
        }
        default: {
            cout << "Error in Continuum/Mechanics/Unit_Cell/write.cpp, function write_steps : Please provide a consistent loading type to write the list of steps in the .inp file" << endl;
            break;
        }
    }
    param_mats.close();
}

void write_nodes_file(const std::vector<Node> &nodes, std::ofstream &out_set){
    
    out_set << "*Node\n";
    for (auto n : nodes) {
        out_set << n;
    }
}
    
void write_node_set(const std::string &name, const Node &node, std::ofstream &out_set){
    
    out_set << "*Nset, nset=" << name << ", unsorted\n";
    out_set << node.number << endl;
}
    
void write_nodes_set(const std::string &name, const std::vector<Node> &set, std::ofstream &out_set){

    unsigned int compteur = 0;
    out_set << "*Nset, nset=" << name << ", unsorted\n";
    for (auto it = set.begin(); it != set.end()-1; ++it) {
        out_set << it->number << ", ";
        if(compteur == 15) {
            out_set << "\n";
            compteur = 0;
        }
        compteur++;
    }
    out_set << set.back().number << endl;
}

void append_perio_nodes(const cubic_mesh &cm_perio, const string &path_data, const string &outputfile) {

    
    std::string buffer;
    std::string path_inputfile = path_data + "/" + outputfile;
    std::ifstream in_set;
    in_set.open(path_inputfile, ios::in);
    bool empty_line = false;
    
    if(in_set) {
        while (!in_set.eof())
        {
            getline (in_set,buffer);
            if (buffer.empty())
                empty_line = true;
        }
    }
    else {
        cout << "Error: cannot open the file " << outputfile << " that details the nodes in the folder :" << path_data << endl;
        return;
    }
    in_set.close();
    
    std::ofstream out_set;
    
    out_set.open(path_inputfile, ios::app);
    if (empty_line == false) {
        out_set << "\n";
    }
    
    for (auto n : *cm_perio.Face_listXp)
        out_set << n;
    for (auto n : *cm_perio.Face_listYp)
        out_set << n;
    for (auto n : *cm_perio.Face_listZp)
        out_set << n;

    for (auto n : *cm_perio.Edge_listXpYm)
        out_set << n;
    for (auto n : *cm_perio.Edge_listXpYp)
        out_set << n;
    for (auto n : *cm_perio.Edge_listXmYp)
        out_set << n;

    for (auto n : *cm_perio.Edge_listXpZm)
        out_set << n;
    for (auto n : *cm_perio.Edge_listXpZp)
        out_set << n;
    for (auto n : *cm_perio.Edge_listXmZp)
        out_set << n;

    for (auto n : *cm_perio.Edge_listYpZm)
        out_set << n;
    for (auto n : *cm_perio.Edge_listYpZp)
        out_set << n;
    for (auto n : *cm_perio.Edge_listYmZp)
        out_set << n;

/*    out_set << *cm_perio.Corner_listXmYpZm;
    out_set << *cm_perio.Corner_listXpYmZm;
    out_set << *cm_perio.Corner_listXpYpZm;
    out_set << *cm_perio.Corner_listXmYmZp;
    out_set << *cm_perio.Corner_listXmYpZp;
    out_set << *cm_perio.Corner_listXpYmZp;
    out_set << *cm_perio.Corner_listXpYpZp;*/
}
    
void write_eq(ostream& s, const equation &eq) {
    
    s << "*Equation\n";
    s << eq.components.size() << "\n";
    auto it = eq.components.begin();
    unsigned int count = 0;
    for(it = eq.components.begin(); it != eq.components.end(); ++it) {
        if (it == eq.components.end() - 1) {
            s << it->node.number << ", " << it->dof << ", " <<  std::fixed << std::setprecision(6) << it->coef;
            s << "\n";
        }
        else {
            s << it->node.number << ", " << it->dof << ", " <<  std::fixed << std::setprecision(6) << it->coef << ", ";
            if (count == 2) {
                s << "\n";
                count = 0;
            }
            count++;
        }
    }
}
    
void write_PBC(const cubic_mesh &cm, const unsigned int &nb_nodes, const string &path_data, const string &outputfile){
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream out_set;
    
    out_set.open(filename, ios::out);
    
    out_set << "************************************\n";
    out_set << "*** SETS FOR BOUNDARY CONDITIONS ***\n";
    out_set << "************************************\n";

    for (unsigned int i=0; i<cm.list_of_corners.size(); i++) {
        if (cm.list_of_corners[i]->number != 0) {
            write_node_set(cm.set_name_corners[i], *cm.list_of_corners[i], out_set);
        }
    }
    for (unsigned int i=0; i<cm.list_of_edges.size(); i++) {
        if (cm.list_of_edges[i]->size() != 0) {
            write_nodes_set(cm.set_name_edges[i], *cm.list_of_edges[i], out_set);
        }
    }
    for (unsigned int i=0; i<cm.list_of_faces.size(); i++) {
        if (cm.list_of_faces[i]->size() != 0) {
            write_nodes_set(cm.set_name_faces[i], *cm.list_of_faces[i], out_set);
        }
    }

    out_set << "************************************\n";
    out_set << "** ALL NODES SET FOR OUTPUT RESULT *\n";
    out_set << "************************************\n";
    out_set << "*Nset, nset=Allnodes, generate\n";
    out_set << "1, " << nb_nodes << ", 1";
    out_set.close();
}

void write_TIE(const cubic_mesh &cm, const cubic_mesh &cm_perio, const string &path_data, const string &outputfile){
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream out_set;
    
    out_set.open(filename, ios::out);
    
    out_set << "************************************\n";
    out_set << "*** SETS FOR TIE CONDITIONS ***\n";
    out_set << "************************************" << endl;
    
    std::vector<Node> set;
    
    //Face Xp
    set.insert(std::end(set), std::begin(*cm_perio.Face_listXp), std::end(*cm_perio.Face_listXp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYm), std::end(*cm_perio.Edge_listXpYm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYp), std::end(*cm_perio.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZm), std::end(*cm_perio.Edge_listXpZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZp), std::end(*cm_perio.Edge_listXpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
    write_nodes_set("set_surface_perio_Xp", set, out_set);
    
    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listYp), std::end(*cm_perio.Face_listYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmYp), std::end(*cm_perio.Edge_listXmYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpYp), std::end(*cm_perio.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZm), std::end(*cm_perio.Edge_listYpZm));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZp), std::end(*cm_perio.Edge_listYpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZm);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
    write_nodes_set("set_surface_perio_Yp", set, out_set);

    set.clear();
    set.insert(std::end(set), std::begin(*cm_perio.Face_listZp), std::end(*cm_perio.Face_listZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXmZp), std::end(*cm_perio.Edge_listXmZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listXpZp), std::end(*cm_perio.Edge_listXpZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYmZp), std::end(*cm_perio.Edge_listYmZp));
    set.insert(std::end(set), std::begin(*cm_perio.Edge_listYpZp), std::end(*cm_perio.Edge_listYpZp));
    set.insert(std::end(set), *cm_perio.Corner_listXmYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXmYpZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYmZp);
    set.insert(std::end(set), *cm_perio.Corner_listXpYpZp);
    write_nodes_set("set_surface_perio_Zp", set, out_set);
    
    //Face Xp
    set.clear();
    set.insert(std::end(set), std::begin(*cm.Face_listXp), std::end(*cm.Face_listXp));
    set.insert(std::end(set), std::begin(*cm.Edge_listXpYm), std::end(*cm.Edge_listXpYm));
    set.insert(std::end(set), std::begin(*cm.Edge_listXpYp), std::end(*cm.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm.Edge_listXpZm), std::end(*cm.Edge_listXpZm));
//    set.insert(std::end(set), std::begin(*cm.Edge_listXpZp), std::end(*cm.Edge_listXpZp));
//    set.insert(std::end(set), *cm.Corner_listXpYmZm);
//    set.insert(std::end(set), *cm.Corner_listXpYpZm);
//    set.insert(std::end(set), *cm.Corner_listXpYmZp);
//    set.insert(std::end(set), *cm.Corner_listXpYpZp);
    write_nodes_set("set_surface_Xp", set, out_set);
    
    set.clear();
    set.insert(std::end(set), std::begin(*cm.Face_listYp), std::end(*cm.Face_listYp));
    set.insert(std::end(set), std::begin(*cm.Edge_listXmYp), std::end(*cm.Edge_listXmYp));
//    set.insert(std::end(set), std::begin(*cm.Edge_listXpYp), std::end(*cm.Edge_listXpYp));
    set.insert(std::end(set), std::begin(*cm.Edge_listYpZm), std::end(*cm.Edge_listYpZm));
    set.insert(std::end(set), std::begin(*cm.Edge_listYpZp), std::end(*cm.Edge_listYpZp));
//    set.insert(std::end(set), *cm.Corner_listXmYpZm);
//    set.insert(std::end(set), *cm.Corner_listXpYpZm);
//    set.insert(std::end(set), *cm.Corner_listXmYpZp);
//    set.insert(std::end(set), *cm.Corner_listXpYpZp);
    write_nodes_set("set_surface_Yp", set, out_set);
    
    set.clear();
    set.insert(std::end(set), std::begin(*cm.Face_listZp), std::end(*cm.Face_listZp));
    set.insert(std::end(set), std::begin(*cm.Edge_listXmZp), std::end(*cm.Edge_listXmZp));
    set.insert(std::end(set), std::begin(*cm.Edge_listXpZp), std::end(*cm.Edge_listXpZp));
    set.insert(std::end(set), std::begin(*cm.Edge_listYmZp), std::end(*cm.Edge_listYmZp));
//    set.insert(std::end(set), std::begin(*cm.Edge_listYpZp), std::end(*cm.Edge_listYpZp));
//    set.insert(std::end(set), *cm.Corner_listXmYmZp);
//    set.insert(std::end(set), *cm.Corner_listXmYpZp);
//    set.insert(std::end(set), *cm.Corner_listXpYmZp);
//    set.insert(std::end(set), *cm.Corner_listXpYpZp);
    write_nodes_set("set_surface_Zp", set, out_set);

    out_set << "*Surface, type=NODE, name=m_set_surface_Xp\n";
    out_set << "set_surface_perio_Xp, 1.\n";
    
    out_set << "*Surface, type=NODE, name=m_set_surface_Yp\n";
    out_set << "set_surface_perio_Yp, 1.\n";
    
    out_set << "*Surface, type=NODE, name=m_set_surface_Zp\n";
    out_set << "set_surface_perio_Zp, 1.\n";
    
    out_set << "*Surface, type=NODE, name=s_set_surface_Xp\n";
    out_set << "set_surface_Xp, 1.\n";

    out_set << "*Surface, type=NODE, name=s_set_surface_Yp\n";
    out_set << "set_surface_Yp, 1.\n";

    out_set << "*Surface, type=NODE, name=s_set_surface_Zp\n";
    out_set << "set_surface_Zp, 1.\n";
    
    out_set << "** Constraint: Constraint-Xp\n";
    out_set << "*Tie, name=Constraint-Xp, adjust=yes\n";
    out_set << "s_set_surface_Xp, m_set_surface_Xp\n";

    out_set << "** Constraint: Constraint-Yp\n";
    out_set << "*Tie, name=Constraint-Yp, adjust=yes\n";
    out_set << "s_set_surface_Yp, m_set_surface_Yp\n";

    out_set << "** Constraint: Constraint-Zp\n";
    out_set << "*Tie, name=Constraint-Zp, adjust=yes\n";
    out_set << "s_set_surface_Zp, m_set_surface_Zp\n";
    
    out_set.close();
}

/*void write_CDN(const cubic_mesh &cm, const string &path_data, const string &outputfile){
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream out_set;
    
    out_set.open(filename, ios::out);
    
    out_set << "************************************\n";
    out_set << "************ CENTRE NODE ***********\n";
    out_set << "************************************\n";
    out_set << "**\n";
    out_set << "*NSet, NSet=CentreNode, Unsorted\n";
    out_set << cm.center_node.number << "\n";
    out_set << "**\n";
    out_set << "************************************\n";
    out_set << "****** CONSTRAIN DRIVER NODES ******\n";
    out_set << "************************************\n";
    out_set << "**\n";
    out_set << "*Node\n";
    out_set << "1010011, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD11\n";
    out_set << "1010011\n";
    out_set << "*Node\n";
    out_set << "1020022, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD22\n";
    out_set << "1020022\n";
    out_set << "*Node\n";
    out_set << "1030033, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD33\n";
    out_set << "1030033\n";
    out_set << "*Node\n";
    out_set << "1040012, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD12\n";
    out_set << "1040012\n";
    out_set << "*Node\n";
    out_set << "1050013, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD13\n";
    out_set << "1050013\n";
    out_set << "*Node\n";
    out_set << "1060023, 0, 0, 0\n";
    out_set << "*NSet, NSet=CD23\n";
    out_set << "1060023\n";
    out_set << "**\n";
    out_set << "*NSet, NSet=CD_nodes, Unsorted\n";
    out_set << "1010011, 1020022, 1030033, 1040012, 1050013, 1060023\n";
    out_set << "**\n";
    out_set << "************************************\n";
    out_set << "********** MPC EQUATIONS ***********\n";
    out_set << "************************************\n";
    out_set << "**\n";
    out_set << "************** FACES ***************\n";
    out_set << "**\n";
    if (cm.Face_listXm->size() != 0) {
        out_set << "** equations between " << cm.set_name_faces[1] << " & " << cm.set_name_faces[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[1] << ", 1, 1.0, " << cm.set_name_faces[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[1] << ", 2, 1.0, " << cm.set_name_faces[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[1] << ", 3, 1.0, " << cm.set_name_faces[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "**\n";
    }
    if (cm.Face_listYm->size() != 0) {
        out_set << "** equations between " << cm.set_name_faces[3] << " & " << cm.set_name_faces[2] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[3] << ", 1, 1.0, " << cm.set_name_faces[2] << ", 1, -1.0, CD12, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[3] << ", 2, 1.0, " << cm.set_name_faces[2] << ", 2, -1.0, CD22, 1, " << -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[3] << ", 3, 1.0, " << cm.set_name_faces[2] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
    }
    if (cm.Face_listZm->size() != 0) {
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[5] << ", 1, 1.0, " << cm.set_name_faces[4] << ", 1, -1.0, CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[5] << ", 2, 1.0, " << cm.set_name_faces[4] << ", 2, -1.0, CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_faces[5] << ", 3, 1.0, " << cm.set_name_faces[4] << ", 3, -1.0, CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    out_set << "************** EDGES ***************\n";
    out_set << "**\n";
    if (cm.Edge_listXmYm->size() != 0) {
        out_set << "** equations between " << cm.set_name_edges[1] << " & " << cm.set_name_edges[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[1] << ", 1, 1.0, " << cm.set_name_edges[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[1] << ", 2, 1.0, " << cm.set_name_edges[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[1] << ", 3, 1.0, " << cm.set_name_edges[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "**\n";
        
        out_set << "** equations between " << cm.set_name_edges[2] << " & " << cm.set_name_edges[0] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[2] << ", 1, 1.0, " << cm.set_name_edges[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << ", CD12, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[2] << ", 2, 1.0, " << cm.set_name_edges[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << ", CD22, 1, " << -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[2] << ", 3, 1.0, " << cm.set_name_edges[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << ", CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
        
        out_set << "** equations between " << cm.set_name_edges[3] << " & " << cm.set_name_edges[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[3] << ", 1, 1.0, " << cm.set_name_edges[0] << ", 1, -1.0, CD12, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[3] << ", 2, 1.0, " << cm.set_name_edges[0] << ", 2, -1.0, CD22, 1, " << -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[3] << ", 3, 1.0, " << cm.set_name_edges[0] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
    }
    if (cm.Edge_listXmZm->size() != 0) {
        out_set << "** equations between " << cm.set_name_edges[5] << " & " << cm.set_name_edges[4] << "\n";        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[5] << ", 1, 1.0, " << cm.set_name_edges[4] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[5] << ", 2, 1.0, " << cm.set_name_edges[4] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[5] << ", 3, 1.0, " << cm.set_name_edges[4] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "**\n";
        out_set << "** equations between " << cm.set_name_edges[6] << " & " << cm.set_name_edges[4] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[6] << ", 1, 1.0, " << cm.set_name_edges[4] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << ", CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[6] << ", 2, 1.0, " << cm.set_name_edges[4] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << ", CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[6] << ", 3, 1.0, " << cm.set_name_edges[4] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << ", CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
        out_set << "** equations between " << cm.set_name_edges[7] << " & " << cm.set_name_edges[4] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[7] << ", 1, 1.0, " << cm.set_name_edges[4] << ", 1, -1.0, CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[7] << ", 2, 1.0, " << cm.set_name_edges[4] << ", 2, -1.0, CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[7] << ", 3, 1.0, " << cm.set_name_edges[4] << ", 3, -1.0, CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    if (cm.Edge_listYmZm->size() != 0) {
        out_set << "** equations between " << cm.set_name_edges[9] << " & " << cm.set_name_edges[8] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[9] << ", 1, 1.0, " << cm.set_name_edges[8] << ", 1, -1.0, CD12, 1, "<< -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[9] << ", 2, 1.0, " << cm.set_name_edges[8] << ", 2, -1.0, CD22, 1, "<< -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[9] << ", 3, 1.0, " << cm.set_name_edges[8] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
        out_set << "** equations between " << cm.set_name_edges[10] << " & " << cm.set_name_edges[8] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[10] << ", 1, 1.0, " << cm.set_name_edges[8] << ", 1, -1.0, CD12, 1, " << -1.*cm.Dy/2. << ", CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[10] << ", 2, 1.0, " << cm.set_name_edges[8] << ", 2, -1.0, CD22, 1, " << -1.*cm.Dy << ", CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_edges[10] << ", 3, 1.0, " << cm.set_name_edges[8] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << ", CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
        out_set << "** equations between " << cm.set_name_edges[11] << " & " << cm.set_name_edges[8] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[11] << ", 1, 1.0, " << cm.set_name_edges[8] << ", 1, -1.0, CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[11] << ", 2, 1.0, " << cm.set_name_edges[8] << ", 2, -1.0, CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_edges[11] << ", 3, 1.0, " << cm.set_name_edges[8] << ", 3, -1.0, CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    out_set << "************* CORNERS **************\n";
    out_set << "**\n";
    if((cm.Corner_listXpYmZm->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[1] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[1] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[1] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[1] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXmYpZm->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[3] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[3] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD12, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[3] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD22, 1, " << -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[3] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXmYmZp->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[4] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[4] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[4] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "3\n";
        out_set << cm.set_name_corners[4] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXpYpZm->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[2] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[2] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << ", CD12, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[2] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << ", CD22, 1, " << -1.*cm.Dy << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[2] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << ", CD23, 1, " << -1.*cm.Dy/2. << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXpYmZp->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[5] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[5] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << ", CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[5] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << ", CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[5] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << ", CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXmYpZp->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[7] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[7] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD12, 1, " << -1.*cm.Dy/2. << ", CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[7] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD22, 1, " << -1.*cm.Dy << ", CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "4\n";
        out_set << cm.set_name_corners[7] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD23, 1, " << -1.*cm.Dy/2. << ", CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    if((cm.Corner_listXpYpZp->number != 0)&&(cm.Corner_listXmYmZm->number != 0)) {
        out_set << "** equations between " << cm.set_name_corners[6] << " & " << cm.set_name_corners[0] << "\n";
        out_set << "*Equation\n";
        out_set << "5\n";
        out_set << cm.set_name_corners[6] << ", 1, 1.0, " << cm.set_name_corners[0] << ", 1, -1.0, CD11, 1, " << -1.*cm.Dx << ", CD12, 1, " << -1.*cm.Dy/2. << ",\n";
        out_set << "CD13, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "5\n";
        out_set << cm.set_name_corners[6] << ", 2, 1.0, " << cm.set_name_corners[0] << ", 2, -1.0, CD12, 1, " << -1.*cm.Dx/2. << ", CD22, 1, " << -1.*cm.Dy << ",\n";
        out_set << "CD23, 1, " << -1.*cm.Dz/2. << "\n";
        out_set << "*Equation\n";
        out_set << "5\n";
        out_set << cm.set_name_corners[6] << ", 3, 1.0, " << cm.set_name_corners[0] << ", 3, -1.0, CD13, 1, " << -1.*cm.Dx/2. << ", CD23, 1, " << -1.*cm.Dy/2. << ",\n";
        out_set << "CD33, 1, " << -1.*cm.Dz << "\n";
        out_set << "**\n";
    }
    out_set << "**\n";
    out_set.close();
}*/
    
void write_NonPerio_CDN(const cubic_mesh &cm, const cubic_mesh &cm_perio, const unsigned int &loading_type, const string &path_data, const string &outputfile){

    //Find the neighbours
    std::vector<std::vector<Node> > neigh;
    std::vector<mat> neigh_dist;
    unsigned int N=6;
    std::vector<mat> weights;
    std::string name_set;
    std::vector<Node> set;
    std::vector<Node> neigh_temp;
    mat weight_temp;
    std::vector<std::string> name_in(N);
    double min_dis = 1.E-6;
    
    //Construct all the equations
    cubic_equation cubic_eq(cm, cm_perio, loading_type);
    equation eq;
    
    std::vector<int> list_dofs;
    //Mechanical
    if (loading_type == 1) {
        list_dofs = {1,2,3};
    }
    //Thermomechanical
    else if (loading_type == 2) {
        list_dofs = {1,2,3};
    }
    //Thermal
    else if (loading_type == 3) {
        list_dofs = {11};
    }
    else
        cout << "error in Continuum_mechanics/Unit_cell/write.cpp : loading type should take values in the range (1,3)" << endl;
    
    std::string filename = path_data + "/" + outputfile;
    std::ofstream out_set;
    out_set.open(filename, ios::out);
    
    out_set << "************************************\n";
    out_set << "************ CENTRE NODE ***********\n";
    out_set << "************************************\n";
    out_set << "**\n";
    out_set << "*NSet, NSet=CentreNode, Unsorted\n";
    out_set << cm.center_node.number << "\n";
    out_set << "**\n";
    out_set << "************************************\n";
    out_set << "****** CONSTRAIN DRIVER NODES ******\n";
    out_set << "************************************\n";
    out_set << "**\n";
    
    for(unsigned int i=0; i<cubic_eq.CD.size();i++) {
        out_set << "*Node\n";
        out_set << cubic_eq.CD_nodes[i];
        write_node_set(cubic_eq.CD_set_name[i],  cubic_eq.CD_nodes[i], out_set);
    }
    write_nodes_set("CD_nodes", cubic_eq.CD_nodes, out_set);
    
    out_set << "**\n";
    out_set << "************************************\n";
    out_set << "********** MPC EQUATIONS ***********\n";
    out_set << "************************************\n";
    
    //Face Xp
    set_faceXp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Face_listXp, N);
    
    weights.resize(cm.Face_listXp->size());
    for (unsigned int i=0; i<cm.Face_listXp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Face_listXp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Face_listXp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i]; //extract a matrix
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 1, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_faces[1], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //Face Yp
    set_faceYp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Face_listYp, N);
    weights.resize(cm.Face_listYp->size());
    for (unsigned int i=0; i<cm.Face_listYp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Face_listYp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Face_listYp->size(), min_dis);
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 2, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_faces[3], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //Face Zp
    set_faceZp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Face_listZp, N);
    weights.resize(cm.Face_listZp->size());
    for (unsigned int i=0; i<cm.Face_listZp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Face_listZp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Face_listZp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 3, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_faces[5], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    /////////Edges////
    N=2;
    //EdgeXpYm
    set_EdgeXpYm(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXpYm, N);
    weights.resize(cm.Edge_listXpYm->size());
    for (unsigned int i=0; i<cm.Edge_listXpYm->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXpYm->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXpYm->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[1], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeXpYp
    set_EdgeXpYp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXpYp, N);
    weights.resize(cm.Edge_listXpYp->size());
    for (unsigned int i=0; i<cm.Edge_listXpYp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXpYp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXpYp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[2], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeXmYp
    set_EdgeXmYp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXmYp, N);
    weights.resize(cm.Edge_listXmYp->size());
    for (unsigned int i=0; i<cm.Edge_listXmYp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXmYp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXmYp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[3], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeXpZm
    set_EdgeXpZm(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXpZm, N);
    weights.resize(cm.Edge_listXpZm->size());
    for (unsigned int i=0; i<cm.Edge_listXpZm->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXpZm->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXpZm->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[5], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeXpZp
    set_EdgeXpZp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXpZp, N);
    weights.resize(cm.Edge_listXpZp->size());
    for (unsigned int i=0; i<cm.Edge_listXpZp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXpZp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXpZp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[6], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeXmZp
    set_EdgeXmZp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listXmZp, N);
    weights.resize(cm.Edge_listXmZp->size());
    for (unsigned int i=0; i<cm.Edge_listXmZp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listXmZp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listXmZp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[7], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeYpZm
    set_EdgeYpZm(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listYpZm, N);
    weights.resize(cm.Edge_listYpZm->size());
    for (unsigned int i=0; i<cm.Edge_listYpZm->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listYpZm->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listYpZm->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[9], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeYpZp
    set_EdgeYpZp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listYpZp, N);
    weights.resize(cm.Edge_listYpZp->size());
    for (unsigned int i=0; i<cm.Edge_listYpZp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listYpZp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listYpZp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[10], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //EdgeYmZp
    set_EdgeYmZp(set, cm_perio);
    //Find the neighbours
    find_neighbours_set(neigh, neigh_dist, set, *cm.Edge_listYmZp, N);
    weights.resize(cm.Edge_listYmZp->size());
    for (unsigned int i=0; i<cm.Edge_listYmZp->size(); i++) {
        weights[i] = zeros(neigh[i].size()-1,4);
    }
    //set_weights(weights, neigh_dist, cm.Edge_listYmZp->size(), min_dis);
    set_weights(weights, neigh_dist, cm.Edge_listYmZp->size(), min_dis);
    
    for (unsigned int i=0; i<neigh.size(); i++) {
        
        neigh_temp = neigh[i]; // extract a vector
        weight_temp = weights[i];
        for (auto j:list_dofs) {
            eq = set_equation(neigh_temp, weight_temp, 0, loading_type, j);
            //Replace in equation the nodes with the equations from cubic_equation
            replace_perio_eq(eq, cubic_eq, weight_temp, cm, cm.set_name_edges[11], loading_type, j);
            write_eq(out_set, eq);
        }
    }
    
    //Corner_listXpYmZm
    if(cm.Corner_listXpYmZm->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXpYmZm[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXpYpZm
    if(cm.Corner_listXpYpZm->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXpYpZm[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXmYpZm
    if(cm.Corner_listXmYpZm->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXmYpZm[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXmYmZp
    if(cm.Corner_listXmYmZp->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXmYmZp[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXpYmZp
    if(cm.Corner_listXpYmZp->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXpYmZp[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXpYpZp
    if(cm.Corner_listXpYpZp->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXpYpZp[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
    //Corner_listXmYpZp
    if(cm.Corner_listXmYpZp->number > 0) {
        for (auto j:list_dofs) {
            eq = cubic_eq.Corner_listXmYpZp[index_from_dof(j, loading_type)];
            write_eq(out_set, eq);
        }
    }
}
    
} //namespace simcoon
