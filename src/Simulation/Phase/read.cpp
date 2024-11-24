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
///@brief To read from NphasesX.dat and NlayerX.dat
///@version 1.0

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <map>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/read.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/eshelby.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void get_phase_charateristics(phase_characteristics &rve, const string &path_data) {

    string inputfile; //file # that stores the microstructure properties
    
    std::map<string, int> list_umat;
    list_umat = {{"ELISO",1},{"ELIST",2},{"ELORT",3},{"MIHEN",100},{"MIMTN",101},{"MISCN",103},{"MIPLN",104}};
    
    int method = list_umat[rve.sptr_matprops->umat_name];
    
    //first we read the behavior of the phases & we construct the tensors if necessary
    switch (method) {
            
        case 100: case 101: case 103: {
            
            //Definition of the static vectors x,wx,y,wy
            ellipsoid_multi::mp = rve.sptr_matprops->props(2);
            ellipsoid_multi::np = rve.sptr_matprops->props(3);
            ellipsoid_multi::x.set_size(ellipsoid_multi::mp);
            ellipsoid_multi::wx.set_size(ellipsoid_multi::mp);
            ellipsoid_multi::y.set_size(ellipsoid_multi::np);
            ellipsoid_multi::wy.set_size(ellipsoid_multi::np);
            points(ellipsoid_multi::x, ellipsoid_multi::wx, ellipsoid_multi::y, ellipsoid_multi::wy,ellipsoid_multi::mp, ellipsoid_multi::np);
            
            inputfile = "Nellipsoids" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_ellipsoid(rve, path_data, inputfile);
            break;
        }
        case 104: {
            inputfile = "Nlayers" + to_string(int(rve.sptr_matprops->props(1))) + ".dat";
            read_layer(rve, path_data, inputfile);
            break;
        }
    }
    
    for (unsigned int i=0; i<rve.sub_phases.size(); i++) {
        get_phase_charateristics(rve.sub_phases[i], path_data);
    }
}

void read_phase(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    unsigned int nphases = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nphases++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the phase characteristics in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nphases--;
    //Assert that the file has been filled correctly
    assert(nphases == rve.sptr_matprops->props(0));
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    rve.sub_phases_construct(nphases,0,1);
    
    int nprops = 0;
    int nstatev = 0;
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nprops >> nstatev;
        
        r.sptr_matprops->resize(nprops);
        natural_basis nb;
        r.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        r.sptr_sv_local->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> buffer;
        }
    }
    paramphases.close();
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        
        paramphases >> r.sptr_matprops->number >> r.sptr_matprops->umat_name >> r.sptr_matprops->save >>  r.sptr_shape->concentration >> r.sptr_matprops->psi_mat >> r.sptr_matprops->theta_mat >> r.sptr_matprops->phi_mat >> buffer >> buffer;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> r.sptr_matprops->props(j);
        }
        
        r.sptr_matprops->psi_mat*=(sim_pi/180.);
        r.sptr_matprops->theta_mat*=(sim_pi/180.);
        r.sptr_matprops->phi_mat*=(sim_pi/180.);
    }
    
    paramphases.close();
}

void read_layer(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    unsigned int nphases = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    std::shared_ptr<layer> sptr_layer;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nphases++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the layer characteristics in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nphases--;
    //Assert that the file has been filled correctly
    assert(nphases == rve.sptr_matprops->props(0));
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    rve.sub_phases_construct(nphases,1,1);
    
    int nprops = 0;
    int nstatev = 0;
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nprops >> nstatev;
        
        r.sptr_matprops->resize(nprops);
        natural_basis nb;
        r.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        r.sptr_sv_local->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> buffer;
        }
    }
    paramphases.close();
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        
        sptr_layer = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        paramphases >> r.sptr_matprops->number >> r.sptr_matprops->umat_name >> r.sptr_matprops->save >>  r.sptr_shape->concentration >> r.sptr_matprops->psi_mat >> r.sptr_matprops->theta_mat >> r.sptr_matprops->phi_mat >> sptr_layer->psi_geom >> sptr_layer->theta_geom >> sptr_layer->phi_geom >> buffer >> buffer;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> r.sptr_matprops->props(j);
        }
        
        r.sptr_matprops->psi_mat*=(sim_pi/180.);
        r.sptr_matprops->theta_mat*=(sim_pi/180.);
        r.sptr_matprops->phi_mat*=(sim_pi/180.);
        
        sptr_layer->psi_geom*=(sim_pi/180.);
        sptr_layer->theta_geom*=(sim_pi/180.);
        sptr_layer->phi_geom*=(sim_pi/180.);
    }
    
    paramphases.close();
}
    
void read_ellipsoid(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    unsigned int nphases = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    std::shared_ptr<ellipsoid> sptr_ellipsoid;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nphases++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the elipsoid characteristics in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nphases--;
    //Assert that the file has been filled correctly
    assert(nphases == rve.sptr_matprops->props(0));
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    rve.sub_phases_construct(nphases,2,1);
    
    int nprops = 0;
    int nstatev = 0;
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nprops >> nstatev;
        
        r.sptr_matprops->resize(nprops);
        natural_basis nb;
        r.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        r.sptr_sv_local->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> buffer;
        }
    }
    paramphases.close();
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        
        sptr_ellipsoid = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        paramphases >> r.sptr_matprops->number >> sptr_ellipsoid->coatingof >> r.sptr_matprops->umat_name >> r.sptr_matprops->save >>  sptr_ellipsoid->concentration >> r.sptr_matprops->psi_mat >> r.sptr_matprops->theta_mat >> r.sptr_matprops->phi_mat >> sptr_ellipsoid->a1 >> sptr_ellipsoid->a2 >>sptr_ellipsoid->a3 >> sptr_ellipsoid->psi_geom >> sptr_ellipsoid->theta_geom >> sptr_ellipsoid->phi_geom >> buffer >> buffer;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> r.sptr_matprops->props(j);
        }
        
        r.sptr_matprops->psi_mat*=(sim_pi/180.);
        r.sptr_matprops->theta_mat*=(sim_pi/180.);
        r.sptr_matprops->phi_mat*=(sim_pi/180.);
        
        sptr_ellipsoid->psi_geom*=(sim_pi/180.);
        sptr_ellipsoid->theta_geom*=(sim_pi/180.);
        sptr_ellipsoid->phi_geom*=(sim_pi/180.);
    }
    paramphases.close();
    
    //Fill the coatedby parameter
    std::shared_ptr<ellipsoid> sptr_ellipsoid_c;
    for (unsigned int i=0; i<nphases; i++) {
    		sptr_ellipsoid = std::dynamic_pointer_cast<ellipsoid>(rve.sub_phases[i].sptr_shape);
    		if (sptr_ellipsoid->coatingof != 0) {
    			sptr_ellipsoid_c = std::dynamic_pointer_cast<ellipsoid>(rve.sub_phases[sptr_ellipsoid->coatingof].sptr_shape);
    			sptr_ellipsoid_c->coatedby = i;
    		}
    }
}

    
void read_cylinder(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    unsigned int nphases = 0;
    std::string buffer;
    std::string path_inputfile = path_data + "/" + inputfile;
    std::ifstream paramphases;
    std::shared_ptr<cylinder> sptr_cylinder;
    
    paramphases.open(path_inputfile, ios::in);
    if(paramphases) {
        while (!paramphases.eof())
        {
            getline (paramphases,buffer);
            if (buffer != "") {
                nphases++;
            }
        }
    }
    else {
        cout << "Error: cannot open the file " << inputfile << " that details the cylinder characteristics in the folder :" << path_data << endl;
        return;
    }
    paramphases.close();
    nphases--;
    //Assert that the file has been filled correctly
    assert(nphases == rve.sptr_matprops->props(0));
    
    //Generate the sub_phase vector and har-create the objects pointed buy the shared_ptrs
    rve.sub_phases_construct(nphases,3,1);
    
    int nprops = 0;
    int nstatev = 0;
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer  >> buffer >> buffer;
    
    for(auto r : rve.sub_phases) {
        paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> nprops >> nstatev;
        
        r.sptr_matprops->resize(nprops);
        natural_basis nb;
        r.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        r.sptr_sv_local->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), rve.sptr_sv_global->T, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> buffer;
        }
    }
    paramphases.close();
    
    paramphases.open(path_inputfile, ios::in);
    paramphases >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer >> buffer;
    
    std::shared_ptr<material_characteristics> sptr_matprops1;
    for(auto r : rve.sub_phases) {
        
        sptr_cylinder = std::dynamic_pointer_cast<cylinder>(r.sptr_shape);
        
        paramphases >> r.sptr_matprops->number >> sptr_cylinder->coatingof >> r.sptr_matprops->umat_name >> r.sptr_matprops->save >>  sptr_cylinder->concentration >> r.sptr_matprops->psi_mat >> r.sptr_matprops->theta_mat >> r.sptr_matprops->phi_mat >> sptr_cylinder->L >> sptr_cylinder->R >> sptr_cylinder->psi_geom >> sptr_cylinder->theta_geom >> sptr_cylinder->phi_geom >> buffer >> buffer;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases >> r.sptr_matprops->props(j);
        }
        
        r.sptr_matprops->psi_mat*=(sim_pi/180.);
        r.sptr_matprops->theta_mat*=(sim_pi/180.);
        r.sptr_matprops->phi_mat*=(sim_pi/180.);
        
        sptr_cylinder->psi_geom*=(sim_pi/180.);
        sptr_cylinder->theta_geom*=(sim_pi/180.);
        sptr_cylinder->phi_geom*=(sim_pi/180.);
    }
    paramphases.close();
    
    //Fill the coatedby parameter
    std::shared_ptr<cylinder> sptr_cylinder_c;
    for (unsigned int i=0; i<nphases; i++) {
    		sptr_cylinder = std::dynamic_pointer_cast<cylinder>(rve.sub_phases[i].sptr_shape);
    		if (sptr_cylinder->coatingof != 0) {
    			sptr_cylinder_c = std::dynamic_pointer_cast<cylinder>(rve.sub_phases[sptr_cylinder->coatingof].sptr_shape);
    			sptr_cylinder_c->coatedby = i;
    		}
    }
}


} //namespace simcoon
