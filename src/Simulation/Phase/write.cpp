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

#include <assert.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/ellipsoid.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

void write_phase(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    std::string buffer;
    std::string filename = path_data + "/" + inputfile;
    std::ofstream paramphases;
    
    paramphases.open(filename, ios::out);
    
    paramphases << "Number\t" << "umat\t" << "save\t" << "c\t" << "psi_mat\t" << "theta_mat\t" << "phi_mat\t" << "nprops\t" << "nstatev\t" << "props\n";
    
    for(auto r : rve.sub_phases) {
        
        r.sptr_matprops->psi_mat*=(180./simcoon::pi);
        r.sptr_matprops->theta_mat*=(180./simcoon::pi);
        r.sptr_matprops->phi_mat*=(180./simcoon::pi);
        
        paramphases << r.sptr_matprops->number << "\t" << r.sptr_matprops->umat_name << "\t" << r.sptr_matprops->save << "\t" <<  r.sptr_shape->concentration << "\t" << r.sptr_matprops->psi_mat << "\t" << r.sptr_matprops->theta_mat << "\t" << r.sptr_matprops->phi_mat << "\t" << r.sptr_matprops->nprops << "\t" << r.sptr_sv_global->nstatev;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases << "\t" << r.sptr_matprops->props(j);
        }
        paramphases << "\n";
    }
    paramphases.close();
}
    
void write_layer(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    std::string buffer;
    std::string filename = path_data + "/" + inputfile;
    std::ofstream paramphases;
    
    paramphases.open(filename, ios::out);
    
    paramphases << "Number\t" << "umat\t" << "save\t" << "c\t" << "psi_mat\t" << "theta_mat\t" << "phi_mat\t" << "psi_geom\t" << "theta_geom\t"	<< "phi_geom\t" << "nprops\t" << "nstatev\t" << "props\n";
    
    for(auto r : rve.sub_phases) {

        r.sptr_matprops->psi_mat*=(180./simcoon::pi);
        r.sptr_matprops->theta_mat*=(180./simcoon::pi);
        r.sptr_matprops->phi_mat*=(180./simcoon::pi);

        auto sptr_layer = std::dynamic_pointer_cast<layer>(r.sptr_shape);
        sptr_layer->psi_geom*=(180./simcoon::pi);
        sptr_layer->theta_geom*=(180./simcoon::pi);
        sptr_layer->phi_geom*=(180./simcoon::pi);
        
        paramphases << r.sptr_matprops->number << "\t" << r.sptr_matprops->umat_name << "\t" << r.sptr_matprops->save << "\t" <<  sptr_layer->concentration << "\t" << r.sptr_matprops->psi_mat << "\t" << r.sptr_matprops->theta_mat << "\t" << r.sptr_matprops->phi_mat << "\t" << sptr_layer->psi_geom << "\t" << sptr_layer->theta_geom << "\t" << sptr_layer->phi_geom << "\t" << r.sptr_matprops->nprops << "\t" << r.sptr_sv_global->nstatev;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases << "\t" << r.sptr_matprops->props(j);
        }
        paramphases << "\n";
    }
    paramphases.close();
}

void write_ellipsoid(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    std::string buffer;
    std::string filename = path_data + "/" + inputfile;
    std::ofstream paramphases;
    
    paramphases.open(filename, ios::out);
    
    paramphases << "Number\t" << "Coatingof\t" << "umat\t" << "save\t" << "c\t" << "psi_mat\t" << "theta_mat\t" << "phi_mat\t" << "a1\t" << "a2\t" << "a3\t" << "psi_geom\t" << "theta_geom\t"	<< "phi_geom\t" << "nprops\t" << "nstatev\t" << "props\n";
    
    for(auto r : rve.sub_phases) {
        
        r.sptr_matprops->psi_mat*=(180./simcoon::pi);
        r.sptr_matprops->theta_mat*=(180./simcoon::pi);
        r.sptr_matprops->phi_mat*=(180./simcoon::pi);
        
        auto sptr_ellipsoid = std::dynamic_pointer_cast<ellipsoid>(r.sptr_shape);
        sptr_ellipsoid->psi_geom*=(180./simcoon::pi);
        sptr_ellipsoid->theta_geom*=(180./simcoon::pi);
        sptr_ellipsoid->phi_geom*=(180./simcoon::pi);
        
        paramphases << r.sptr_matprops->number << "\t" << sptr_ellipsoid->coatingof << "\t" << r.sptr_matprops->umat_name << "\t" << r.sptr_matprops->save << "\t" << sptr_ellipsoid->concentration << "\t" << r.sptr_matprops->psi_mat << "\t" << r.sptr_matprops->theta_mat << "\t" << r.sptr_matprops->phi_mat << "\t" << sptr_ellipsoid->a1 << "\t" << sptr_ellipsoid->a2 << "\t" << sptr_ellipsoid->a3 << "\t" << sptr_ellipsoid->psi_geom << "\t" << sptr_ellipsoid->theta_geom << "\t" << sptr_ellipsoid->phi_geom << "\t" << r.sptr_matprops->nprops << "\t" << r.sptr_sv_global->nstatev;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases << "\t" << r.sptr_matprops->props(j);
        }
        paramphases << "\n";
    }
    paramphases.close();
}

void write_cylinder(phase_characteristics &rve, const string &path_data, const string &inputfile) {
    
    std::string buffer;
    std::string filename = path_data + "/" + inputfile;
    std::ofstream paramphases;
    
    paramphases.open(filename, ios::out);
    
    paramphases << "Number\t" << "Coatingof\t" << "umat\t" << "save\t" << "c\t" << "psi_mat\t" << "theta_mat\t" << "phi_mat\t" << "L\t" << "R\t" << "psi_geom\t" << "theta_geom\t"	<< "phi_geom\t" << "nprops\t" << "nstatev\t" << "props\n";
    
    for(auto r : rve.sub_phases) {
        
        r.sptr_matprops->psi_mat*=(180./simcoon::pi);
        r.sptr_matprops->theta_mat*=(180./simcoon::pi);
        r.sptr_matprops->phi_mat*=(180./simcoon::pi);
        
        auto sptr_cylinder = std::dynamic_pointer_cast<cylinder>(r.sptr_shape);
        sptr_cylinder->psi_geom*=(180./simcoon::pi);
        sptr_cylinder->theta_geom*=(180./simcoon::pi);
        sptr_cylinder->phi_geom*=(180./simcoon::pi);
        
        
        paramphases << r.sptr_matprops->number << "\t" << sptr_cylinder->coatingof << "\t" << r.sptr_matprops->umat_name << "\t" << r.sptr_matprops->save << "\t" <<  sptr_cylinder->concentration << "\t" << r.sptr_matprops->psi_mat << "\t" << r.sptr_matprops->theta_mat << "\t" << r.sptr_matprops->phi_mat << "\t" << sptr_cylinder->L << "\t" << sptr_cylinder->R << "\t" << sptr_cylinder->psi_geom << "\t" << sptr_cylinder->theta_geom << "\t" << sptr_cylinder->phi_geom << "\t" << r.sptr_matprops->nprops << "\t" << r.sptr_sv_global->nstatev;
        
        for(int j=0; j<r.sptr_matprops->nprops; j++) {
            paramphases << "\t" << r.sptr_matprops->props(j);
        }
        paramphases << "\n";
    }
    paramphases.close();
}

} //namespace simcoon
