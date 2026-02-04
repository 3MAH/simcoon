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

///@file phase_characteristics.cpp
///@brief Characteristics of a phase, which hereditates from:
///- material_characteristics
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <memory>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Geometry/geometry.hpp>
#include <simcoon/Simulation/Geometry/layer.hpp>
#include <simcoon/Simulation/Geometry/cylinder.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/hyperelastic.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/phase_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/layer_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/ellipsoid_multi.hpp>
#include <simcoon/Continuum_mechanics/Homogenization/cylinder_multi.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Simulation/Solver/output.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

//=====Private methods for phase_characteristics===================================

//=====Public methods for phase_characteristics============================================

/*!
  \brief default constructor
*/
//-------------------------------------------------------------
phase_characteristics::phase_characteristics()
//-------------------------------------------------------------
{
    shape_type = 0;
    sv_type = 0;
    sptr_matprops = std::make_shared<material_characteristics>();
    
    //Note : the construction of sptr_shape = std::make_shared.. and sptr_sv = std::make_shared.. is made in construct(int, int)
   
}

/*!
  \brief Constructor with parameters
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
phase_characteristics::phase_characteristics(const int &mshape_type, const int &msv_type, const std::shared_ptr<geometry> &msptr_shape, const std::shared_ptr<phase_multi> &msptr_multi, const std::shared_ptr<material_characteristics> &msptr_matprops, const std::shared_ptr<state_variables> &msptr_sv_global, const std::shared_ptr<state_variables> &msptr_sv_local, const std::shared_ptr<std::ofstream> &msptr_out_global, const std::shared_ptr<std::ofstream> &msptr_out_local, const std::string &msub_phases_file)
//-------------------------------------------------------------
{
    shape_type = mshape_type;
    sv_type = msv_type;
    sptr_shape = msptr_shape;
    sptr_multi = msptr_multi;
    sptr_matprops = msptr_matprops;
    sptr_sv_global = msptr_sv_global;
    sptr_sv_local = msptr_sv_local;
    sptr_out_global = msptr_out_global;
    sptr_out_local = msptr_out_local;
    sub_phases_file = msub_phases_file;
}

/*!
  \brief Copy constructor
  \param s phase_characteristics object to duplicate
*/
    
//------------------------------------------------------
phase_characteristics::phase_characteristics(const phase_characteristics& pc)
//------------------------------------------------------
{
    shape_type = pc.shape_type;
    sv_type = pc.sv_type;
    sptr_shape = pc.sptr_shape;
    sptr_multi = pc.sptr_multi;
    sptr_matprops = pc.sptr_matprops;
    sptr_sv_global = pc.sptr_sv_global;
    sptr_sv_local = pc.sptr_sv_local;
    sptr_out_global = pc.sptr_out_global;
    sptr_out_local = pc.sptr_out_local;
    
    sub_phases = pc.sub_phases;
    sub_phases_file = pc.sub_phases_file;
}
    

/*!
  \brief Destructor

  Deletes phase_characteristics, the shared ptr related object will be detroyed automatically if they is no pointer to it.
*/
//-------------------------------------
phase_characteristics::~phase_characteristics() {}
//-------------------------------------

/*!
  \brief Standard operator = for phase_characteristics
*/
  
    
//-------------------------------------------------------------
void phase_characteristics::construct(const int &mshape_type, const int &msv_type)
//-------------------------------------------------------------
{
    assert(msv_type > 0);
    
    shape_type = mshape_type;
    sv_type = msv_type;
    
    //Switch case for the geometry of the phase
    switch (shape_type) {
        case 0: {
            sptr_shape = std::make_shared<geometry>();
            sptr_multi = std::make_shared<phase_multi>();
            break;
        }
        case 1: {
            sptr_shape = std::make_shared<layer>();
            sptr_multi = std::make_shared<layer_multi>();
            break;
        }
        case 2: {
            sptr_shape = std::make_shared<ellipsoid>();
            sptr_multi = std::make_shared<ellipsoid_multi>();
            break;
        }
        case 3: {
            sptr_shape = std::make_shared<cylinder>();
            sptr_multi = std::make_shared<cylinder_multi>();
            break;
        }
        default: {
            cout << "error: The geometry type does not correspond (0 for general, 1 for layer, 2 for ellipsoid, 3 for cylinder)\n";
            exit(0);
            break;
        }
    }
    
    //Switch case for the state_variables type of the phase
    switch (sv_type) {
        case 0: {
            sptr_sv_global = std::make_shared<state_variables>();
            sptr_sv_local = std::make_shared<state_variables>();
        }
        case 1: {
            sptr_sv_global = std::make_shared<state_variables_M>();
            sptr_sv_local = std::make_shared<state_variables_M>();
            break;
        }
        case 2: {
            sptr_sv_global = std::make_shared<state_variables_T>();
            sptr_sv_local = std::make_shared<state_variables_T>();
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }

}

//-------------------------------------------------------------
void phase_characteristics::sub_phases_construct(const int &nphases, const int &mshape_type, const int &msv_type)
//-------------------------------------------------------------
{
    assert(nphases > 0);
    assert(msv_type > 0);
    
    //First we resize the vector of sub_phases
    sub_phases.resize(nphases);

    //Second, for each component of the sub_phase vector we construct the phases
    for (int i=0; i<nphases; i++) {
        sub_phases[i].construct(mshape_type, msv_type);
    }
}
    
//----------------------------------------------------------------------
void phase_characteristics::to_start()
//----------------------------------------------------------------------
{
    
    switch (sv_type) {
        case 1: {
            auto sv_M_g = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
            auto sv_M_l = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_local);
            sv_M_g->to_start();
            sv_M_l->to_start();
            break;
        }
        case 2: {
            auto sv_T_g = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
            auto sv_T_l = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
            sv_T_g->to_start();
            sv_T_l->to_start();
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }
    sptr_multi->to_start();
    for(auto r : sub_phases) {
        r.to_start();
    }
    
}
    
//----------------------------------------------------------------------
void phase_characteristics::set_start(const int &corate_type)
//----------------------------------------------------------------------
{
    switch (sv_type) {
        case 1: {
            auto sv_M_g = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
            auto sv_M_l = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_local);
            sv_M_g->set_start(corate_type);
            sv_M_l->set_start(corate_type);
            break;
        }
        case 2: {
            auto sv_T_g = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
            auto sv_T_l = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
            sv_T_g->set_start(corate_type);
            sv_T_l->set_start(corate_type);
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }
    sptr_multi->set_start();
    for(auto r : sub_phases) {
        r.set_start(corate_type);
    }
}

//-------------------------------------------------------------
void phase_characteristics::local2global()
//-------------------------------------------------------------
{
    //Switch case for the state_variables type of the phase
    switch (sv_type) {
        case 1: {
            auto sv_M_g = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
            auto sv_M_l = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_local);
            sv_M_g->rotate_l2g(*sv_M_l, sptr_matprops->psi_mat, sptr_matprops->theta_mat, sptr_matprops->phi_mat);
            break;
        }
        case 2: {
            auto sv_T_g = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
            auto sv_T_l = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
            sv_T_g->rotate_l2g(*sv_T_l, sptr_matprops->psi_mat, sptr_matprops->theta_mat, sptr_matprops->phi_mat);
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }

}

//-------------------------------------------------------------
void phase_characteristics::global2local()
//-------------------------------------------------------------
{
    //Switch case for the state_variables type of the phase
    switch (sv_type) {
        case 1: {
            auto sv_M_g = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
            auto sv_M_l = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_local);
            sv_M_l->rotate_g2l(*sv_M_g, sptr_matprops->psi_mat, sptr_matprops->theta_mat, sptr_matprops->phi_mat);
            break;
        }
        case 2: {
            auto sv_T_g = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
            auto sv_T_l = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
            sv_T_l->rotate_g2l(*sv_T_g, sptr_matprops->psi_mat, sptr_matprops->theta_mat, sptr_matprops->phi_mat);
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }
}

//----------------------------------------------------------------------
phase_characteristics& phase_characteristics::operator = (const phase_characteristics& pc)
//----------------------------------------------------------------------
{
    shape_type = pc.shape_type;
    sv_type = pc.sv_type;
    sptr_shape = pc.sptr_shape;
    sptr_matprops = pc.sptr_matprops;
    sptr_sv_global = pc.sptr_sv_global;
    sptr_sv_local = pc.sptr_sv_local;
    sptr_out_global = pc.sptr_out_global;
    sptr_out_local = pc.sptr_out_local;
    
    sub_phases = pc.sub_phases;
    sub_phases_file = pc.sub_phases_file;
    
    return *this;
}

//----------------------------------------------------------------------
void phase_characteristics::define_output(const std::string &path, const std::string &outputfile, const std::string &coordsys)
//----------------------------------------------------------------------
{

    std::string ext_filename = outputfile.substr(outputfile.length()-4,outputfile.length());
    std::string filename = outputfile.substr(0,outputfile.length()-4); //to remove the extension
//    if(sptr_matprops->number > 0)
    filename = filename + '-' + std::to_string(sptr_matprops->number) + ext_filename;
    std::string path_filename = path + "/" + filename;
//    else
//        filename = filename + ext_filename;
    
//    std::ofstream of_file(filename);
    if(coordsys == "global") {
        sptr_out_global = make_shared<ofstream>(path_filename);
    }
    else if(coordsys == "local") {
        sptr_out_local = make_shared<ofstream>(path_filename);
    }
    
    for(unsigned int i=0; i<sub_phases.size(); i++) {
        sub_phases[i].define_output(path, filename, coordsys);
    }
    
}
    
//----------------------------------------------------------------------
void phase_characteristics::output(const solver_output &so, const int &kblock, const int &kcycle, const int&kstep, const int &kinc, const double & Time, const std::string &coordsys)
//----------------------------------------------------------------------
{

    if(coordsys == "global") {
        *sptr_out_global << kblock+1 << "\t";
        *sptr_out_global << kcycle+1 << "\t";
        *sptr_out_global << kstep+1 << "\t";
        *sptr_out_global << kinc+1 << "\t";
        *sptr_out_global << Time << "\t\t";
        
        //Switch case for the state_variables type of the phase
        if (so.o_nb_T) {
    
            switch (sv_type) {
                case 1: {
                    *sptr_out_global << sptr_sv_global->T  << "\t";
                    *sptr_out_global << 0 << "\t";                //This is for the flux Q
                    *sptr_out_global << 0 << "\t";                //This is for the rpl
                    break;
                }
                case 2: {
                    //We need to cast sv
                    std::shared_ptr<state_variables_T> sv_T = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
                    *sptr_out_global << sv_T->T  << "\t";
                    *sptr_out_global << sv_T->Q << "\t";                //This is for the flux
                    *sptr_out_global << sv_T->r << "\t";                //This is for the r
                    break;
                }
                default: {
                    cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
                    exit(0);
                    break;
                }
            }
        }
    
        //output
        if (so.o_nb_strain) {
            switch (so.o_strain_type) {
                case 0: {
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_global << sptr_sv_global->Etot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 1: {
                        vec E_biot;
                        try {
                            E_biot = t2v_strain(sqrtmat_sympd(2.*v2t_strain(sptr_sv_global->Etot)+eye(3,3)) - eye(3,3));
                        } catch (const std::runtime_error &e) {
                            cerr << "Error in sqrtmat_sympd : " << e.what() << endl;
                            throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside phase_characteristics::output.");
                        }                                    
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_global << E_biot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 2: {
                    vec F_vec = vectorise(sptr_sv_global->F1.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_global << F_vec(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 3: {
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_global << sptr_sv_global->etot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 4: {
                    mat R_temp = zeros(3,3);
                    mat V = zeros(3,3);                    
                    VR_decomposition(V, R_temp, sptr_sv_global->F1);
                    vec lambda_bar = isochoric_pstretch_from_V(V);
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_global << lambda_bar(so.o_strain(z)) << "\t";
                    }
                    break;
                }                
                default: {
                    cout << "Error in phase_characteristics::output : The output strain type is not valid (0 : Green-Lagrange, 1 for logarithmic) : " << so.o_strain_type << endl;
                    exit(0);
                }
                
            }
        }
        if (so.o_nb_stress) {
            switch (so.o_stress_type) {
                case 0: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_global << sptr_sv_global->PKII(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 1: {
                    mat PKI = Kirchoff2PKI(v2t_stress(sptr_sv_global->tau), sptr_sv_global->F1);
                    mat Nominal_stress = PKI.t();
                    vec Nominal_stress_vec = vectorise(Nominal_stress.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_global << Nominal_stress_vec(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 2: {
                    mat PKI = Kirchoff2PKI(v2t_stress(sptr_sv_global->tau), sptr_sv_global->F1);
                    vec PK1_vec = vectorise(PKI.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_global << PK1_vec(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 3: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_global << sptr_sv_global->tau(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 4: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_global << sptr_sv_global->sigma(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                default: {
                    cout << "Error in phase_characteristics::output : The output stres type is not valid (0 : Piola-Kirchoff II, 1 for Kirchoff, 2 for Cauchy) : " << so.o_stress_type << endl;
                    exit(0);
                }
            }
        }

        switch (so.o_rotation_type) {
            case 1: {
                vec R_vec = vectorise(sptr_sv_global->R.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_global << R_vec(z) << "\t";
                }
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        *sptr_out_global << sptr_sv_global->nb.g_i[i](j) << "\t";
                    }
                }
                break;
            }
            case 2: {
                vec DR_vec = vectorise(sptr_sv_global->DR.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_global << DR_vec(z) << "\t";
                }
                break;
            }
            case 3: {
                vec R_vec = vectorise(sptr_sv_global->R.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_global << R_vec(z) << "\t";
                }
                vec DR_vec = vectorise(sptr_sv_global->DR.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_global << DR_vec(z) << "\t";
                }
                break;
            }
            default: {
                break;
            }
        }

        switch (so.o_tangent_modulus) {
            case 1: {
                switch (sv_type) {
                    case 1: {
                        std::shared_ptr<state_variables_M> sv_M = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
                        vec Lt_vec = vectorise(sv_M->Lt.t()); //The transpose is to obtain a vec with row-wise concatenation
                        for (int z=0; z<36; z++)
                            *sptr_out_global << Lt_vec(z) << "\t";
                        break;
                    }
                    case 2: {
                        //We need to cast sv
                        std::shared_ptr<state_variables_T> sv_T = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
                        vec dSdE_vec = vectorise(sv_T->dSdE.t()); //The transpose is to obtain a vec with row-wise concatenation
                        for (int z=0; z<36; z++)
                            *sptr_out_global << dSdE_vec(z) << "\t";
                        for (int z=0; z<6; z++)
                            *sptr_out_global << sv_T->dSdT(z,0) << "\t";
                        for (int z=0; z<6; z++)
                            *sptr_out_global << sv_T->drdE(z,0) << "\t";
                        *sptr_out_global << sv_T->drdT(0,0) << "\t";
                        break;
                    }
                    default: {
                        cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
                        exit(0);
                        break;
                    }
                }
            }
            default: {
                break;
            }
        }
        
        switch (sv_type) {
            case 1: {
                std::shared_ptr<state_variables_M> sv_M = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_global);
                *sptr_out_global << sv_M->Wm(0)  << "\t";
                *sptr_out_global << sv_M->Wm(1)  << "\t";
                *sptr_out_global << sv_M->Wm(2)  << "\t";
                *sptr_out_global << sv_M->Wm(3)  << "\t";
                break;
            }
            case 2: {
                //We need to cast sv
                std::shared_ptr<state_variables_T> sv_T = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_global);
                *sptr_out_global << sv_T->Wm(0)  << "\t";
                *sptr_out_global << sv_T->Wm(1)  << "\t";
                *sptr_out_global << sv_T->Wm(2)  << "\t";
                *sptr_out_global << sv_T->Wm(3)  << "\t";
                *sptr_out_global << sv_T->Wt(0)  << "\t";
                *sptr_out_global << sv_T->Wt(1)  << "\t";
                *sptr_out_global << sv_T->Wt(2)  << "\t";
                break;
            }
            default: {
                cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
                exit(0);
                break;
            }
        }
        
        *sptr_out_global << "\t";
        if(so.o_nw_statev != 0){
            if (so.o_wanted_statev(0) < 0) {
                for(int k = 0 ; k < sptr_sv_global->nstatev ; k++)
                *sptr_out_global << sptr_sv_global->statev(k) << "\t";
            }
            else{
                for(int k = 0 ; k < so.o_nw_statev ; k++){
                    for (int l = so.o_wanted_statev(k); l < (so.o_range_statev(k)+1); l++){
                        *sptr_out_global << sptr_sv_global->statev(l) << "\t";
                    }
                }
            }
        }
        *sptr_out_global << endl;
        
        for(auto r : sub_phases) {
            r.output(so, kblock, kcycle, kstep, kinc, Time, "global");
        }
    }
    else if(coordsys == "local") {
    
        *sptr_out_local << kblock+1 << "\t";
        *sptr_out_local << kcycle+1 << "\t";
        *sptr_out_local << kstep+1 << "\t";
        *sptr_out_local << kinc+1 << "\t";
        *sptr_out_local << Time << "\t\t";
        
        //Switch case for the state_variables type of the phase
        if (so.o_nb_T) {
            
            switch (sv_type) {
                    case 1: {
                        *sptr_out_local << sptr_sv_local->T  << "\t";
                        *sptr_out_local << 0 << "\t";                //This is for the flux Q
                        *sptr_out_local << 0 << "\t";                //This is for the r
                        break;
                    }
                    case 2: {
                        //We need to cast sv
                        std::shared_ptr<state_variables_T> sv_T = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
                        *sptr_out_local << sv_T->T  << "\t";
                        *sptr_out_local << sv_T->Q << "\t";                //This is for the flux
                        *sptr_out_local << sv_T->r << "\t";                //This is for the r
                        break;
                    }
                default: {
                    cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
                    exit(0);
                    break;
                }
            }
        }
    
        //output
        if (so.o_nb_strain) {
            switch (so.o_strain_type) {
                case 0: {
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_local << sptr_sv_local->Etot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 1: {
                        vec E_biot;
                        try {
                            E_biot = t2v_strain(sqrtmat_sympd(2.*v2t_strain(sptr_sv_local->Etot)+eye(3,3)) - eye(3,3));
                        } catch (const std::runtime_error &e) {
                            cerr << "Error in sqrtmat_sympd : " << e.what() << endl;
                            throw simcoon::exception_sqrtmat_sympd("Error in sqrtmat_sympd function inside phase_characteristics::output.");
                        }                                    
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_local << E_biot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 2: {
                    vec F_vec = vectorise(sptr_sv_local->F1.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_local << F_vec(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 3: {
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_local << sptr_sv_local->etot(so.o_strain(z)) << "\t";
                    }
                    break;
                }
                case 4: {
                    mat R_temp = zeros(3,3);
                    mat V = zeros(3,3);                    
                    VR_decomposition(V, R_temp, sptr_sv_local->F1);
                    vec lambda_bar = isochoric_pstretch_from_V(V);
                    for (int z=0; z<so.o_nb_strain; z++) {
                        *sptr_out_local << lambda_bar(so.o_strain(z)) << "\t";
                    }
                    break;
                }                
                default: {
                    cout << "Error in phase_characteristics::output : The output strain type is not valid (0 : Green-Lagrange, 1 for logarithmic) : " << so.o_strain_type << endl;
                    exit(0);
                }
                
            }
        }
        if (so.o_nb_stress) {
            switch (so.o_stress_type) {
                case 0: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_local << sptr_sv_local->PKII(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 1: {
                    mat PKI = Kirchoff2PKI(v2t_stress(sptr_sv_local->tau), sptr_sv_local->F1);
                    mat Nominal_stress = PKI.t();
                    vec Nominal_stress_vec = vectorise(Nominal_stress.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_local << Nominal_stress_vec(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 2: {
                    mat PKI = Kirchoff2PKI(v2t_stress(sptr_sv_local->tau), sptr_sv_global->F1);
                    vec PK1_vec = vectorise(PKI.t()); //The transpose is to obtain a vec with row-wise concatenation
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_local << PK1_vec(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 3: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_local << sptr_sv_local->tau(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                case 4: {
                    for (int z=0; z<so.o_nb_stress; z++) {
                        *sptr_out_local << sptr_sv_local->sigma(so.o_stress(z)) << "\t";
                    }
                    break;
                }
                default: {
                    cout << "Error in phase_characteristics::output : The output stres type is not valid (0 : Piola-Kirchoff II, 1 for Kirchoff, 2 for Cauchy) : " << so.o_stress_type << endl;
                    exit(0);
                }
            }
        }
        
        switch (so.o_rotation_type) {
            case 1: {
                vec R_vec = vectorise(sptr_sv_local->R.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_local << R_vec(z) << "\t";
                }
                break;
            }
            case 2: {
                vec DR_vec = vectorise(sptr_sv_local->DR.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_local << DR_vec(z) << "\t";
                }
                break;
            }
            case 3: {
                vec R_vec = vectorise(sptr_sv_local->R.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_local << R_vec(z) << "\t";
                }
                vec DR_vec = vectorise(sptr_sv_local->DR.t()); //The transpose is to obtain a vec with row-wise concatenation
                for (int z=0; z<9; z++) {
                    *sptr_out_local << DR_vec(z) << "\t";
                }
                break;
            }
            default: {
                break;
            }
        }

        switch (sv_type) {
            case 1: {
                std::shared_ptr<state_variables_M> sv_M = std::dynamic_pointer_cast<state_variables_M>(sptr_sv_local);
                *sptr_out_local << sv_M->Wm(0)  << "\t";
                *sptr_out_local << sv_M->Wm(1)  << "\t";
                *sptr_out_local << sv_M->Wm(2)  << "\t";
                *sptr_out_local << sv_M->Wm(3)  << "\t";
                break;
            }
            case 2: {
                //We need to cast sv
                std::shared_ptr<state_variables_T> sv_T = std::dynamic_pointer_cast<state_variables_T>(sptr_sv_local);
                *sptr_out_local << sv_T->Wm(0)  << "\t";
                *sptr_out_local << sv_T->Wm(1)  << "\t";
                *sptr_out_local << sv_T->Wm(2)  << "\t";
                *sptr_out_local << sv_T->Wm(3)  << "\t";
                *sptr_out_local << sv_T->Wt(0)  << "\t";
                *sptr_out_local << sv_T->Wt(1)  << "\t";
                *sptr_out_local << sv_T->Wt(2)  << "\t";
                break;
            }
            default: {
                cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
                exit(0);
                break;
            }
        }
        
        *sptr_out_local << "\t";
        if(so.o_nw_statev != 0){
            if (so.o_wanted_statev(0) < 0) {
                for(int k = 0 ; k < sptr_sv_local->nstatev ; k++)
                *sptr_out_local << sptr_sv_local->statev(k) << "\t";
            }
            else{
                for(int k = 0 ; k < so.o_nw_statev ; k++){
                    for (int l = so.o_wanted_statev(k); l < (so.o_range_statev(k)+1); l++){
                        *sptr_out_local << sptr_sv_local->statev(l) << "\t";
                    }
                }
            }
        }
        *sptr_out_local << endl;
        
        for(auto r : sub_phases) {
            r.output(so, kblock, kcycle, kstep, kinc, Time, "local");
        }
        
    }
        

    
}
    
    
//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const phase_characteristics& pc)
//--------------------------------------------------------------------------
{
	s << "Display phase characteristics:\n";
	s << "Display geometry characteristics:\n";
    s << *pc.sptr_shape;
	s << "Display material characteristics:\n";
    s << *pc.sptr_matprops;
    s << "Display global state variables:\n";
    s << *pc.sptr_sv_global;
    s << "Display local state variables:\n";
    s << *pc.sptr_sv_local;
    
    for(auto r : pc.sub_phases) {
        s << r;
    }
    
    s << "\n\n";

	return s;
}
    
//----------------------------------------------------------------------
void phase_characteristics::copy(const phase_characteristics& pc)
//----------------------------------------------------------------------
{
    shape_type =  pc.shape_type;
    sv_type = pc.sv_type;
    sptr_matprops = std::make_shared<material_characteristics>(*pc.sptr_matprops);
    
    //Switch case for the geometry of the phase
    switch (shape_type) {
        case 0: {
            sptr_shape = std::make_shared<geometry>(*pc.sptr_shape);
            sptr_multi = std::make_shared<phase_multi>(*pc.sptr_multi);
            break;
        }
        case 1: {
            std::shared_ptr<layer> lay = std::dynamic_pointer_cast<layer>(pc.sptr_shape);
            std::shared_ptr<layer_multi> lay_multi = std::dynamic_pointer_cast<layer_multi>(pc.sptr_multi);
            sptr_shape = std::make_shared<layer>(*lay);
            sptr_multi = std::make_shared<layer_multi>(*lay_multi);
            break;
        }
        case 2: {
            std::shared_ptr<ellipsoid> elli = std::dynamic_pointer_cast<ellipsoid>(pc.sptr_shape);
            std::shared_ptr<ellipsoid_multi> elli_multi = std::dynamic_pointer_cast<ellipsoid_multi>(pc.sptr_multi);
            sptr_shape = std::make_shared<ellipsoid>(*elli);
            sptr_multi = std::make_shared<ellipsoid_multi>(*elli_multi);
            break;
        }
        case 3: {
            std::shared_ptr<cylinder> cyl = std::dynamic_pointer_cast<cylinder>(pc.sptr_shape);
            std::shared_ptr<cylinder_multi> cyl_multi = std::dynamic_pointer_cast<cylinder_multi>(pc.sptr_multi);
            sptr_shape = std::make_shared<cylinder>(*cyl);
            sptr_multi = std::make_shared<cylinder_multi>(*cyl_multi);
            break;
        }
        default: {
            cout << "error: The geometry type does not correspond (0 for general, 1 for layer, 2 for ellipsoid, 3 for cylinder)\n";
            exit(0);
            break;
        }
    }
    
    //Switch case for the state_variables type of the phase
    switch (sv_type) {
        case 0: {
            sptr_sv_global = std::make_shared<state_variables>(*pc.sptr_sv_global);
            sptr_sv_local = std::make_shared<state_variables>(*pc.sptr_sv_local);
            break;
        }
        case 1: {
            std::shared_ptr<state_variables_M> sv_M_g = std::dynamic_pointer_cast<state_variables_M>(pc.sptr_sv_global);
            std::shared_ptr<state_variables_M> sv_M_l = std::dynamic_pointer_cast<state_variables_M>(pc.sptr_sv_local);
            sptr_sv_global = std::make_shared<state_variables_M>(*sv_M_g);
            sptr_sv_local = std::make_shared<state_variables_M>(*sv_M_l);
            break;
        }
        case 2: {
            std::shared_ptr<state_variables_T> sv_T_g = std::dynamic_pointer_cast<state_variables_T>(pc.sptr_sv_global);
            std::shared_ptr<state_variables_T> sv_T_l = std::dynamic_pointer_cast<state_variables_T>(pc.sptr_sv_local);
            sptr_sv_global = std::make_shared<state_variables_T>(*sv_T_g);
            sptr_sv_local = std::make_shared<state_variables_T>(*sv_T_l);
            break;
        }
        default: {
            cout << "error: The state_variable type does not correspond (1 for Mechanical, 2 for Thermomechanical)\n";
            exit(0);
            break;
        }
    }

    sub_phases.clear();
    for (unsigned int i=0; i<pc.sub_phases.size(); i++) {
        phase_characteristics temp;
        temp.copy(pc.sub_phases[i]);
        sub_phases.push_back(temp);
    }
    sub_phases_file = pc.sub_phases_file;
}

} //namespace simcoon
