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

///@file solver.cpp
///@brief solver: solve the mechanical thermomechanical equilibrium			//
//	for a homogeneous loading path, allowing repeatable steps
///@version 1.9

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include <memory>
#include <filesystem>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Simulation/Phase/material_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/stress.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Functions/natural_basis.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Simulation/Solver/solver_sink.hpp>

using namespace std;
using namespace arma;

namespace {

// Single definition of the step-cut decision: bisect the increment unless it
// is already at its minimal fraction. Returns whether the cut was applied —
// call sites decide what a refusal means (throw a typed error, rethrow, ...).
inline bool try_step_cut(const double &Dtinc_cur, const double &Dn_mini,
                         const double &div_tnew_dt, double &tnew_dt) {
    if (fabs(Dtinc_cur - Dn_mini) > simcoon::iota) {
        tnew_dt = div_tnew_dt;
        return true;
    }
    return false;
}

// Step-cut-or-rethrow policy shared by every recoverable-failure catch in the
// solver Newton loops. Must only be called from inside a catch block.
inline void step_cut_or_rethrow(const double &Dtinc_cur, const double &Dn_mini,
                                const double &div_tnew_dt, double &tnew_dt,
                                int &compteur) {
    if (try_step_cut(Dtinc_cur, Dn_mini, div_tnew_dt, tnew_dt)) {
        compteur = 0;
    } else {
        throw;
    }
}

}  // namespace

namespace simcoon{

int solver_run(std::vector<block> &blocks, const double &T_init, const solver_output &so, const string &umat_name, const vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const solver_params &ctrl, solver_results_sink &sink) {

    if (ctrl.tangent_mode < simcoon::tangent_none || ctrl.tangent_mode > simcoon::tangent_algorithmic) {
        throw std::invalid_argument("solver: tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
                                    + std::to_string(ctrl.tangent_mode) + " (3 = closest-point is reserved)");
    }
    // tabular (mode 3) steps carry an ABSOLUTE time column: repeating them (ncycle > 1)
    // is ill-defined (and generate() consumes the '2' hold flags on the first pass)
    for (const auto &bl : blocks) {
        if (bl.ncycle > 1) {
            for (const auto &sptr : bl.steps) {
                if (sptr->mode == 3) {
                    throw simcoon::exception_solver("block " + std::to_string(bl.number) + ": tabular (mode 3) steps cannot be cycled (ncycle > 1); unroll the cycles into explicit steps");
                }
            }
        }
    }

    // aliases keep the extracted historical solver() body below textually identical
    const double &div_tnew_dt_solver = ctrl.div_tnew_dt;
    const double &mul_tnew_dt_solver = ctrl.mul_tnew_dt;
    const int &miniter_solver = ctrl.miniter;
    const int &maxiter_solver = ctrl.maxiter;
    const int &inforce_solver = ctrl.inforce;
    const double &precision_solver = ctrl.precision;
    const double &lambda_solver = ctrl.lambda;
    const int &tangent_mode = ctrl.tangent_mode;

	///Usefull UMAT variables
	int ndi = 3;
	int nshr = 3;
    phase_characteristics rve;  // Representative volume element

    unsigned int size_meca = 0; //6 for small perturbation, 9 for finite deformation
	bool start = true;
	double Time = 0.;
	double DTime = 0.;
    double tnew_dt = 1.;

    mat C = zeros(6,6); //Stiffness dS/dE
    mat c = zeros(6,6); //stifness dtau/deps
    mat DR = eye(3,3);
    mat R = eye(3,3);

//    mat dSdE = zeros(6,6);
//    mat dSdT = zeros(1,6);
    mat dQdE = zeros(6,1);
    mat dQdT = zeros(1,1);

    ///Material properties
    rve.sptr_matprops->update(0, umat_name, 1, psi_rve, theta_rve, phi_rve, props.n_elem, props);

    //Output
    int o_ncount = 0;
    double o_tcount = 0.;

    double error = 0.;
    vec residual;
    vec Delta;
    int nK = 0; // The size of the problem to solve
    mat K;
    mat invK;
    int compteur = 0.;
    
    int inc = 0.;
    double tinc=0.;
    double Dtinc=0.;
    double Dtinc_cur=0.;
    double q_conv = 0.;        //q_conv parameter for 0D convexion, Q_conv = qconv (T-T_init), with q_conv = rho*c_p\tau, tau being a time constant for convexion thermal mechanical conditions
    
    /// Block loop
    for(unsigned int i = 0 ; i < blocks.size() ; i++){

        switch(blocks[i].type) {
            case 1: { //Mechanical
                
                /// resize the problem to solve
                residual = zeros(6);
                Delta = zeros(6);
                K = zeros(6,6);
                invK = zeros(6,6);

                if(blocks[i].control_type <= 4) {
                    size_meca = 6;
                }
                else {
                    size_meca = 9;
                }
//                if((blocks[i].control_type == 1)||(blocks[i].control_type == 2))
//                    size_meca = 6;
//                else if(blocks[i].control_type == 3)
//                    size_meca = 9;
                
                shared_ptr<state_variables_M> sv_M;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    natural_basis nb;
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                else {
                    //sv_M is reassigned properly
                    sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                }
                sv_M->tangent_mode = tangent_mode;
                sv_M->L = zeros(6,6);
                sv_M->Lt = zeros(6,6);
                
                //At start, the rotation increment is null
                DTime = 0.;
                sv_M->DEtot = zeros(6);
                sv_M->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_M(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, corate_type, tnew_dt);
                
                shared_ptr<step_meca> sptr_meca;
                if(solver_type == 1) {
                    //RNL
                    sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[0].steps[0]);
                    assert(blocks[i].control_type == 1);
                    sptr_meca->generate(Time, sv_M->Etot, sv_M->sigma, sv_M->T);
                    
                    Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);

                    //jacobian inversion
                    bool inv_success = inv(invK, K);
                    if (!inv_success) {
                        throw simcoon::exception_solver("Singular Jacobian matrix during solver initialization.");
                    }
                }
                else if ((solver_type < 0)||(solver_type > 2)) {
                    cout << "Error, the solver type is not properly defined";
                    return 1;
                }

                if(start) {
                    sink.init(rve);
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(corate_type); //DEtot = 0 and DT = 0 and DR = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++){
                    
                        sptr_meca = std::dynamic_pointer_cast<step_meca>(blocks[i].steps[j]);
                        if (blocks[i].control_type == 1) {
                            sptr_meca->generate(Time, sv_M->Etot, sv_M->sigma, sv_M->T);
                        }
                        else if (blocks[i].control_type == 2) {
                            sptr_meca->generate(Time, sv_M->Etot, sv_M->PKII, sv_M->T);
                        }
                        else if (blocks[i].control_type == 3) {
                            sptr_meca->generate(Time, sv_M->etot, sv_M->sigma, sv_M->T);
//                            sptr_meca->generate(Time, sv_M->etot, sv_M->tau, sv_M->T);
                        }
                        else if (blocks[i].control_type == 4) {
                            vec Biot_vec = t2v_stress(sv_M->Biot_stress());
                            sptr_meca->generate(Time, t2v_strain(sv_M->U0), Biot_vec, sv_M->T); 
                        }                        
                        else if((blocks[i].control_type == 5)||(blocks[i].control_type == 6)) {
                            sptr_meca->generate_kin(Time, sv_M->F0, sv_M->T);
                        }
                        else {
                            throw simcoon::exception_solver("error in Simulation/Solver/solver.cpp: control_type should be a int value in a range of 1 to 6");
                        }
                    
                        nK = sum(sptr_meca->cBC_meca);
                        
                        inc = 0;
                        while(inc < sptr_meca->ninc) {
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if(sptr_meca->cBC_meca(k)) {
                                        sptr_meca->mecas(inc,k) -= residual(k);
                                    }
                                }
                            }
                            
                            while (tinc<1.) {
                                
                                try {
                                sptr_meca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur, inforce_solver);
                                                             
                                if(nK == 0){
                                    
                                    if (blocks[i].control_type == 1) {
                                        sv_M->DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        sv_M->DR = eye(3,3);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                    }
                                    else if (blocks[i].control_type == 2) {
                                        sv_M->DEtot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm
                                        DTime = Dtinc*sptr_meca->times(inc);
                                        mat HW_inv;
                                        bool inv_success = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                        if (!inv_success) {
                                            throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (control_type=2).");
                                        }
                                        DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);

                                        sv_M->F0 = ER_to_F(v2t_strain(sv_M->Etot), sptr_meca->BC_R);
                                        sv_M->F1 = ER_to_F(v2t_strain(sv_M->Etot + sv_M->DEtot), sptr_meca->BC_R*DR);
                                        
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);

                                        sv_M->Detot = t2v_strain(Delta_log_strain_corate(sv_M->F0, sv_M->F1, sv_M->DR, D, Omega, DTime, corate_type));
                                        //mat e_tot_log = t2v_strain(0.5*logmat_sympd(L_Cauchy_Green(sv_M->F1)));
                                        //mat E_dot2 = (1./DTime)*v2t_strain(sv_M->DEtot);
                                    }
                                    else if (blocks[i].control_type == 3) {
                                        sv_M->Detot = Dtinc*sptr_meca->mecas.row(inc).t();
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm
                                        DTime = Dtinc*sptr_meca->times(inc);

                                        mat HW_inv;
                                        bool inv_success = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                        if (!inv_success) {
                                            throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (control_type=3).");
                                        }
                                        DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);

                                        sv_M->F0 = eR_to_F(v2t_strain(sv_M->etot), sptr_meca->BC_R);
                                        sv_M->F1 = eR_to_F(v2t_strain(sv_M->etot + sv_M->Detot), sptr_meca->BC_R*DR);

                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);

                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;
                                        
                                        if (DTime > simcoon::iota)
                                            D = sv_M->Detot/DTime;
                                        else
                                            D = zeros(3,3);
                                    }
                                    else if (blocks[i].control_type == 4) {

                                        sv_M->U1 = sv_M->U0 + v2t_strain(Dtinc*sptr_meca->mecas.row(inc).t());
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        //Application of the Hughes-Winget (1980) algorithm

                                        DTime = Dtinc*sptr_meca->times(inc);
                                        mat HW_inv;
                                        bool inv_success = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                        if (!inv_success) {
                                            throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (control_type=4).");
                                        }
                                        DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);

                                        sv_M->F0 = sptr_meca->BC_R*sv_M->U0;
                                        sv_M->F1 = (sptr_meca->BC_R*DR)*(sv_M->U1);
                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;
                                                                                
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);

                                        sv_M->Detot = t2v_strain(Delta_log_strain_corate(sv_M->F0, sv_M->F1, sv_M->DR, D, Omega, DTime, corate_type));
                                    }                                    
                                    else {
                                        sv_M->F1 = v2t(sptr_meca->BC_mecas.row(inc).t());
                                        sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                        DTime = Dtinc*sptr_meca->times(inc);
                                    
                                        mat D = zeros(3,3);
                                        mat Omega = zeros(3,3);
                                        corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);
                                        sv_M->Detot = t2v_strain(Delta_log_strain_corate(sv_M->F0, sv_M->F1, sv_M->DR, D, Omega, DTime, corate_type));
                                        sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;
                                    }
                                    rve.to_start();
                                    run_umat_M(rve, sv_M->DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, corate_type, tnew_dt);
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                
                                    error = 1.;
                                    
                                    if (blocks[i].control_type == 1) {
                                    
                                        sv_M->DEtot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 2) {
                                        sv_M->DEtot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->PKII(k) - sv_M->PKII_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 3) {
                                        sv_M->Detot = zeros(6);
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_M->Detot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }
                                    else if (blocks[i].control_type == 4) {
                                        vec Biot_stress = t2v_stress(sv_M->Biot_stress());
                                        vec Biot_stress_start = t2v_stress(sv_M->Biot_stress_start());
                                        vec DU_vec = t2v_strain(sv_M->U1 - sv_M->U0);                                        
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_meca->cBC_meca(k)) {
                                                residual(k) = Biot_stress(k) - Biot_stress_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(DU_vec(k) - Dtinc*sptr_meca->mecas(inc,k));
                                            }
                                        }
                                    }                                    
                                    else {
                                        throw simcoon::exception_solver("error, control types 5 and 6 are intended for use in strain-controlled loading only");
                                    }
                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {

                                        if(solver_type != 1){
                                            // classic
                                            ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                            //we use the ddsdde (Lt) from the previous increment
                                            if (blocks[i].control_type == 1) {
                                                Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);
                                            }
                                            else if (blocks[i].control_type == 2) {
                                                // ct2: residual on PKII(E) -> tangent dS/dE via DtauDe_corate_2_DSDE (corate-matched
                                                // box-tangent pull-back; see objective_rates.hpp). NOT Dsigma_LieDD_2_DSDE (spurious J + wrong spin).
                                                C = DtauDe_corate_2_DSDE(sv_M->Lt, corate_type, sv_M->F1, v2t_stress(sv_M->tau));
                                                Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                            }
                                            else if (blocks[i].control_type == 3) {
                                                // ct3: residual on Kirchhoff tau (conjugate to ln V) -> Jacobian is the box tangent Lt directly.
                                                Lt_2_K(sv_M->Lt, K, sptr_meca->cBC_meca, lambda_solver);
                                            }
                                            else if (blocks[i].control_type == 4) {
                                                // ct4: residual on Biot stress -> chain dS/dE (corate-matched) into d(Biot)/dU via DSDE_DBiotStressDU.
                                                mat DSDE = DtauDe_corate_2_DSDE(sv_M->Lt, corate_type, sv_M->F1, v2t_stress(sv_M->tau));
                                                mat R = zeros(3,3);
                                                mat U = zeros(3,3);
                                                RU_decomposition(R,U,sv_M->F1);
                                                C = DSDE_DBiotStressDU(DSDE, U, v2t_stress(sv_M->PKII));
                                                Lt_2_K(C, K, sptr_meca->cBC_meca, lambda_solver);
                                            }

                                            ///jacobian inversion
                                            bool inv_success = inv(invK, K);
                                            if (!inv_success) {
                                                // Degenerate tangent at the trial state (e.g. a
                                                // branch-flip excursion under stress control): bisect
                                                // the increment; at the minimal fraction fall through
                                                // to the existing inforce path rather than throwing
                                                // (a hard throw here overrides the inforce contract and
                                                // is platform-fragile — LAPACK-backend-dependent
                                                // singularity detection).
                                                try_step_cut(Dtinc_cur, sptr_meca->Dn_mini, div_tnew_dt_solver, tnew_dt);
                                                compteur = maxiter_solver;
                                                break;
                                            }

                                            /// Prediction of the component of the strain tensor
                                            Delta = -invK * residual;
                                        }
                                        else if(solver_type == 1) {
                                            //RNL
                                            vec sigma_in_red = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    sigma_in_red(k) = sv_M->sigma_in(k) - sv_M->sigma_in_start(k);
                                                }
                                                else {
                                                    sigma_in_red(k) = 0.;
                                                }
                                            }
                                            Delta = -invK * residual;
                                        }

                                        if (blocks[i].control_type == 1) {
                                            sv_M->DR = eye(3,3);
                                            sv_M->DEtot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            DTime = Dtinc*sptr_meca->times(inc);
                                        }
                                        else if (blocks[i].control_type == 2) {
                                        
                                            sv_M->DEtot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);
                                            mat HW_inv;
                                            bool inv_success_hw = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                            if (!inv_success_hw) {
                                                throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (NR iteration, control_type=2).");
                                            }
                                            DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);

                                            sv_M->F0 = ER_to_F(v2t_strain(sv_M->Etot), sptr_meca->BC_R);
                                            sv_M->F1 = ER_to_F(v2t_strain(sv_M->Etot + sv_M->DEtot), sptr_meca->BC_R*DR);
                                    
                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);
                                            sv_M->Detot = t2v_strain(Delta_log_strain_corate(sv_M->F0, sv_M->F1, sv_M->DR, D, Omega, DTime, corate_type));
                                        }
                                        else if (blocks[i].control_type == 3) {

                                            sv_M->Detot += Delta;
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);
                                            mat HW_inv;
                                            bool inv_success_hw = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                            if (!inv_success_hw) {
                                                throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (NR iteration, control_type=3).");
                                            }
                                            DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);

                                            sv_M->F0 = eR_to_F(v2t_strain(sv_M->etot), sptr_meca->BC_R);
                                            sv_M->F1 = eR_to_F(v2t_strain(sv_M->etot + sv_M->Detot), sptr_meca->BC_R*DR);

                                            sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;

                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);
                                            if (DTime > simcoon::iota)
                                                D = sv_M->Detot/DTime;
                                            else
                                                D = zeros(3,3);
                                        }
                                        else if (blocks[i].control_type == 4) {
                                    
                                            sv_M->U1 += v2t_strain(Delta);
                                            sv_M->DT = Dtinc*sptr_meca->Ts(inc);
                                            //Application of the Hughes-Winget (1980) algorithm
                                            DTime = Dtinc*sptr_meca->times(inc);
                                            mat HW_inv;
                                            bool inv_success_hw = inv(HW_inv, eye(3,3)-0.5*DTime*sptr_meca->BC_w);
                                            if (!inv_success_hw) {
                                                throw simcoon::exception_solver("Singular matrix in Hughes-Winget rotation update (NR iteration, control_type=4).");
                                            }
                                            DR = HW_inv*(eye(3,3) + 0.5*sptr_meca->BC_w*DTime);
                                            sv_M->F0 = sptr_meca->BC_R*sv_M->U0;
                                            sv_M->F1 = (sptr_meca->BC_R*DR)*sv_M->U1;   // BC_R*DR (right): consistent with the accumulator BC_R = BC_R*DR (step_meca) and the ct4 nK==0 / ct2 / ct3 sites; was DR*BC_R
                                            sv_M->DEtot = t2v_strain(Green_Lagrange(sv_M->F1)) - sv_M->Etot;;
                                            mat D = zeros(3,3);
                                            mat Omega = zeros(3,3);
                                            corate_kinematics(corate_type, sv_M->DR, D, Omega, sv_M->F0, sv_M->F1, DTime);
                                            sv_M->Detot = t2v_strain(Delta_log_strain_corate(sv_M->F0, sv_M->F1, sv_M->DR, D, Omega, DTime, corate_type));
                                        }      
                                        rve.to_start();
                                        run_umat_M(rve, sv_M->DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, corate_type, tnew_dt);

                                        if (blocks[i].control_type == 1) {
                                        
                                            //sv_M->DEtot = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 2) {
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = sv_M->PKII(k) - sv_M->PKII_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->DEtot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 3) {
                                            //sv_M->DEtot = zeros(6);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = sv_M->tau(k) - sv_M->tau_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(sv_M->Detot(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }
                                        else if (blocks[i].control_type == 4) {

                                            vec Biot_stress = t2v_stress(sv_M->Biot_stress());
                                            vec Biot_stress_start = t2v_stress(sv_M->Biot_stress_start());                                        
                                            vec DU_vec = t2v_strain(sv_M->U1 - sv_M->U0);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_meca->cBC_meca(k)) {
                                                    residual(k) = Biot_stress(k) - Biot_stress_start(k) - Dtinc*sptr_meca->mecas(inc,k);
                                                }
                                                else {
                                                    residual(k) = lambda_solver*(DU_vec(k) - Dtinc*sptr_meca->mecas(inc,k));
                                                }
                                            }
                                        }                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_meca->Dn_mini) > simcoon::iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
/*                                if((fabs(Dtinc_cur - sptr_meca->Dn_mini) < simcoon::iota)&&(tnew_dt < 1.)) {
//                                    cout << "The subroutine has required a step reduction lower than the minimal indicated at" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                
                                if((error > 1000.*precision_solver)&&(Dtinc_cur == sptr_meca->Dn_mini)) {
//                                    cout << "The error has exceeded 100 times the precision, the simulation has stopped at " << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }*/
                                
                                if(error > precision_solver) {
                                    if(Dtinc_cur == sptr_meca->Dn_mini) {
                                        if(inforce_solver == 1) {
                                            
                                            cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_meca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
//                                            cout << "The next increment has integrated the error to avoid propagation\n";
                                            //The solver has been inforced!
                                            tnew_dt = 1.;
                                            
                                            if (inc+1<sptr_meca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_meca->cBC_meca(k)) {
                                                        sptr_meca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                }
                                            }
                                        }
                                        else if (inforce_solver == 2) {
                                            tnew_dt = 1.;
                                            
                                            if (inc+1<sptr_meca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_meca->cBC_meca(k)) {
                                                        sptr_meca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                }
                                            }
                                        }
                                        else return 1;

                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }

                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                }
                                // Recoverable failures of the trial state — F-reconstruction
                                // (sqrtmat/expmat), singular inverse in the kinematics or the
                                // tangent assembly, singular return-map Jacobian (FB det) —
                                // are bisected away; anything else propagates.
                                catch (const simcoon::exception_sqrtmat_sympd &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_meca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_expmat_sym &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_meca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_inv &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_meca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_det &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_meca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }

                                sptr_meca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime, DR, corate_type);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {

                                sink.record(rve, so, i, n, j, inc, Time);

                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }

                            tinc = 0.;
                            inc++;
                         }

                    }

                }
                break;
            }
            case 2: { //Thermomechanical
                
                /// resize the problem to solve
                residual = zeros(7);
                Delta = zeros(7);
                K = zeros(7,7);
                invK = zeros(7,7);
                
                shared_ptr<state_variables_T> sv_T;
                
                if(start) {
                    rve.construct(0,blocks[i].type);
                    natural_basis nb;
                    rve.sptr_sv_global->update(zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), zeros(6), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), eye(3,3), T_init, 0., nstatev, zeros(nstatev), zeros(nstatev), nb);
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                else {
                    //Dynamic cast from some other (possible state_variable_M/T)
                    /*sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
                     rve.construct(0,blocks[i].type);
                     rve.sptr_sv_global->update(sv_M->Etot, sv_M->DEtot, sv_M->sigma, sv_M->sigma_start, sv_M->T, sv_M->DT, sv_M->sse, sv_M->spd, nstatev, sv_M->statev, sv_M->statev_start);*/
                    //sv_M is reassigned properly
                    sv_T = std::dynamic_pointer_cast<state_variables_T>(rve.sptr_sv_global);
                }
                sv_T->tangent_mode = tangent_mode;
                
                sv_T->dSdE = zeros(6,6);
                sv_T->dSdT = zeros(6,1);
                dQdE = zeros(1,6);
                dQdT = zeros(1,1);
                
                DR = eye(3,3);
                DTime = 0.;
                sv_T->DEtot = zeros(6);
                sv_T->DT = 0.;
                
                //Run the umat for the first time in the block. So that we get the proper tangent properties
                run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                
                sv_T->Q = -1.*sv_T->r;    //Since DTime=0;
                dQdT = lambda_solver;  //To avoid any singularity in the system                
                
                shared_ptr<step_thermomeca> sptr_thermomeca;
                if(solver_type == 1) {
                    //RNL
                    sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[0].steps[0]);
                    sptr_thermomeca->generate(Time, sv_T->Etot, sv_T->sigma, sv_T->T);

                    Lth_2_K(sv_T->dSdE, sv_T->dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);

                    //jacobian inversion
                    bool inv_success = inv(invK, K);
                    if (!inv_success) {
                        throw simcoon::exception_solver("Singular Jacobian matrix during thermomechanical solver initialization.");
                    }
                }
                else if ((solver_type < 0)||(solver_type > 2)) {
                    cout << "Error, the solver type is not properly defined";
                    return 1;
                }

                if(start) {
                    sink.init(rve);
                }
                //Set the start values of sigma_start=sigma and statev_start=statev for all phases
                rve.set_start(corate_type); //DEtot = 0 and DT = 0 so we can use it safely here
                start = false;
                
                /// Cycle loop
                for(unsigned int n = 0; n < blocks[i].ncycle; n++){
                    
                    /// Step loop
                    for(unsigned int j = 0; j < blocks[i].nstep; j++){
                        
                        
                        shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(blocks[i].steps[j]);
                        sptr_thermomeca->generate(Time, sv_T->Etot, sv_T->sigma, sv_T->T);
                        
                        nK = sum(sptr_thermomeca->cBC_meca);
                        
                        inc = 0;
                        if(sptr_thermomeca->cBC_T == 3)
                            q_conv = sptr_thermomeca->BC_T;
                        
                        while(inc < sptr_thermomeca->ninc) {
                            
                            
                            if(error > precision_solver) {
                                for(int k = 0 ; k < 6 ; k++)
                                {
                                    if (sptr_thermomeca->cBC_meca(k)) {
                                        sptr_thermomeca->mecas(inc,k) -= residual(k);
                                    }
                                }
                                if (sptr_thermomeca->cBC_T) {
                                    sptr_thermomeca->Ts(inc) -= residual(6);
                                }
                            }
                            
                            while (tinc<1.) {

                                try {
                                sptr_thermomeca->compute_inc(tnew_dt, inc, tinc, Dtinc, Dtinc_cur, inforce_solver);
                                
                                if(nK + sptr_thermomeca->cBC_T == 0){
                                    
                                    sv_T->DEtot = Dtinc*sptr_thermomeca->mecas.row(inc).t();
                                    sv_T->DT = Dtinc*sptr_thermomeca->Ts(inc);
                                    DTime = Dtinc*sptr_thermomeca->times(inc);
                                    
                                    run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                                    sv_T->Q = -1.*sv_T->r;
                                    
                                }
                                else{
                                    /// ********************** SOLVING THE MIXED PROBLEM NRSTRUCT ***********************************
                                    ///Saving stress and stress set point at the beginning of the loop
                                    
                                    error = 1.;
                                    
                                    sv_T->DEtot = zeros(6);
                                    sv_T->DT = 0.;
                                    
                                    //Construction of the initial residual
                                    for(int k = 0 ; k < 6 ; k++)
                                    {
                                        if (sptr_thermomeca->cBC_meca(k)) {
                                            residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                        }
                                        else {
                                            residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                        }
                                    }
                                    if (sptr_thermomeca->cBC_T == 1) {
                                        residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                    }
                                    else if(sptr_thermomeca->cBC_T == 0) {
                                        residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                    }
                                    else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                        residual(6) = sv_T->Q + q_conv*(sv_T->T-T_init);
                                    }
                                    else {
                                        cout << "error : The Thermal BC is not recognized\n";
                                        return 1;
                                    }

                                    while((error > precision_solver)&&(compteur < maxiter_solver)) {
                                        
                                        if(solver_type != 1){
                                            // classic
                                            ///Prediction of the strain increment using the tangent modulus given from the umat_ function
                                            //we use the ddsdde (Lt) from the previous increment
                                            Lth_2_K(sv_T->dSdE, sv_T->dSdT, dQdE, dQdT, K, sptr_thermomeca->cBC_meca, sptr_thermomeca->cBC_T, lambda_solver);

                                            ///jacobian inversion
                                            bool inv_success = inv(invK, K);
                                            if (!inv_success) {
                                                // Same cut-then-inforce policy as the mechanical loop
                                                // (no throw at Dn_mini — respects the inforce contract).
                                                try_step_cut(Dtinc_cur, sptr_thermomeca->Dn_mini, div_tnew_dt_solver, tnew_dt);
                                                compteur = maxiter_solver;
                                                break;
                                            }

                                            /// Prediction of the component of the strain tensor
                                            Delta = -invK * residual;
                                        }
                                        else if(solver_type == 1) {
                                            //RNL
                                            vec sigma_in_red = zeros(7);
                                            for(int k = 0 ; k < 6 ; k++)
                                            {
                                                if (sptr_thermomeca->cBC_meca(k)) {
                                                    sigma_in_red(k) = sv_T->sigma_in(k) - sv_T->sigma_in_start(k);
                                                }
                                                else {
                                                    sigma_in_red(k) = 0.;
                                                }
                                            }
                                            sigma_in_red(6) = -1.*sv_T->r_in;
                                            Delta = -invK * residual;
                                        }
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            sv_T->DEtot(k) += Delta(k);
                                        }
                                        sv_T->DT += Delta(6);
                                        DTime = Dtinc*sptr_thermomeca->times(inc);
                                        
                                        rve.to_start();
                                        run_umat_T(rve, DR, Time, DTime, ndi, nshr, start, solver_type, blocks[i].control_type, tnew_dt);
                                        
                                        if (DTime < 1.E-12) {
                                            sv_T->Q = -1.*sv_T->r;    //Since DTime=0;
                                            
                                            dQdE = -1.*sv_T->drdE.t();
                                            dQdT = lambda_solver;  //To avoid any singularity in the system
                     
                                        }
                                        else{
                                            //Attention here, we solve Phi = Q - Q_conv = Q - q_conv*(sv_T->T-T_init)) = Q - (rho*c_p*(1./tau)*(sv_T->T-T_init)) = 0
                                            //So the derivative / T has the extra q_conv = rho*c_p*(1./tau)
                                            //It is actually icluded here in dQdT
                                            sv_T->Q = -1.*sv_T->r;
                                            
                                            dQdE = -sv_T->drdE.t();
                                            
                                            if (sptr_thermomeca->cBC_T < 3) {
                                                dQdT = -1.*sv_T->drdT;
                                            }
                                            else if(sptr_thermomeca->cBC_T == 3) {
                                                dQdT = -1.*sv_T->drdT + q_conv;
                                            }
                                            
                                        }
                                                                                
                                        for(int k = 0 ; k < 6 ; k++)
                                        {
                                            if (sptr_thermomeca->cBC_meca(k)) {
                                                residual(k) = sv_T->sigma(k) - sv_T->sigma_start(k) - Dtinc*sptr_thermomeca->mecas(inc,k);
                                            }
                                            else {
                                                residual(k) = lambda_solver*(sv_T->DEtot(k) - Dtinc*sptr_thermomeca->mecas(inc,k));
                                            }
                                        }
                                        if (sptr_thermomeca->cBC_T == 1) {
                                            residual(6) = sv_T->Q - sptr_thermomeca->Ts(inc);
                                        }
                                        else if(sptr_thermomeca->cBC_T == 0) {
                                            residual(6) = lambda_solver*(sv_T->DT - Dtinc*sptr_thermomeca->Ts(inc));
                                        }
                                        else if(sptr_thermomeca->cBC_T == 3) { //Special case of 0D convexion that depends on temperature assumption
                                            residual(6) = sv_T->Q + q_conv*(sv_T->T-T_init);
                                        }
                                        else {
                                            cout << "error : The Thermal BC is not recognized\n";
                                            return 1;
                                        }
                                        
                                        compteur++;
                                        error = norm(residual, 2.);
                                        
                                        if(tnew_dt < 1.) {
                                            if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) > simcoon::iota)||(inforce_solver == 0)) {
                                                compteur = maxiter_solver;
                                            }
                                        }
                                        
                                    }
                                    
                                }
                                
/*                                if((fabs(Dtinc_cur - sptr_thermomeca->Dn_mini) < simcoon::iota)&&(tnew_dt < 1.)) {
                                    cout << "The subroutine has required a step reduction lower than the minimal indicated at" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                
                                if((error > 1000.*precision_solver)&&(Dtinc_cur == sptr_thermomeca->Dn_mini)) {
                                    cout << "The error has exceeded 1000 times the precision, the simulation has stopped at " << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << "\n";
                                    //The solver has been inforced!
                                    return;
                                }
                                */
                                
                                if(error > precision_solver) {
                                    if(Dtinc_cur == sptr_thermomeca->Dn_mini) {
                                        if(inforce_solver == 1) {
                                        
                                            cout << "The solver has been inforced to proceed (Solver issue) at step:" << sptr_thermomeca->number << " inc: " << inc << " and fraction:" << tinc << ", with the error: " << error << "\n";
                                            cout << "The next increment has integrated the error to avoid propagation\n";
                                            //The solver has been inforced!
                                            tnew_dt = 1.;
                                        
                                            if (inc+1<sptr_thermomeca->ninc) {
                                                for(int k = 0 ; k < 6 ; k++)
                                                {
                                                    if(sptr_thermomeca->cBC_meca(k)) {
                                                        sptr_thermomeca->mecas(inc+1,k) -= residual(k);
                                                    }
                                                    if (sptr_thermomeca->cBC_T) {
                                                        sptr_thermomeca->Ts(inc+1) -= residual(6);
                                                    }
                                                    
                                                }
                                            }
                                        }
                                        else return 1;

                                    }
                                    else {
                                        tnew_dt = div_tnew_dt_solver;
                                    }
                                }

                                if((compteur < miniter_solver)&&(tnew_dt >= 1.)) {
                                    tnew_dt = mul_tnew_dt_solver;
                                }
                                compteur = 0;
                                }
                                // Same recoverable-failure set and step-cut policy as the
                                // mechanical Newton loop above.
                                catch (const simcoon::exception_sqrtmat_sympd &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_thermomeca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_expmat_sym &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_thermomeca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_inv &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_thermomeca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }
                                catch (const simcoon::exception_det &) {
                                    step_cut_or_rethrow(Dtinc_cur, sptr_thermomeca->Dn_mini, div_tnew_dt_solver, tnew_dt, compteur);
                                }

                                sptr_thermomeca->assess_inc(tnew_dt, tinc, Dtinc, rve ,Time, DTime, DR, corate_type);
                                //start variables ready for the next increment
                                
                            }
                            
                            //At the end of each increment, check if results should be written
                            if (so.o_type(i) == 1) {
                                o_ncount++;
                            }
                            if (so.o_type(i) == 2) {
                                o_tcount+=DTime;
                            }
                            
                            //Write the results
                            if (((so.o_type(i) == 1)&&(o_ncount == so.o_nfreq(i)))||(((so.o_type(i) == 2)&&(fabs(o_tcount - so.o_tfreq(i)) < 1.E-12)))) {

                                sink.record(rve, so, i, n, j, inc, Time);

                                if (so.o_type(i) == 1) {
                                    o_ncount = 0;
                                }
                                if (so.o_type(i) == 2) {
                                    o_tcount = 0.;
                                }
                            }

                            tinc = 0.;
                            inc++;
                        }

                    }

                }
                break;
            }
            default: {
                cout << "the block type is not defined!\n";
                break;
            }
        }
        //end of blocks loops
    }

    return 0;
}

void solver(const string &umat_name, const vec &props, const unsigned int &nstatev, const double &psi_rve, const double &theta_rve, const double &phi_rve, const int &solver_type, const int &corate_type, const double &div_tnew_dt_solver, const double &mul_tnew_dt_solver, const int &miniter_solver, const int &maxiter_solver, const int &inforce_solver, const double &precision_solver, const double &lambda_solver, const std::string &path_data, const std::string &path_results, const std::string &pathfile, const std::string &outputfile, const int &tangent_mode) {
    if (tangent_mode < simcoon::tangent_none || tangent_mode > simcoon::tangent_algorithmic) {
        throw std::invalid_argument("solver: tangent_mode must be 0 (none), 1 (continuum) or 2 (algorithmic); got "
                                    + std::to_string(tangent_mode) + " (3 = closest-point is reserved)");
    }

    //Check if the required directories exist:
    if(!filesystem::is_directory(path_data)) {
        cout << "error: the folder for the data, " << path_data << ", is not present" << endl;
        return;
    }
    if(!filesystem::is_directory(path_results)) {
        cout << "The folder for the results, " << path_results << ", is not present and has been created" << endl;
        filesystem::create_directory(path_results);
    }

    std::string ext_filename = outputfile.substr(outputfile.length()-4,outputfile.length());
    std::string filename = outputfile.substr(0,outputfile.length()-4); //to remove the extension

    std::string outputfile_global = filename + "_global" + ext_filename;
    std::string outputfile_local = filename + "_local" + ext_filename;

    std::string output_info_file = "output.dat";

    std::vector<block> blocks;  //loading blocks
    double T_init = 0.;

    //Read the loading path
    read_path(blocks, T_init, path_data, pathfile);

    solver_output so(blocks.size());
    read_output(so, blocks.size(), nstatev, path_data, output_info_file);

    //Check output and step files
    check_path_output(blocks, so);

    solver_params ctrl;
    ctrl.div_tnew_dt = div_tnew_dt_solver;
    ctrl.mul_tnew_dt = mul_tnew_dt_solver;
    ctrl.miniter = miniter_solver;
    ctrl.maxiter = maxiter_solver;
    ctrl.inforce = inforce_solver;
    ctrl.precision = precision_solver;
    ctrl.lambda = lambda_solver;
    ctrl.tangent_mode = tangent_mode;

    solver_file_sink sink(path_results, outputfile_global, outputfile_local);
    //status intentionally ignored: the historical file-driven solver() returned void on early aborts
    solver_run(blocks, T_init, so, umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, corate_type, ctrl, sink);
}
    
} //namespace simcoon
