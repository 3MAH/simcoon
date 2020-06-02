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

///@file Tconstitutive.cpp
///@brief Test for Constitutive tensors in Voigt notation
///@version 1.0

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "EPKCP"
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Simulation/Solver/read.hpp>
#include <simcoon/Simulation/Solver/block.hpp>
#include <simcoon/Simulation/Solver/step.hpp>
#include <simcoon/Simulation/Solver/step_meca.hpp>
#include <simcoon/Simulation/Solver/step_thermomeca.hpp>
#include <simcoon/Simulation/Solver/solver.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

BOOST_AUTO_TEST_CASE( EPKCP_solver )
{
    string path_data = "data";
    string path_results = "results";
    string outputfile = "results_job.txt";
    string pathfile = "path.txt";
    string materialfile = "material.dat";
    string sol_essentials = "solver_essentials.inp";
    string sol_control = "solver_control.inp";
    
    string umat_name;
    unsigned int nprops = 0;
    unsigned int nstatev = 0;
    vec props;
    
    double psi_rve = 0.;
    double theta_rve = 0.;
    double phi_rve = 0.;
    
    int solver_type = 0;
    double div_tnew_dt_solver = 0.;
    double mul_tnew_dt_solver = 0.;
    int miniter_solver = 0;
    int maxiter_solver = 0;
    int inforce_solver = 0;
    double precision_solver = 0.;
    double lambda_solver = 0.;
    
    solver_essentials(solver_type, path_data, sol_essentials);
    solver_control(div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, sol_control);
    
    read_matprops(umat_name, nprops, props, nstatev, psi_rve, theta_rve, phi_rve, path_data, materialfile);
    solver(umat_name, props, nstatev, psi_rve, theta_rve, phi_rve, solver_type, div_tnew_dt_solver, mul_tnew_dt_solver, miniter_solver, maxiter_solver, inforce_solver, precision_solver, lambda_solver, path_data, path_results, pathfile, outputfile);
    
    string path_comparison = "comparison/results_job_global-0.txt";
    string path_outputfile = path_results + "/" + "results_job_global-0.txt";
    
    mat C;
    C.load(path_comparison);
    mat R;
    R.load(path_outputfile);
    
    for (unsigned int i=0; i<C.n_rows; i++) {
        for (unsigned int j=0; j<C.n_cols; j++) {
                BOOST_CHECK( fabs(C(i,j) - R(i,j)) < 1.E-6 );
        }
    }
}
