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

///@file Elastic_orthotropic.cpp
///@brief User subroutine for ortothropic elastic materials in 3D case

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Umat/Finite/hypoelastic_orthotropic.hpp>

using namespace std;
using namespace arma;

namespace simcoon{

///@brief The hypoelastic orthotropic UMAT requires 12 constants:
///@brief props[0-2] : 3 Young moduli (Ex, Ey, Ez)
///@brief props[3-5] : 3 Poisson ratios (nu_xy, nu_xz, nu_yz)
///@brief props[6-8] : 3 shear moduli (Gxy, Gxz, Gyz)
///@brief props[9-11] : 3 CTEs (alpha_x, alpha_y, alpha_z)

///@brief No statev is required for thermoelastic constitutive law

void umat_hypoelasticity_ortho(const string &umat_name, const vec &Etot, const vec &DEtot, const mat &F0, const mat &F1, vec &sigma, mat &Lt, mat &L, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT, const double &Time, const double &DTime, double &Wm, double &Wm_r, double &Wm_ir, double &Wm_d, const int &ndi, const int &nshr, const bool &start, double &tnew_dt, const int &corate_type, const int &tangent_mode)
{

    UNUSED(umat_name);
    UNUSED(Etot);
    UNUSED(DR);
    UNUSED(nprops);
    UNUSED(nstatev);
    UNUSED(statev);
    UNUSED(Time);
    UNUSED(DTime);
    UNUSED(nshr);
    UNUSED(tnew_dt);
    
    double T_init = statev(0);
    
    //From the props to the material properties
    double Ex = props(0);
    double Ey = props(1);
    double Ez = props(2);
    double nuxy = props(3);
    double nuxz = props(4);
    double nuyz = props(5);
    double Gxy = props(6);
    double Gxz = props(7);
    double Gyz = props(8);
    double alphax = props(9);
    double alphay = props(10);
    double alphaz = props(11);

    //Elastic stiffness tensor
    L = L_ortho(Ex,Ey,Ez,nuxy,nuxz,nuyz,Gxy,Gxz,Gyz, "EnuG");    
    
    ///@brief Initialization
    if(start)
    {
        T_init = T;
        sigma = zeros(6);
        
        Wm = 0.;
        Wm_r = 0.;
        Wm_ir = 0.;
        Wm_d = 0.;
    }
    
	vec sigma_start = sigma;
	
	//definition of the CTE tensor
	vec alpha = zeros(6);
	alpha(0) = alphax;
	alpha(1) = alphay;
	alpha(2) = alphaz;
    
	//Compute the elastic strain and the related stress	
    vec DEel = DEtot - alpha*DT;
    sigma = el_pred(sigma_start, L, DEel);
    
    // This kernel integrates a corotational CAUCHY rate, so L is dsigma/dD. The consumer reads
    // Lt as the canonical box d(tau_hat)/dDe, and with tau = J sigma and dJ/dt = J tr(D),
    //
    //     d(tau_circ)/dD = J [ d(sigma_circ)/dD + sigma (x) I ],
    //
    // for ANY rate linear in the tensor and its transport -- the J and the stress-proportional
    // term are both needed. `Lt = L` used to hand over the unconverted dsigma/dD: wrong by
    // exactly that, and invisible at J ~ 1, which is why it survived.
    //
    // sigma (x) I is sigma on the STRESS index pair and I on the STRAIN pair, and it is NOT
    // symmetrised here: the 0.5*(X + X^T) in abaqus_jacobian is a concession to Abaqus's
    // symmetric equation solver, not the exact operator. In engineering Voigt that is
    // sigma * I_v^T with I_v = {1,1,1,0,0,0}.
    //
    // No corate conversion: L is already expressed in the corotated frame the solver handed
    // this kernel, so it is in-rate whatever corate_type is -- the same reason the small-strain
    // boxes need none.
    UNUSED(corate_type);
    double J_hypo;
    try {
        J_hypo = det(F1);
    } catch (const std::runtime_error &e) {
        cerr << "Error in det: " << e.what() << endl;
        throw simcoon::exception_det("Error in det function inside umat_hypoelasticity_ortho.");
    }
    const vec I_v = {1., 1., 1., 0., 0., 0.};
    Lt = J_hypo*(L + sigma*I_v.t());
        
    //Computation of the mechanical and thermal work quantities
    Wm += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_r += 0.5*sum((sigma_start+sigma)%DEtot);
    Wm_ir += 0.;
    Wm_d += 0.;
    
    statev(0) = T_init;
}

} //namespace simcoon
