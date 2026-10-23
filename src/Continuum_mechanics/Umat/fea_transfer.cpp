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

///@file fea_transfer.cpp
///@brief Functions to transfer data between FEA software formats and simcoon internal format
///@version 1.0

#include <iostream>
#include <algorithm>
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Maths/rotation.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/kinematics.hpp>
#include <simcoon/Continuum_mechanics/Functions/objective_rates.hpp>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>
#include <simcoon/Simulation/Phase/state_variables.hpp>
#include <simcoon/Simulation/Phase/state_variables_T.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

//=============================================================================
// ABAQUS UMAT TRANSFER FUNCTIONS
//=============================================================================

namespace {

// Abaqus arrays hold the element's ntens components in order (ndi direct, then nshr shears 12,
// 13, 23): entry i is the simcoon component K(i). The others have zero strain (ndi = 3) or are
// stress-free, solved by the material point (ndi < 3, condensed_state.hpp)
struct abaqus_components {
    int ndi, ntens;
    abaqus_components(const int &ndi_, const int &nshr) : ndi(ndi_), ntens(ndi_ + nshr) {}
    int operator()(const int &i) const { return (i < ndi) ? i : 3 + i - ndi; }
};

// Abaqus vector -> simcoon 6 entries (vec, 6x1 or 1x6 mat), zero outside K
template <class V> void abaqus2vec(const double *a, const abaqus_components &K, V &v) {
    v.zeros();
    for (int i = 0; i < K.ntens; i++) v(K(i)) = a[i];
}

template <class V> void vec2abaqus(const V &v, const abaqus_components &K, double *a) {
    for (int i = 0; i < K.ntens; i++) a[i] = v(K(i));
}

// column-major ntens x ntens <-> the K x K block of a 6x6
void abaqus2mat(const double *a, const abaqus_components &K, mat &M) {
    M.zeros(6, 6);
    for (int j = 0; j < K.ntens; j++)
        for (int i = 0; i < K.ntens; i++) M(K(i), K(j)) = a[j*K.ntens + i];
}

void mat2abaqus(const mat &M, const abaqus_components &K, double *a) {
    for (int j = 0; j < K.ntens; j++)
        for (int i = 0; i < K.ntens; i++) a[j*K.ntens + i] = M(K(i), K(j));
}

// The stress-free strain (ndi < 3, condensed_state.hpp) is kept in STATEV after the energies,
// in simcoon order. Read at the start of the increment, rotated with it (DROT) as Abaqus rotates
// STRAN: the material point then solves its change, not the whole value
void read_free_strain(const double *statev, const int &offset, const uvec &C, const mat &DR, vec &Etot) {
    if (C.is_empty()) return;
    vec E_C = zeros(6);
    for (uword j = 0; j < C.n_elem; j++) E_C(C(j)) = statev[offset + j];
    Etot(C) = rotate_strain(E_C, DR)(C);
}

void write_free_strain(const vec &Etot_end, const int &offset, const uvec &C, double *statev) {
    if (Etot_end.n_elem != 6) return;
    for (uword j = 0; j < C.n_elem; j++) statev[offset + j] = Etot_end(C(j));
}

}  // namespace

void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, vec &sigma, mat &Lt, vec &Wm, vec &statev_smart)
{
    const abaqus_components K(ndi, nshr);
    abaqus2vec(stress, K, sigma);
    abaqus2mat(ddsdde, K, Lt);

    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+4];
    }
}

void abaqus2smart_M(const double *stress, const double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &statev_smart, mat &DR, bool &start)
{
    const abaqus_components K(ndi, nshr);
    abaqus2vec(stress, K, sigma);
    abaqus2vec(stran, K, Etot);
    abaqus2vec(dstran, K, DEtot);
    abaqus2mat(ddsdde, K, Lt);

    ///@brief rotation matrix
    DR = mat(drot, 3, 3);

    ///@brief Temperature and temperature increment creation
    T = temperature;
    DT = Dtemperature;

    ///@brief Time
    Time = time[1];
    DTime = dtime;

    ///@brief Initialization
    start = (Time < 1E-12);

    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    const uvec C = make_condensed_state(ndi, nshr).C;
    read_free_strain(statev, 4, C, DR, Etot);
    const int offset = 4 + C.n_elem;
    for (int i=0; i<nstatev-offset; i++) {
        statev_smart(i) = statev[i+offset];
    }
}

void abaqus2smart_T(const double *stress, const double *ddsdde, const double *ddsddt, const double *drplde, const double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &Wt, vec &statev_smart, mat &DR, bool &start) {

    const abaqus_components K(ndi, nshr);
    abaqus2vec(stress, K, sigma);
    abaqus2vec(stran, K, Etot);
    abaqus2vec(dstran, K, DEtot);
    abaqus2mat(ddsdde, K, dSdE);
    abaqus2vec(ddsddt, K, dSdT);
    abaqus2vec(drplde, K, drpldE);
    drpldT(0) = drpldt;

    ///@brief rotation matrix
    DR = mat(drot, 3, 3);

    ///@brief Temperature and temperature increment creation
    T = temperature;
    DT = Dtemperature;

    ///@brief Time
    Time = time[1];
    DTime = dtime;

    ///@brief Initialization
    start = (Time < 1E-12);

    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<3; i++) {
        Wt(i) = statev[i+4];
    }
    const uvec C = make_condensed_state(ndi, nshr).C;
    read_free_strain(statev, 7, C, DR, Etot);
    const int offset = 7 + C.n_elem;
    for (int i=0; i<nstatev-offset; i++) {
        statev_smart(i) = statev[i+offset];
    }
}

bool abaqus_nlgeom(const double *dfgrd1) {
    // Abaqus sets DFGRD1 to the identity when the step has no geometric nonlinearity
    const double identity[9] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    for (int i = 0; i < 9; i++) {
        if (dfgrd1[i] != identity[i]) return true;
    }
    return false;
}

mat abaqus_jacobian(const mat &Lt, const vec &sigma) {
    // symmetric part of sigma (x) I: what Abaqus's default (symmetric) solver keeps
    const vec I = {1., 1., 1., 0., 0., 0.};
    const mat sigma_I = sigma * I.t();
    return Lt + 0.5*(sigma_I + sigma_I.t());
}

int abaqus_kinematic_form(const char *cmname) {
    // CMNAME is a space-padded CHARACTER*80 (no terminator); a C string stops at '\0', so
    // each character is read only past a non-terminating one. `_R` must be the whole
    // suffix: any other tail (EPICP_REF, ...) is ignored, as before the suffix existed
    auto end = [](char c) { return c == ' ' || c == '\0'; };
    return (cmname[5] == '_' && (cmname[6] == 'R' || cmname[6] == 'r') && end(cmname[7])) ? 3 : corate_type_total;
}

int abaqus_kinematic_statev(const int &corate_type) {
    return (corate_type == corate_type_total) ? 0 : 15;
}

bool abaqus_finite_route(const std::string &umat_name, const int &ndi, const int &corate_type, const bool &thermomechanical, const double *dfgrd1) {
    if (ndi != 3 || !abaqus_nlgeom(dfgrd1))
        return false;
    if (thermomechanical)   // the finite thermomechanical dispatch serves every law of the small one
        return true;
    const auto &names = finite_umat_names();
    auto it = names.find(umat_name);
    if (it == names.end() || it->second <= 1)   // UMEXT, UMABA: foreign plugins keep their own convention
        return false;
    return corate_type != corate_type_total || total_formulation_serves(umat_name);
}

void abaqus2smart_finite(state_variables &sv, const double *drot, const double *dfgrd0, const double *dfgrd1, const double *kin_statev, const int &corate_type, const double &DTime) {
    sv.F0 = mat(dfgrd0, 3, 3);
    sv.F1 = mat(dfgrd1, 3, 3);
    // Abaqus rotated STRESS by DROT for its own corotational route: back to the global basis
    const mat DROT(drot, 3, 3);
    sv.sigma = det(sv.F0)*rotate_stress(sv.sigma, mat(DROT.t()));
    sv.tau = sv.sigma;
    sv.Etot = t2v_strain(Green_Lagrange(sv.F0));
    sv.DEtot = t2v_strain(Green_Lagrange(sv.F1)) - sv.Etot;
    if (corate_type == corate_type_total) {   // a state function of F
        sv.etot = t2v_strain(0.5*logmat_sympd(mat(sv.F0.t()*sv.F0)));
        sv.R = eye(3, 3);
    }
    else {
        sv.etot = vec(kin_statev, 6);
        sv.R = mat(kin_statev + 6, 3, 3);
        if (!any(vectorise(sv.R)))   // STATEV starts at zero: an undeformed start
            sv.R.eye();
    }
    mat D, Omega;
    const double dt = (DTime > simcoon::iota) ? DTime : 1.;   // D dt is what the increment uses
    corate_kinematics(corate_type, sv.DR, D, Omega, sv.F0, sv.F1, dt);
    sv.Detot = t2v_strain(Delta_log_strain_corate(sv.F0, sv.F1, sv.DR, D, Omega, dt, corate_type));
}

void smart2abaqus_finite_M(const state_variables &sv, const mat &Lt, const int &corate_type, vec &sigma, mat &ddsdde, double *kin_statev) {
    const mat tau = v2t_stress(sv.tau);
    const mat DSDE = DtauDe_corate_2_DSDE(Lt, corate_type, sv.F1, tau);
    sigma = sv.tau/det(sv.F1);
    ddsdde = DSDE_2_Dsigma_JaumannDD(DSDE, sv.F1, v2t_stress(sigma));
    if (corate_type != corate_type_total) {   // the commit of set_start
        const vec etot_end = rotate_strain(sv.etot, sv.DR) + sv.Detot;
        const mat R_end = sv.DR*sv.R;
        std::copy(etot_end.begin(), etot_end.end(), kin_statev);
        std::copy(R_end.begin(), R_end.end(), kin_statev + 6);
    }
}

void smart2abaqus_finite_T(const state_variables_T &sv, const int &corate_type, const double &DTime, vec &sigma, mat &ddsdde, vec &ddsddt, double &rpl, vec &drplde, double &drpldt, double *kin_statev) {
    smart2abaqus_finite_M(sv, sv.dSdE, corate_type, sigma, ddsdde, kin_statev);
    const double J = det(sv.F1);
    const double per_volume_time = (DTime > simcoon::iota) ? 1./(J*DTime) : 0.;
    ddsddt = vectorise(sv.dSdT)/J;
    rpl = sv.r*per_volume_time;
    drpldt = sv.drdT(0, 0)*per_volume_time;
    const vec drdE = box_covector_2_dE(vectorise(sv.drdE), sv.F1, corate_type);
    const vec I = {1., 1., 1., 0., 0., 0.};
    drplde = (t2v_stress(sv.F1*v2t_stress(drdE)*sv.F1.t()) - sv.r*I)*per_volume_time;
}

void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const vec &Wm, const mat &Lt_kernel, const bool &nlgeom, const vec &Etot_end)
{
    // the tangent arrives condensed (select_umat_M): DDSDDE is its K x K block; sigma (x) I only
    // for NLGEOM with 3 direct components, plane stress and 1D take it as it is
    const mat Lt = (nlgeom && ndi == 3) ? abaqus_jacobian(Lt_kernel, sigma) : Lt_kernel;
    const abaqus_components K(ndi, nshr);
    vec2abaqus(sigma, K, stress);
    mat2abaqus(Lt, K, ddsdde);

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    const uvec C = make_condensed_state(ndi, nshr).C;
    write_free_strain(Etot_end, 4, C, statev);
    const int offset = 4 + C.n_elem;
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+offset] = statev_smart(i);
    }
}

void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, double *dstran, double *time, double &dtime, double &temperature, double &Dtemperature, int &nprops, double *props,  int &nstatev, double *statev, const int &ndi, const int &nshr, double *drot, const vec &sigma, const mat &Lt_kernel, const vec &Etot, const vec &DEtot, const double &T, const double &DT, const double &Time, const double &DTime, const vec &props_smart, const vec &Wm, const vec &statev_smart, const mat &DR, bool &start)
{
    UNUSED(nstatev);
    // simcoon -> the INPUT arrays of an Abaqus UMAT (UMABA plugin): ddsdde is overwritten by
    // that UMAT, so the kernel tangent goes out as it is
    const abaqus_components K(ndi, nshr);
    vec2abaqus(sigma, K, stress);
    vec2abaqus(Etot, K, stran);
    vec2abaqus(DEtot, K, dstran);
    mat2abaqus(Lt_kernel, K, ddsdde);

    ///@brief rotation matrix
    std::copy(DR.begin(), DR.end(), drot);

    ///@brief Temperature and temperature increment creation
    temperature = T;
    Dtemperature = DT;

    ///@brief Time
    time[1] = Time;
    dtime = DTime;

    ///@brief Initialization
    start = (Time < 1E-12);

    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props[i] = props_smart(i);
    }

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+4] = statev_smart(i);
    }
}

void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, double &rpl, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const double &r, const vec &Wm, const vec &Wt, const mat &dSdE_kernel, const mat &dSdT, const mat &drpldE, const mat &drpldT, const bool &nlgeom, const vec &Etot_end) {

    // dSdE, dSdT, drpldE and drpldT arrive condensed on the element's components (select_umat_T)
    const mat dSdE = (nlgeom && ndi == 3) ? abaqus_jacobian(dSdE_kernel, sigma) : dSdE_kernel;
    const abaqus_components K(ndi, nshr);
    vec2abaqus(sigma, K, stress);
    mat2abaqus(dSdE, K, ddsdde);
    vec2abaqus(dSdT, K, ddsddt);
    vec2abaqus(drpldE, K, drplde);
    drpldt = drpldT(0,0);
    rpl = r;

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (int i=0; i<3; i++) {
        statev[i+4] = Wt(i);
    }
    const uvec C = make_condensed_state(ndi, nshr).C;
    write_free_strain(Etot_end, 7, C, statev);
    const int offset = 7 + C.n_elem;
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+offset] = statev_smart(i);
    }
}

//=============================================================================
// ANSYS USERMAT TRANSFER FUNCTIONS
//=============================================================================

void ansys2smart_M(const double *stress, const double *dstran,
    const double &sedEl, const double &sedPl, const double &epseq,
    const double *statev, const double *props,
    const double &Time, const double &DTime,
    const double &temperature, const double &Dtemperature,
    const int &ncomp, const int &nprops, const int &nstatev,
    vec &sigma, vec &DEtot, vec &Wm, 
    vec &statev_smart, vec &props_smart,
    double &T, double &DT)
{
    // Ansys uses Voigt notation: (11, 22, 33, 12, 23, 13)
    // simcoon uses:              (11, 22, 33, 12, 13, 23)
    // So components 4 and 5 are swapped
    
    if (ncomp == 6) {  // 3D
        // Stress: swap components 4 and 5
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(2) = stress[2];
        sigma(3) = stress[3];
        sigma(4) = stress[5];  // Ansys 13 → simcoon index 4
        sigma(5) = stress[4];  // Ansys 23 → simcoon index 5
        
        // Strain increment: swap components 4 and 5
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(2) = dstran[2];
        DEtot(3) = dstran[3];
        DEtot(4) = dstran[5];  // Ansys 13 → simcoon index 4
        DEtot(5) = dstran[4];  // Ansys 23 → simcoon index 5
    }
    else if (ncomp == 4) {  // 2D (plane strain/axisymmetric)
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(2) = stress[2];
        sigma(3) = stress[3];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(2) = dstran[2];
        DEtot(3) = dstran[3];
    }
    else if (ncomp == 3) {  // 2D plane stress
        sigma(0) = stress[0];
        sigma(1) = stress[1];
        sigma(3) = stress[2];
        
        DEtot(0) = dstran[0];
        DEtot(1) = dstran[1];
        DEtot(3) = dstran[2];
    }
    else if (ncomp == 1) {  // 1D
        sigma(0) = stress[0];
        DEtot(0) = dstran[0];
    }
    
    // Temperature
    T = temperature;
    DT = Dtemperature;
    
    // Material properties
    for (int i = 0; i < nprops; i++) {
        props_smart(i) = props[i];
    }
    
    // Work quantities from sedEl and sedPl
    Wm(0) = sedEl + sedPl;  // Total strain energy
    Wm(1) = sedEl;          // Elastic (reversible) strain energy
    Wm(2) = sedPl;          // Plastic (irreversible) strain energy
    Wm(3) = 0.;             // Damage energy (not directly available)
    
    // State variables (first 4 are Wm, then user statev)
    int nstatev_smart = nstatev - 4;
    for (int i = 0; i < nstatev_smart; i++) {
        statev_smart(i) = statev[i + 4];
    }
}

void smart2ansys_M(double *stress, double *ddsdde,
    double &sedEl, double &sedPl, double *statev,
    const int &ncomp,
    const vec &sigma, const mat &Lt,
    const vec &Wm, const vec &statev_smart)
{
    // Ansys uses Voigt notation: (11, 22, 33, 12, 23, 13)
    // simcoon uses:              (11, 22, 33, 12, 13, 23)
    // So components 4 and 5 are swapped
    
    if (ncomp == 6) {  // 3D
        // Stress: swap components 4 and 5
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(2);
        stress[3] = sigma(3);
        stress[4] = sigma(5);  // simcoon 23 → Ansys index 4
        stress[5] = sigma(4);  // simcoon 13 → Ansys index 5
        
        // Tangent matrix: need to swap rows 4↔5 and columns 4↔5
        // Create permutation: simcoon indices [0,1,2,3,4,5] → Ansys [0,1,2,3,5,4]
        int perm[6] = {0, 1, 2, 3, 5, 4};
        
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                // Ansys uses row-major storage
                ddsdde[i * 6 + j] = Lt(perm[i], perm[j]);
            }
        }
    }
    else if (ncomp == 4) {  // 2D (plane strain/axisymmetric)
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(2);
        stress[3] = sigma(3);
        
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                ddsdde[i * 4 + j] = Lt(i, j);
            }
        }
    }
    else if (ncomp == 3) {  // 2D plane stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);
        
        // Condense tangent for plane stress
        assert(Lt(2,2) > 0.);
        ddsdde[0] = Lt(0,0) - Lt(0,2)*Lt(2,0)/Lt(2,2);
        ddsdde[4] = Lt(1,1) - Lt(1,2)*Lt(2,1)/Lt(2,2);
        ddsdde[3] = Lt(0,1) - Lt(0,2)*Lt(2,1)/Lt(2,2);
        ddsdde[1] = Lt(1,0) - Lt(1,2)*Lt(2,0)/Lt(2,2);
        ddsdde[6] = Lt(0,3) - Lt(0,2)*Lt(2,3)/Lt(2,2);
        ddsdde[7] = Lt(1,3) - Lt(1,2)*Lt(2,3)/Lt(2,2);
        ddsdde[2] = Lt(3,0) - Lt(3,2)*Lt(2,0)/Lt(2,2);
        ddsdde[5] = Lt(3,1) - Lt(3,2)*Lt(2,1)/Lt(2,2);
        ddsdde[8] = Lt(3,3) - Lt(3,2)*Lt(2,3)/Lt(2,2);
    }
    else if (ncomp == 1) {  // 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    
    // Work quantities
    sedEl = Wm(1);  // Elastic strain energy
    sedPl = Wm(2);  // Plastic strain energy
    
    // State variables
    for (int i = 0; i < 4; i++) {
        statev[i] = Wm(i);
    }
    for (unsigned int i = 0; i < statev_smart.n_elem; i++) {
        statev[i + 4] = statev_smart(i);
    }
}
    
} // namespace simcoon
