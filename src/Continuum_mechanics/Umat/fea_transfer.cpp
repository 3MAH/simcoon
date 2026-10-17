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
#include <assert.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>
#include <simcoon/Continuum_mechanics/Umat/condensed_state.hpp>

using namespace std;
using namespace arma;

namespace simcoon {

//=============================================================================
// ABAQUS UMAT TRANSFER FUNCTIONS
//=============================================================================

namespace {

// The element's ntens components (Abaqus order: ndi direct, then nshr shears 12, 13, 23) are the
// kept set K of the (ndi, nshr) state, in the same order: entry i of an Abaqus array is the
// simcoon component K(i). The other components have zero strain (ndi = 3) or are stress-free
// and condensed by the material point (ndi < 3, select_umat_M / select_umat_T).
uvec abaqus_components(const int &ndi, const int &nshr) {
    return make_condensed_state(ndi, nshr).K;
}

// Abaqus vector (ntens) -> simcoon 6-vector, zero outside K
void abaqus2vec(const double *a, const uvec &K, vec &v) {
    v.zeros(6);
    for (uword i = 0; i < K.n_elem; i++) v(K(i)) = a[i];
}

// 6-entry simcoon vector or 1x6 / 6x1 matrix -> Abaqus vector (ntens)
template <class V> void vec2abaqus(const V &v, const uvec &K, double *a) {
    for (uword i = 0; i < K.n_elem; i++) a[i] = v(K(i));
}

// Abaqus matrix (ntens x ntens, column-major) -> simcoon 6x6, zero outside K x K
void abaqus2mat(const double *a, const uvec &K, mat &M) {
    const uword n = K.n_elem;
    M.zeros(6, 6);
    for (uword j = 0; j < n; j++)
        for (uword i = 0; i < n; i++) M(K(i), K(j)) = a[j*n + i];
}

// simcoon 6x6 -> Abaqus matrix (ntens x ntens, column-major), the K x K block
void mat2abaqus(const mat &M, const uvec &K, double *a) {
    const uword n = K.n_elem;
    for (uword j = 0; j < n; j++)
        for (uword i = 0; i < n; i++) a[j*n + i] = M(K(i), K(j));
}

void abaqus2rotation(const double *drot, mat &DR) {
    DR.set_size(3, 3);
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) DR(i, j) = drot[j*3 + i];
}

}  // namespace

void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, vec &sigma, mat &Lt, vec &Wm, vec &statev_smart)
{
    const uvec K = abaqus_components(ndi, nshr);
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
    // The condensed components arrive at zero: the material point solves them (no state kept)
    const uvec K = abaqus_components(ndi, nshr);
    abaqus2vec(stress, K, sigma);
    abaqus2vec(stran, K, Etot);
    abaqus2vec(dstran, K, DEtot);
    abaqus2mat(ddsdde, K, Lt);

    ///@brief rotation matrix
    abaqus2rotation(drot, DR);

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
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+4];
    }
}

void abaqus2smart_T(const double *stress, const double *ddsdde, const double *ddsddt, const double *drplde, const double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &Wt, vec &statev_smart, mat &DR, bool &start) {

    const uvec K = abaqus_components(ndi, nshr);
    abaqus2vec(stress, K, sigma);
    abaqus2vec(stran, K, Etot);
    abaqus2vec(dstran, K, DEtot);
    abaqus2mat(ddsdde, K, dSdE);
    dSdT.zeros();
    drpldE.zeros();
    for (uword i = 0; i < K.n_elem; i++) {
        dSdT(K(i)) = ddsddt[i];
        drpldE(K(i)) = drplde[i];
    }
    drpldT(0) = drpldt;

    ///@brief rotation matrix
    abaqus2rotation(drot, DR);

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
    for (int i=0; i<nstatev-7; i++) {
        statev_smart(i) = statev[i+7];
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

void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const vec &Wm, const mat &Lt_kernel, const bool &nlgeom)
{
    // The tangent arrives condensed on the element's components (select_umat_M): the K x K
    // block is DDSDDE. The sigma (x) I term belongs to NLGEOM steps and to the 3 direct
    // components: plane stress and 1D take the condensed kernel tangent as it is (the element
    // owns the thickness)
    const mat Lt = (nlgeom && ndi == 3) ? abaqus_jacobian(Lt_kernel, sigma) : Lt_kernel;
    const uvec K = abaqus_components(ndi, nshr);
    vec2abaqus(sigma, K, stress);
    mat2abaqus(Lt, K, ddsdde);

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+4] = statev_smart(i);
    }
}

void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, double *dstran, double *time, double &dtime, double &temperature, double &Dtemperature, int &nprops, double *props,  int &nstatev, double *statev, const int &ndi, const int &nshr, double *drot, const vec &sigma, const mat &Lt_kernel, const vec &Etot, const vec &DEtot, const double &T, const double &DT, const double &Time, const double &DTime, const vec &props_smart, const vec &Wm, const vec &statev_smart, const mat &DR, bool &start)
{
    UNUSED(nstatev);
    // simcoon -> the INPUT arrays of an Abaqus UMAT (UMABA plugin): ddsdde is overwritten by
    // that UMAT, so the kernel tangent goes out as it is
    const uvec K = abaqus_components(ndi, nshr);
    vec2abaqus(sigma, K, stress);
    vec2abaqus(Etot, K, stran);
    vec2abaqus(DEtot, K, dstran);
    mat2abaqus(Lt_kernel, K, ddsdde);

    ///@brief rotation matrix
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) drot[j*3 + i] = DR(i, j);

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

void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, double &rpl, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const double &r, const vec &Wm, const vec &Wt, const mat &dSdE_kernel, const mat &dSdT, const mat &drpldE, const mat &drpldT, const bool &nlgeom) {

    // dSdE, dSdT, drpldE and drpldT arrive condensed on the element's components (select_umat_T)
    const mat dSdE = (nlgeom && ndi == 3) ? abaqus_jacobian(dSdE_kernel, sigma) : dSdE_kernel;
    const uvec K = abaqus_components(ndi, nshr);
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
    for (unsigned int i=0; i<statev_smart.n_elem; i++) {
        statev[i+7] = statev_smart(i);
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
