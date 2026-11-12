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
#include <cctype>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
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

// ==== FE component maps (fea_transfer.hpp) ====

fea_components ansys_components(const int &nDirect, const int &nShear) {
    if (nDirect == 3 && nShear == 3) return {6, {0, 1, 2, 3, 5, 4}};
    if (nDirect == 3 && nShear == 1) return {4, {0, 1, 2, 3, 0, 0}};
    if (nDirect == 2 && nShear == 1) return {3, {0, 1, 3, 0, 0, 0}};
    if (nDirect == 1 && nShear == 0) return {1, {0, 0, 0, 0, 0, 0}};
    throw std::invalid_argument("Ansys USERMAT: nDirect = " + std::to_string(nDirect) + ", nShear = "
                                + std::to_string(nShear) + " is not supported (3D, plane strain, "
                                "axisymmetric, plane stress and 1D are)");
}

void fe2mat(const double *a, const fea_components &K, mat &M) {
    M.zeros(6, 6);
    for (int j = 0; j < K.ntens; j++)
        for (int i = 0; i < K.ntens; i++) M(K(i), K(j)) = a[j*K.ntens + i];
}

void mat2fe(const mat &M, const fea_components &K, double *a) {
    for (int j = 0; j < K.ntens; j++)
        for (int i = 0; i < K.ntens; i++) a[j*K.ntens + i] = M(K(i), K(j));
}

// ==== ABAQUS UMAT TRANSFER FUNCTIONS ====

namespace {

// Abaqus: ndi direct, then the nshr shears 12, 13, 23, the kept set of the condensed state
fea_components abaqus_components(const condensed_state &cs) {
    fea_components c{static_cast<int>(cs.K.n_elem), {0, 0, 0, 0, 0, 0}};
    for (int i = 0; i < c.ntens; i++) c.K[i] = static_cast<int>(cs.K(i));
    return c;
}

// STATEV after nheader: the stress-free strain (simcoon order), then the law; it is rotated with DR
// like the FE strain, so the material point solves its change, not its whole value
void unpack_statev(const double *statev, const int &nstatev, const int &nheader, const uvec &C, const mat &DR, vec &Etot, vec &statev_smart) {
    if (!C.is_empty()) {
        vec E_C = zeros(6);
        for (uword j = 0; j < C.n_elem; j++) E_C(C(j)) = statev[nheader + j];
        Etot(C) = rotate_strain(E_C, DR)(C);
    }
    const int offset = nheader + C.n_elem;
    for (int i = 0; i < nstatev - offset; i++) statev_smart(i) = statev[i + offset];
}

void pack_statev(const vec &Etot_end, const vec &statev_smart, const int &nheader, const uvec &C, double *statev) {
    if (Etot_end.n_elem == 6) {
        for (uword j = 0; j < C.n_elem; j++) statev[nheader + j] = Etot_end(C(j));
    }
    const int offset = nheader + C.n_elem;
    for (uword i = 0; i < statev_smart.n_elem; i++) statev[i + offset] = statev_smart(i);
}

}  // namespace

void abaqus2smart_M_light(const double *stress, const double *ddsdde, const int &nstatev, double *statev, const int &ndi, const int &nshr, vec &sigma, mat &Lt, vec &Wm, vec &statev_smart)
{
    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    fe2vec(stress, K, sigma);
    fe2mat(ddsdde, K, Lt);

    for (int i=0; i<4; i++) {
        Wm(i) = statev[i];
    }
    for (int i=0; i<nstatev-4; i++) {
        statev_smart(i) = statev[i+4];
    }
}

void abaqus2smart_M(const double *stress, const double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &statev_smart, mat &DR, bool &start)
{
    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    fe2vec(stress, K, sigma);
    fe2vec(stran, K, Etot);
    fe2vec(dstran, K, DEtot);
    fe2mat(ddsdde, K, Lt);

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
    unpack_statev(statev, nstatev, 4, cs.C, DR, Etot, statev_smart);
}

void abaqus2smart_T(const double *stress, const double *ddsdde, const double *ddsddt, const double *drplde, const double &drpldt, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops, const double *props, const int &nstatev, double *statev, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &dSdE, vec &dSdT, vec &drpldE, double &drpldT, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &Wm, vec &Wt, vec &statev_smart, mat &DR, bool &start) {

    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    fe2vec(stress, K, sigma);
    fe2vec(stran, K, Etot);
    fe2vec(dstran, K, DEtot);
    fe2mat(ddsdde, K, dSdE);
    fe2vec(ddsddt, K, dSdT);
    fe2vec(drplde, K, drpldE);
    drpldT = drpldt;

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
    unpack_statev(statev, nstatev, 7, cs.C, DR, Etot, statev_smart);
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
    // symmetric part of sigma (x) I: what Abaqus's default (symmetric) solver keeps, built in
    // place (the products and sums of Lt + 0.5*(sigma I^T + I sigma^T), without the temporaries)
    static const double I[6] = {1., 1., 1., 0., 0., 0.};
    mat J(6, 6);
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            const double s1 = sigma(i)*I[j];
            const double s2 = sigma(j)*I[i];
            const double h = (s1 + s2)*0.5;
            J(i, j) = Lt(i, j) + h;
        }
    }
    return J;
}

namespace {
// CMNAME: space-padded CHARACTER*80 (a C string stops at '\0')
std::size_t cmname_length(const char *cmname) {
    std::size_t n = 0;
    while (n < 80 && cmname[n] != '\0' && cmname[n] != ' ') ++n;
    return n;
}
bool cmname_rate_suffix(const char *cmname, const std::size_t n) {
    return n >= 2 && cmname[n - 2] == '_' && (cmname[n - 1] == 'R' || cmname[n - 1] == 'r');
}
}  // namespace

std::string abaqus_umat_name(const char *cmname) {
    const std::size_t n = cmname_length(cmname);
    std::string name(cmname, cmname_rate_suffix(cmname, n) ? n - 2 : n);
    // Abaqus upper-cases material names: the _closed suffix is matched in any case
    static const std::string closed = "_closed";
    bool closed_suffix = false;
    if (name.size() > closed.size()) {
        const std::size_t k = name.size() - closed.size();
        closed_suffix = true;
        for (std::size_t i = 0; i < closed.size() && closed_suffix; ++i)
            closed_suffix = std::tolower(static_cast<unsigned char>(name[k + i])) == closed[i];
        if (closed_suffix) name.replace(k, closed.size(), closed);
    }
    // a descriptive tail on a law name (EPICP_STEEL): the first five letters, said once per name
    if (name.size() > 5 && !closed_suffix && umat_entry_of(name) == nullptr) {
        static std::mutex said_mutex;
        static std::set<std::string> said;
        const std::lock_guard<std::mutex> lock(said_mutex);
        if (said.insert(name).second) {
            std::cerr << "simcoon UMAT: CMNAME '" << name << "' names no law; its first five letters '"
                      << name.substr(0, 5) << "' select it\n";
        }
        name.resize(5);
    }
    return name;
}

int abaqus_kinematic_form(const char *cmname) {
    return cmname_rate_suffix(cmname, cmname_length(cmname)) ? 3 : corate_type_total;
}

int abaqus_kinematic_statev(const int &corate_type) {
    return (corate_type == corate_type_total) ? 0 : 15;
}

int fea_law_statev(const int &nstatev, const int &nheader, const int &ndi, const int &nshr, const int &corate_type) {
    const int nbridge = nheader + static_cast<int>(make_condensed_state(ndi, nshr).C.n_elem) + abaqus_kinematic_statev(corate_type);
    if (nstatev < nbridge) {
        throw std::invalid_argument("the STATEV size (Abaqus DEPVAR, Ansys TB,STATE) " + std::to_string(nstatev) + " is below the "
                                    + std::to_string(nbridge) + " entries the bridge keeps"
                                    + ((corate_type == corate_type_total) ? std::string() : std::string(" (the _R suffix takes 15)"))
                                    + ", before the law's own");
    }
    return nstatev - nbridge;
}

bool fea_numerical_failure(const std::exception &e) {
    return dynamic_cast<const std::runtime_error *>(&e) != nullptr
        && dynamic_cast<const simcoon::exception_solver *>(&e) == nullptr;
}

bool abaqus_finite_route(const std::string &umat_name, const int &ndi, const int &corate_type, const bool &thermomechanical, const double *dfgrd1) {
    (void)ndi;   // a condensed state solves its free stretches on the finite route (condensed_F)
    if (!abaqus_nlgeom(dfgrd1))
        return false;
    if (thermomechanical)   // the finite thermomechanical dispatch serves every law of the small one under log_R
        return corate_type != corate_type_total || total_formulation_serves_T(umat_name);
    const umat_entry *entry = umat_entry_of(umat_name);
    if (entry == nullptr || entry->finite == UmatKernel::none)   // small-strain only, the external plugins included
        return false;
    return corate_type != corate_type_total || total_formulation_serves(umat_name);
}

void abaqus2smart_finite(state_variables &sv, const double *drot, const double *dfgrd0, const double *dfgrd1, const double *kin_statev, const int &corate_type, const double &DTime, const condensed_state &cs) {
    sv.F0 = mat(dfgrd0, 3, 3);
    sv.F1 = mat(dfgrd1, 3, 3);
    // F0 is the converged state of the last increment: no step cut cures a singular one
    if (det(sv.F0) <= simcoon::iota)
        throw std::invalid_argument("abaqus2smart_finite: DFGRD0 is singular; the converged state of the last "
                                    "increment cannot be cured by a step cut");
    const vec E_C_stored = sv.Etot;   // ndi < 3: the free box strain of the last increment (STATEV slots)
    const bool total = (corate_type == corate_type_total);
    if (total) {   // a state function of F
        mat lnC0;
        if (!logmat_sympd(lnC0, mat(sv.F0.t()*sv.F0)))
            throw simcoon::exception_logmat_sympd("abaqus2smart_finite: ln C0 failed");
        sv.etot = t2v_strain(0.5*lnC0);
        sv.R = eye(3, 3);
    }
    else {
        sv.etot = vec(kin_statev, 6);
        sv.R = mat(kin_statev + 6, 3, 3);
        if (!any(vectorise(sv.R)))   // STATEV starts at zero: an undeformed start
            sv.R.eye();
    }
    if (cs.active()) {   // the host's free stretches are placeholders: F0 from the strain state, F1 with the start ones
        if (total) {
            sv.etot.elem(cs.C) = E_C_stored(cs.C);
            sv.F0 = condensed_F(sv.F0).total(sv.etot);
            mat lnC1;
            if (!logmat_sympd(lnC1, mat(sv.F1.t()*sv.F1)))
                throw simcoon::exception_logmat_sympd("abaqus2smart_finite: ln C1 failed");
            vec H1 = t2v_strain(0.5*lnC1);
            H1(cs.C) = sv.etot.elem(cs.C);
            sv.F1 = condensed_F(sv.F1).total(H1);
        }
        else {
            sv.F0 = condensed_F(sv.F0)(cs, sv.etot);
            sv.F1 = condensed_F(sv.F1)(cs, sv.etot);
        }
    }
    // Abaqus rotated STRESS by DROT for its own corotational route: back to the global basis
    const mat DROT(drot, 3, 3);
    sv.sigma = det(sv.F0)*rotate_stress(sv.sigma, mat(DROT.t()));
    sv.tau = sv.sigma;
    sv.Etot = t2v_strain(Green_Lagrange(sv.F0));
    sv.DEtot = t2v_strain(Green_Lagrange(sv.F1)) - sv.Etot;
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
        const vec etot_end = finite_box_strain_end(sv, corate_type);
        const mat R_end = sv.DR*sv.R;
        std::copy(etot_end.begin(), etot_end.end(), kin_statev);
        std::copy(R_end.begin(), R_end.end(), kin_statev + 6);
    }
}

vec finite_box_strain_end(const state_variables &sv, const int &corate_type) {
    return (corate_type == corate_type_total) ? vec(sv.etot + sv.Detot) : vec(rotate_strain(sv.etot, sv.DR) + sv.Detot);
}

void smart2abaqus_finite_T(const state_variables_T &sv, const int &corate_type, const double &DTime, vec &sigma, mat &ddsdde, vec &ddsddt, double &rpl, vec &drplde, double &drpldt, double *kin_statev) {
    smart2abaqus_finite_M(sv, sv.dSdE, corate_type, sigma, ddsdde, kin_statev);
    const double J = det(sv.F1);
    const double per_volume_time = (DTime > simcoon::iota) ? 1./(J*DTime) : 0.;
    ddsddt = sv.dSdT/J;
    rpl = sv.r*per_volume_time;
    drpldt = sv.drdT*per_volume_time;
    const vec drdE = box_covector_2_dE(sv.drdE, sv.F1, corate_type);
    const vec I = {1., 1., 1., 0., 0., 0.};
    drplde = (t2v_stress(sv.F1*v2t_stress(drdE)*sv.F1.t()) - sv.r*I)*per_volume_time;
}

void smart2abaqus_M(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const vec &Wm, const mat &Lt_kernel, const bool &nlgeom, const vec &Etot_end)
{
    // the tangent arrives condensed (select_umat_M): DDSDDE is its K x K block; sigma (x) I only
    // for NLGEOM with 3 direct components, plane stress and 1D take it as it is
    const mat Lt = (nlgeom && ndi == 3) ? abaqus_jacobian(Lt_kernel, sigma) : Lt_kernel;
    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    vec2fe(sigma, K, stress);
    mat2fe(Lt, K, ddsdde);

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    pack_statev(Etot_end, statev_smart, 4, cs.C, statev);
}

void smart2abaqus_M_full(double *stress, double *ddsdde, double *stran, double *dstran, double *time, double &dtime, double &temperature, double &Dtemperature, int &nprops, double *props,  int &nstatev, double *statev, const int &ndi, const int &nshr, double *drot, const vec &sigma, const mat &Lt_kernel, const vec &Etot, const vec &DEtot, const double &T, const double &DT, const double &Time, const double &DTime, const vec &props_smart, const vec &Wm, const vec &statev_smart, const mat &DR, bool &start)
{
    UNUSED(nstatev);
    // simcoon -> the INPUT arrays of an Abaqus UMAT (UMABA plugin): ddsdde is overwritten by
    // that UMAT, so the kernel tangent goes out as it is
    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    vec2fe(sigma, K, stress);
    vec2fe(Etot, K, stran);
    vec2fe(DEtot, K, dstran);
    mat2fe(Lt_kernel, K, ddsdde);

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

void smart2abaqus_T(double *stress, double *ddsdde, double *ddsddt, double *drplde, double &drpldt, double &rpl, double *statev, const int &ndi, const int &nshr, const vec &sigma, const vec &statev_smart, const double &r, const vec &Wm, const vec &Wt, const mat &dSdE_kernel, const vec &dSdT, const vec &drpldE, const double &drpldT, const bool &nlgeom, const vec &Etot_end) {

    // dSdE, dSdT, drpldE and drpldT arrive condensed on the element's components (select_umat_T)
    const mat dSdE = (nlgeom && ndi == 3) ? abaqus_jacobian(dSdE_kernel, sigma) : dSdE_kernel;
    const condensed_state cs = make_condensed_state(ndi, nshr);
    const fea_components K = abaqus_components(cs);
    vec2fe(sigma, K, stress);
    mat2fe(dSdE, K, ddsdde);
    vec2fe(dSdT, K, ddsddt);
    vec2fe(drpldE, K, drplde);
    drpldt = drpldT;
    rpl = r;

    ///@brief : Pass the state variables
    for (int i=0; i<4; i++) {
        statev[i] = Wm(i);
    }
    for (int i=0; i<3; i++) {
        statev[i+4] = Wt(i);
    }
    pack_statev(Etot_end, statev_smart, 7, cs.C, statev);
}

// ==== ANSYS USERMAT TRANSFER FUNCTIONS ====

void ansys2smart_M(const double *stress, const double *dsdePl, const double *strain, const double *dstrain,
    const double &Time, const double &Temp, const double &dTemp,
    const int &nprops, const double *props, const int &nstatev, const double *statev,
    const int &nDirect, const int &nShear,
    vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT,
    vec &props_smart, vec &Wm, vec &statev_smart, bool &start)
{
    const fea_components K = ansys_components(nDirect, nShear);
    const condensed_state cs = make_condensed_state(nDirect, nShear);
    fe2vec(stress, K, sigma);
    fe2vec(strain, K, Etot);
    fe2vec(dstrain, K, DEtot);
    fe2mat(dsdePl, K, Lt);

    T = Temp;
    DT = dTemp;
    start = (Time < 1E-12);

    for (int i = 0; i < nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i = 0; i < 4; i++) {
        Wm(i) = statev[i];
    }
    unpack_statev(statev, nstatev, 4, cs.C, eye(3, 3), Etot, statev_smart);   // co-rotated: no rotation
}

void smart2ansys_M(double *stress, double *dsdePl, double *statev, double &sedEl, double &sedPl,
    const int &nDirect, const int &nShear,
    const vec &sigma, const mat &Lt,
    const vec &Wm, const vec &statev_smart, const vec &Etot_end)
{
    const fea_components K = ansys_components(nDirect, nShear);
    const condensed_state cs = make_condensed_state(nDirect, nShear);
    vec2fe(sigma, K, stress);
    mat2fe(Lt, K, dsdePl);

    for (int i = 0; i < 4; i++) {
        statev[i] = Wm(i);
    }
    sedEl = Wm(1);
    sedPl = Wm(2) + Wm(3);
    pack_statev(Etot_end, statev_smart, 4, cs.C, statev);
}

} // namespace simcoon
