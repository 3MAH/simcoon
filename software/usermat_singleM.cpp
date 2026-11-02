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

///@file usermat_singleM.cpp
///@brief USERMAT template to run simcoon subroutines using Ansys MAPDL
///@brief Implemented in 1D-2D-3D
///@author Chemisky
///@version 1.0
///@date 01/2026

///@details This file provides a bridge between Ansys MAPDL's USERMAT interface
///         and simcoon's constitutive models. The model is selected by the
///         first material property (props[0]) which should be set to the
///         model code (e.g., 2 for ELISO, 5 for EPICP, etc.).
///
///         Voigt notation difference:
///         - simcoon/Abaqus: (11, 22, 33, 12, 13, 23) -> indices (0,1,2,3,4,5)
///         - Ansys:          (11, 22, 33, 12, 23, 13) -> indices (0,1,2,3,4,5)
///         Indices 4 and 5 are swapped between conventions.
///
///         Key limitations of Ansys USERMAT (from MFront/TFEL documentation):
///         - External state variables are not supported
///         - Internal state variables cannot be initialized to non-zero values
///         - No standard way of defining orthotropic axes
///         - No way of controlling time step increase/decrease
///
///         Material properties in Ansys:
///         - props[0]: Model code (see select_umat_M for available models)
///         - props[1:nProp-1]: Model-specific material properties
///
///         State variables layout (same as Abaqus):
///         - ustatev[0:3]: Work quantities (Wm, Wm_r, Wm_ir, Wm_d)
///         - ustatev[4:4+nC]: stress-free strain, solved by the material point
///           (plane stress nC = 3, 1D nC = 5, otherwise none)
///         - ustatev[4+nC:nStatev-1]: Model-specific state variables
///
///         rotateM (read row by row) is the rotation increment of the state variables; an
///         array that is not a rotation reads as the identity.
///         The conversions are ansys2smart_M / smart2ansys_M (fea_transfer.hpp).

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <map>

#include <simcoon/parameter.hpp>
#include <simcoon/Simulation/Phase/phase_characteristics.hpp>
#include <simcoon/Simulation/Phase/state_variables_M.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_smart.hpp>
#include <simcoon/Continuum_mechanics/Umat/fea_transfer.hpp>

using namespace std;
using namespace arma;
using namespace simcoon;

///@brief Map from model code to 5-character UMAT name
static const std::map<int, std::string> model_code_to_name = {
    {0, "UMEXT"},   // External plugin
    {1, "UMABA"},   // Abaqus plugin
    {2, "ELISO"},   // Isotropic elasticity
    {3, "ELIST"},   // Transversely isotropic elasticity
    {4, "ELORT"},   // Orthotropic elasticity
    {5, "EPICP"},   // Isotropic plasticity (isotropic hardening)
    {6, "EPKCP"},   // Kinematic + isotropic hardening
    {7, "EPCHA"},   // Chaboche cyclic plasticity
    {8, "SMADI"},   // SMA unified model
    {9, "SMAAI"},   // SMA anisotropic model
    {10, "LLDM0"},  // Ladeveze-Le Dantec ply damage
    {11, "ZENER"},  // Zener viscoelastic (single branch)
    {12, "ZENNK"},  // Zener viscoelastic (N branches)
    {13, "PRONK"},  // Prony series viscoelastic
    {14, "EPTRI"},  // Tresca-family plasticity (Hill layout, modular adapter)
    {15, "SMADC"},  // SMA unified model, cubic elasticity
    {16, "SMAAC"},  // SMA anisotropic criterion, cubic elasticity
    {17, "EPHIL"},  // Hill anisotropic plasticity (iso hardening)
    {18, "EPHAC"},  // Hill + Chaboche
    {19, "EPANI"},  // Anisotropic plasticity
    {20, "EPDFA"},  // 
    {21, "EPCHG"},  // 
    {22, "EPHIN"},  // 
    {23, "SMAMO"},  // SMA Morin
    {24, "SMAMC"},  // SMA Morin Casciato
    {25, "SMRDI"},  // SMA unified transformation + reorientation
    {26, "SMRDC"},  // SMA unified T+R, cubic elasticity
    {27, "SMRAI"},  // SMA unified T+R, anisotropic criterion
    {28, "SMRAC"},  // SMA unified T+R, cubic + anisotropic
    {29, "EPJCK"},  // Johnson-Cook (rate and temperature dependent J2 plasticity)
    {200, "MODUL"}, // Modular composition (props encode the mechanisms)
    {100, "MIHEN"}, // Mori-Tanaka (Eshelby)
    {101, "MIMTN"}, // Mori-Tanaka N phases
    {103, "MISCN"}, // Self-consistent N phases
    {104, "MIPLN"}  // Periodic layered
};

///@brief Ansys USERMAT interface
///
///@param matId Material ID number
///@param elemId Element number
///@param kDomIntPt Domain integration point number
///@param kLayer Layer number
///@param kSectPt Section point number
///@param ldStep Load step number
///@param iSubst Substep number
///@param keycut Cutback flag (output: set to 1 when simcoon requests a smaller increment or the
///       trial state fails numerically)
///@param nDirect Number of direct stress components (1, 2, or 3)
///@param nShear Number of shear stress components (0, 1, or 3)
///@param ncomp Total stress components (nDirect + nShear)
///@param nStatev Number of state variables
///@param nProp Number of material properties (first one is model code)
///@param nPt Total number of integration points
///@param Time Current time
///@param dTime Time increment
///@param Temp Temperature at start of increment
///@param dTemp Temperature increment
///@param stress Stress array (input/output, dimension ncomp)
///@param ustatev State variables array (input/output, dimension nStatev)
///@param dsdePl Tangent modulus array (output, dimension ncomp*ncomp, column-major)
///@param sedEl Elastic strain energy density (output)
///@param sedPl Plastic strain energy density (output)
///@param epseq Equivalent plastic strain (output)
///@param strain Total strain array (input, dimension ncomp)
///@param dstrain Strain increment array (input, dimension ncomp)
///@param epsPl Plastic strain array (input, dimension ncomp)
///@param prop Material properties array (input, dimension nProp)
///@param coords Integration point coordinates (input, dimension 3)
///@param rotateM Rotation matrix (input, dimension 3x3, row-major; identity when not a rotation)
///@param defGrad_t Deformation gradient at time t (input, dimension 3x3)
///@param defGrad Deformation gradient at time t+dt (input, dimension 3x3)
///@param tsstif Transverse shear stiffness (for shells)

extern "C" void usermat_(
    int *matId, int *elemId, int *kDomIntPt, int *kLayer, int *kSectPt,
    int *ldStep, int *iSubst, int *keycut,
    int *nDirect, int *nShear, int *ncomp, int *nStatev, int *nProp,
    int *nPt, double *Time, double *dTime, double *Temp, double *dTemp,
    double *stress, double *ustatev, double *dsdePl,
    double *sedEl, double *sedPl, double *epseq,
    double *strain, double *dstrain, double *epsPl,
    double *prop, double *coords, double *rotateM,
    double *defGrad_t, double *defGrad, double *tsstif)
try {
    UNUSED(matId);
    UNUSED(elemId);
    UNUSED(kDomIntPt);
    UNUSED(kLayer);
    UNUSED(kSectPt);
    UNUSED(ldStep);
    UNUSED(iSubst);
    UNUSED(ncomp);
    UNUSED(nPt);
    UNUSED(sedEl);
    UNUSED(sedPl);
    UNUSED(epseq);
    UNUSED(epsPl);
    UNUSED(coords);
    UNUSED(defGrad_t);
    UNUSED(defGrad);
    UNUSED(tsstif);

    *keycut = 0;   // raised below when the kernel asks for a smaller increment

    bool start = false;
    int solver_type = 0;
    double pnewdt = 1.0;  // < 1 after the call: step-cut request, mapped to keycut

    // The model code is the first material property
    if (*nProp < 1)
        throw std::invalid_argument("Ansys USERMAT: the first material property must be the model code");
    const int model_code = static_cast<int>(prop[0]);
    auto it = model_code_to_name.find(model_code);
    if (it == model_code_to_name.end())
        throw std::invalid_argument("Ansys USERMAT: unknown model code " + std::to_string(model_code));
    const std::string umat_name = it->second;

    // USTATEV: Wm, the stress-free strain (nDirect < 3), then the law
    const int nstatev_smart = fea_law_statev(*nStatev, 4, *nDirect, *nShear, corate_type_total);
    const int nprops_smart = *nProp - 1;

    phase_characteristics rve;
    rve.construct(0, 1);
    auto rve_sv_M = std::dynamic_pointer_cast<state_variables_M>(rve.sptr_sv_global);
    rve_sv_M->resize(nstatev_smart);
    rve.sptr_matprops->resize(nprops_smart);
    vec props_smart = zeros(nprops_smart);
    mat DR = eye(3, 3);

    ansys2smart_M(stress, strain, dstrain, *Time, *Temp, *dTemp, nprops_smart, prop + 1, *nStatev, ustatev,
                  *nDirect, *nShear, rotateM, rve_sv_M->sigma, rve_sv_M->Etot, rve_sv_M->DEtot, rve_sv_M->T,
                  rve_sv_M->DT, props_smart, rve_sv_M->Wm, rve_sv_M->statev, DR, start);
    rve.sptr_matprops->update(0, umat_name, 1., 0., 0., 0., nprops_smart, props_smart);

    select_umat_M(rve, DR, *Time, *dTime, *nDirect, *nShear, start, solver_type, pnewdt);

    smart2ansys_M(stress, dsdePl, ustatev, *nDirect, *nShear, rve_sv_M->sigma, rve_sv_M->Lt, rve_sv_M->Wm,
                  rve_sv_M->statev, rve_sv_M->Etot + rve_sv_M->DEtot);

    // simcoon asks for a smaller increment through pnewdt < 1: Ansys bisects on keycut = 1.
    if (pnewdt < 1.) *keycut = 1;
}
catch (const std::runtime_error &e) {
    // a numerical failure of the trial state (a simcoon_error, an Armadillo decomposition): Ansys
    // retries the increment smaller. No exception may cross into the Fortran caller
    std::cerr << "simcoon USERMAT: " << e.what() << "; cutback requested\n";
    *keycut = 1;
}
catch (const std::exception &e) {
    // an input error (unknown model code, USTATEV too small, unsupported element layout): no
    // smaller increment cures it
    std::cerr << "simcoon USERMAT: " << e.what() << "; the analysis stops\n";
    std::abort();
}
