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

///@file umat_registry.cpp
///@brief The table of the constitutive law names
///@version 1.0

#include <map>
#include <string>
#include <vector>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Continuum_mechanics/Umat/umat_registry.hpp>

namespace simcoon{

namespace {

using K = UmatKernel;
using C = CondensedForms;
using T = Tensor2Type;
using SM = StressMeasure;

// a log-strain box with its tensorial statev declared (empty: no tensorial state)
umat_convention box(std::vector<StatevTensor> statev_tensors = {})
{
    return {SM::kirchhoff, true, true, std::move(statev_tensors), false};
}

// a log-strain box whose statev layout is not declared: refused under corates 4 and 5
umat_convention box_undeclared()
{
    return {SM::kirchhoff, true, false, {}, false};
}

// a kernel built on F: objective in the lab frame, native material tangent dS/dE
umat_convention on_F()
{
    return {SM::kirchhoff, false, false, {}, true};
}

// a log-strain box answering in a host code's conventions (Cauchy stress; the finite dispatch converts
// its stress and tangent), statev layout undeclared
umat_convention host_box()
{
    return {SM::cauchy, true, false, {}, false};
}

}  // namespace

// plain types only: no arma static (destructor at DLL unload, Windows)
const std::map<std::string, umat_entry> &umat_registry()
{
    static const std::map<std::string, umat_entry> registry = {
        // name            mechanical            finite                thermomechanical       finite convention                     condensed          L_eff
        // --- external plugins: the small-strain route, and the finite route as log-strain boxes in the
        //     material frame (UMEXT in the simcoon convention, as PYEXT; the host plugins in their own:
        //     Cauchy, a Jaumann (Abaqus) or Green-Naghdi (Ansys) tangent converted by the dispatch) ---
        {"UMEXT",       {K::plugin_ext,        K::plugin_ext,         K::none,                 box_undeclared()}},
        {"UMABA",       {K::plugin_abaqus,     K::plugin_abaqus,      K::none,                 host_box()}},
        {"UMANS",       {K::plugin_ansys,      K::plugin_ansys,       K::none,                 host_box()}},
        // --- linear elasticity: the modular adapter on the mechanical faces, dedicated thermomechanical kernels ---
        {"ELISO",       {K::legacy_modular,    K::legacy_modular,     K::elasticity_iso,       box(),                              C::rate_and_total, K::elasticity_iso}},
        {"ELIST",       {K::legacy_modular,    K::legacy_modular,     K::elasticity_trans_iso, box(),                              C::rate_and_total, K::elasticity_trans_iso}},
        {"ELORT",       {K::legacy_modular,    K::legacy_modular,     K::elasticity_ortho,     box(),                              C::rate_and_total, K::elasticity_ortho}},
        // --- plasticity: the modular adapter; the dedicated kernels under _closed ---
        {"EPICP",       {K::legacy_modular,    K::legacy_modular,     K::legacy_modular,       box({{2, T::strain}}),              C::rate_and_total}},   // [T_init, p, EP]
        {"EPICP_closed",{K::plasticity_iso,    K::plasticity_iso,     K::plasticity_iso,       box({{2, T::strain}})}},
        {"EPKCP",       {K::legacy_modular,    K::legacy_modular,     K::legacy_modular,       box_undeclared(),                   C::rate_and_total}},
        {"EPKCP_closed",{K::none,              K::none,               K::plasticity_kin_iso}},
        {"EPCHA",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box({{2, T::strain},                // EP
                                                                                                   {8, T::strain},                // a_1 (no stored back stress)
                                                                                                   {14, T::strain}}),             // a_2
                                                                                                                                  C::rate_and_total}},
        {"EPCHA_closed",{K::plasticity_chaboche, K::plasticity_chaboche, K::none,              box({{2, T::strain},                // EP
                                                                                                   {8, T::strain},                // a_1
                                                                                                   {14, T::strain},               // a_2
                                                                                                   {20, T::stress},               // X_1
                                                                                                   {26, T::stress}})}},           // X_2
        {"EPJCK",       {K::johnson_cook,      K::johnson_cook,       K::johnson_cook,         box({{2, T::strain}}),              C::rate_and_total}},   // EP, the EPICP_closed layout
        {"EPHIL",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared(),                   C::rate_and_total}},
        {"EPTRI",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        {"EPHAC",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        {"EPANI",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        {"EPDFA",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        {"EPCHG",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        {"EPHIN",       {K::legacy_modular,    K::legacy_modular,     K::none,                 box_undeclared()}},
        // --- shape memory alloys and damage: small-strain route ---
        {"SMADI",       {K::sma_unified,       K::none,               K::sma_unified,          {},                                 C::rate}},
        {"SMADC",       {K::sma_unified,       K::none,               K::sma_unified}},
        {"SMAAI",       {K::sma_unified,       K::none,               K::sma_unified}},
        {"SMAAC",       {K::sma_unified,       K::none,               K::sma_unified}},
        {"SMRDI",       {K::sma_unified_R,     K::none,               K::none}},
        {"SMRDC",       {K::sma_unified_R,     K::none,               K::none}},
        {"SMRAI",       {K::sma_unified_R,     K::none,               K::none}},
        {"SMRAC",       {K::sma_unified_R,     K::none,               K::none}},
        {"SMAMO",       {K::sma_mono,          K::none,               K::none}},
        {"SMAMC",       {K::sma_mono,          K::none,               K::none}},
        {"LLDM0",       {K::damage_LLD0,       K::none,               K::none}},
        // --- viscoelasticity ---
        {"ZENER",       {K::legacy_modular,    K::legacy_modular,     K::legacy_modular,       box_undeclared(),                   C::rate}},   // one Maxwell branch by parameter map
        {"ZENER_closed",{K::zener,             K::none,               K::zener}},
        {"ZENNK",       {K::zener_N,           K::none,               K::zener_N,              {},                                 C::rate}},
        {"PRONK",       {K::legacy_modular,    K::legacy_modular,     K::legacy_modular,       box_undeclared(),                   C::rate}},   // one block per branch: no declared layout
        {"PRONK_closed",{K::prony_N,           K::none,               K::prony_N}},
        // --- finite-strain kernels ---
        {"HYPOO",       {K::none,              K::hypoelasticity_ortho, K::none,               box()}},   // the rate counterpart of ELORT, tau_{n+1} = tau_n + L : DEel
        {"SNTVE",       {K::none,              K::saint_venant,       K::none,                 on_F(),                             C::rate_and_total}},
        {"NEOHI",       {K::none,              K::neo_hookean_incomp, K::none,                 on_F(),                             C::rate_and_total}},
        {"NEOHC",       {K::none,              K::hyper_invariants,   K::none,                 on_F(),                             C::rate_and_total}},
        {"MOORI",       {K::none,              K::hyper_invariants,   K::none,                 on_F(),                             C::rate_and_total}},
        {"YEOHH",       {K::none,              K::hyper_invariants,   K::none,                 on_F(),                             C::rate_and_total}},
        {"ISHAH",       {K::none,              K::hyper_invariants,   K::none,                 on_F()}},
        {"GETHH",       {K::none,              K::hyper_invariants,   K::none,                 on_F()}},
        {"SWANH",       {K::none,              K::hyper_invariants,   K::none,                 on_F()}},
        {"HOLZA",       {K::none,              K::hyper_invariants,   K::none,                 on_F()}},
        {"OGDEN",       {K::none,              K::hyper_pstretch,     K::none,                 on_F()}},
        // --- the modular engine: a Hencky hyperelastic composition under finite strain (corate 3 or the total form) ---
        {"MODUL",       {K::modular,           K::modular,            K::modular,              box_undeclared(),                   C::rate_and_total}},
        // --- a Python law through the registered callback: fed the log strain, returns tau ---
        {"PYEXT",       {K::python_callback,   K::python_callback,    K::none,                 box_undeclared()}},
        // --- mean-field schemes over sub-phases ---
        {"MIHEN",       {K::mean_field_homogeneous,     K::none,      K::none,                 {},                                 C::none,           K::mean_field_homogeneous}},
        {"MIMTN",       {K::mean_field_mori_tanaka,     K::none,      K::none,                 {},                                 C::none,           K::mean_field_mori_tanaka}},
        {"MISCN",       {K::mean_field_self_consistent, K::none,      K::none,                 {},                                 C::none,           K::mean_field_self_consistent}},
        {"MIPLN",       {K::mean_field_periodic_layer,  K::none,      K::none,                 {},                                 C::none,           K::mean_field_periodic_layer}},
    };
    return registry;
}

const umat_entry *umat_entry_of(const std::string &umat_name)
{
    const auto &registry = umat_registry();
    const auto it = registry.find(umat_name);
    return (it != registry.end()) ? &it->second : nullptr;
}

int mean_field_scheme(const UmatKernel &kernel)
{
    switch (kernel) {
        case K::mean_field_homogeneous: return 100;
        case K::mean_field_mori_tanaka: return 101;
        case K::mean_field_self_consistent: return 103;
        case K::mean_field_periodic_layer: return 104;
        default: return 0;
    }
}

bool kernel_reads_ncp(const UmatKernel &kernel, const bool &thermomechanical)
{
    switch (kernel) {
        case K::modular: case K::legacy_modular: return true;
        case K::plasticity_iso: return !thermomechanical;
        default: return false;
    }
}

bool kernel_serves_closest_point(const UmatKernel &kernel, const bool &thermomechanical)
{
    switch (kernel) {
        case K::modular: case K::legacy_modular: case K::plasticity_iso: return true;
        case K::plasticity_kin_iso: return thermomechanical;
        case K::plasticity_chaboche: case K::johnson_cook: case K::sma_unified: case K::sma_unified_R: return !thermomechanical;
        default: return false;
    }
}

} //namespace simcoon
