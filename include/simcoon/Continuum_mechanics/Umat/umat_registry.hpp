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

/**
 * @file umat_registry.hpp
 * @brief The one table of the constitutive law names: for each @c umat_name, the kernel that
 *        serves it on each face of the dispatch, its finite-strain conventions and its
 *        capabilities.
 *
 * Every dispatch of the library reads this table: the small-strain and finite-strain
 * mechanical routes and the thermomechanical one (umat_smart.cpp), the point-wise entries
 * @c sim.umat / @c sim.umat_T of the Python module, the Abaqus and Ansys bridges
 * (fea_transfer.cpp), the elastic stiffness of the mean-field schemes (umat_L_elastic.cpp) and
 * the sub-phase layout of a mean-field name (multiphase.cpp). Adding a law is one row here and
 * one @c case per face that serves it; a row without its case is a thrown @c std::logic_error,
 * never a silent default.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace simcoon{

/**
 * @brief What a kernel's raw @c sigma out-parameter holds.
 *
 * Every NATIVE simcoon kernel is @c kirchhoff: in the logarithmic framework the potential
 * differentiates to \f$ \boldsymbol{\tau} \f$ per REFERENCE volume. It is the default route
 * stress, and \f$ \boldsymbol{\sigma} = \boldsymbol{\tau}/J \f$ is formed at the output
 * boundaries only. @c cauchy is the contract of a kernel whose convention belongs to a host code
 * (the Abaqus and Ansys plugins): the finite dispatch hands it \f$ \sigma_n = \tau_n / J_n \f$ and
 * converts its stress back (\f$ \tau = J \sigma \f$) and its tangent (the host's rate) to the route.
 */
enum class StressMeasure { kirchhoff, cauchy };

// Declared in Functions/tensor.hpp (kept out of this widely included header).
enum class Tensor2Type;

/**
 * @brief One tensorial internal variable of a kernel's @c statev: 6 engineering Voigt
 *        components starting at @c offset, with its Tensor2Type (@c strain: engineering shear;
 *        @c stress). The type selects the transport, through tensor2 (see umat_convention).
 */
struct StatevTensor {
    int offset;
    Tensor2Type type;
};

/**
 * @brief The conventions a kernel's raw outputs are expressed in on the finite-strain route.
 *
 * The tangent rate is deliberately absent: every kernel is handed the solver's
 * @c corate_type and must emit \f$ \mathbf{L}_t \f$ in it. Unlike the stress measure, the
 * rate is therefore a public contract, not a per-kernel declaration.
 *
 * @c material_frame marks a kernel fed the logarithmic strain (a "box" kernel): on the finite
 * route it runs in the frame that follows the material, \f$ \hat{\mathbf{R}}_{n+1} \f$ (the
 * accumulated corate rotation for corates 0-3, the polar rotation of \f$ \mathbf{F}_1 \f$ for 4
 * and 5), so its anisotropy axes rotate with the body and its @c statev is stored in that
 * frame. It is then called with \f$ \Delta\mathbf{R} = \mathbf{I} \f$. Kernels built from
 * \f$ \mathbf{F} \f$ are objective by construction and keep the lab frame.
 *
 * For corates 4 and 5 the frame-relative increment is a rotation-free stretch, and the
 * dispatcher transports the tensors listed in @c statev_tensors with it, each with its variance.
 * @c layout_declared says the list is complete (an empty list: no tensorial state); a material
 * frame kernel that has not declared it is refused under corates 4 and 5.
 *
 * @c material_tangent marks a kernel built on \f$ \mathbf{F} \f$ whose native tangent is the
 * material one, \f$ \partial\mathbf{S}/\partial\mathbf{E} \f$ (a pull-back of its spatial tangent):
 * called with corate_type_material_tangent it returns that tangent instead of a box, and its
 * stress and work do not depend on the corate it is given. A caller that needs
 * \f$ \partial\mathbf{S}/\partial\mathbf{E} \f$ asks for it (select_umat_M_finite) instead of
 * converting the box back.
 *
 * @see output_convention_of
 */
struct umat_convention {
    StressMeasure stress = StressMeasure::kirchhoff;
    bool material_frame = false;
    bool layout_declared = false;
    std::vector<StatevTensor> statev_tensors = {};
    bool material_tangent = false;
};

/**
 * @brief The kernels of the library: the function a face of the dispatch calls for a name.
 *
 * One value names one implementation, whatever the face: @c plasticity_iso is
 * umat_plasticity_iso() on the mechanical faces and umat_plasticity_iso_T() on the
 * thermomechanical one. The dedicated kernels of the names the modular engine serves stay
 * reachable under their @c _closed name (EPICP_closed is @c plasticity_iso, EPICP is
 * @c legacy_modular). The mean-field schemes are kernels of the small-strain face only, and
 * mean_field_scheme() gives the code umat_multi() switches on.
 */
enum class UmatKernel {
    none,                       ///< the face does not serve the name
    plugin_ext,                 ///< UMEXT: umat_plugin_ext_api, loaded at run time
    plugin_abaqus,              ///< UMABA: an Abaqus UMAT behind umat_plugin_aba_api
    plugin_ansys,               ///< UMANS: an Ansys USERMAT behind umat_plugin_ans_api
    elasticity_iso,             ///< umat_elasticity_iso_T; the closed-form stiffness L_iso of L_eff
    elasticity_trans_iso,       ///< umat_elasticity_trans_iso_T; L_isotrans
    elasticity_ortho,           ///< umat_elasticity_ortho_T; L_ortho
    plasticity_iso,             ///< umat_plasticity_iso / _T (EPICP_closed)
    plasticity_kin_iso,         ///< umat_plasticity_kin_iso_T (EPKCP_closed)
    plasticity_chaboche,        ///< umat_plasticity_chaboche (EPCHA_closed)
    johnson_cook,               ///< umat_plasticity_johnson_cook_CCP / _T (EPJCK)
    damage_LLD0,                ///< umat_damage_LLD_0 (LLDM0)
    zener,                      ///< umat_zener_fast / _T (ZENER_closed)
    zener_N,                    ///< umat_zener_Nfast / _T (ZENNK)
    prony_N,                    ///< umat_prony_Nfast / _T (PRONK_closed)
    sma_unified,                ///< umat_sma_unified_T / _T_T (SMADI, SMADC, SMAAI, SMAAC)
    sma_unified_R,              ///< umat_sma_unified_TR (SMRDI, SMRDC, SMRAI, SMRAC)
    sma_mono,                   ///< umat_sma_mono (SMAMO, SMAMC)
    hypoelasticity_ortho,       ///< umat_hypoelasticity_ortho (HYPOO), finite route
    saint_venant,               ///< umat_saint_venant (SNTVE), finite route
    neo_hookean_incomp,         ///< umat_neo_hookean_incomp (NEOHI), finite route
    hyper_invariants,           ///< umat_generic_hyper_invariants (NEOHC, MOORI, YEOHH, ISHAH, GETHH, SWANH, HOLZA)
    hyper_pstretch,             ///< umat_generic_hyper_pstretch (OGDEN)
    mean_field_homogeneous,     ///< umat_multi, scheme 100 (MIHEN)
    mean_field_mori_tanaka,     ///< umat_multi, scheme 101 (MIMTN)
    mean_field_self_consistent, ///< umat_multi, scheme 103 (MISCN)
    mean_field_periodic_layer,  ///< umat_multi, scheme 104 (MIPLN)
    modular,                    ///< umat_modular / _T (MODUL)
    legacy_modular,             ///< umat_legacy_modular / _T: a legacy name translated to the modular engine (legacy_adapters.cpp)
    python_callback             ///< umat_callback_M (PYEXT)
};

/**
 * @brief The condensed stress states (ndi < 3) @c sim.umat has validated for a name against the
 *        3D law at the stress-free strain (test_umat_condensed.py): none, the rate form under
 *        log_R, or the rate and the total forms.
 */
enum class CondensedForms { none, rate, rate_and_total };

/**
 * @brief One row of the registry: what each face of the dispatch calls for a name, and what the
 *        name can do.
 *
 * @c mechanical and @c finite are the small-strain and finite-strain mechanical faces
 * (select_umat_M, select_umat_M_finite): a name on both faces has the same kernel on both, fed
 * the small strain or the logarithmic box strain (the Python @c sim.umat relies on it).
 * @c thermomechanical is the face of select_umat_T on both routes (the finite one admits every
 * thermomechanical name under log_R and the total form). @c convention is read on the finite
 * face only. @c elastic_stiffness is the closed-form stiffness get_L_elastic() builds for the
 * name (L_eff of the micromechanics): the elasticity kernels from the props, a mean-field kernel
 * from the sub-phases, @c none for a law L_eff refuses.
 */
struct umat_entry {
    UmatKernel mechanical = UmatKernel::none;
    UmatKernel finite = UmatKernel::none;
    UmatKernel thermomechanical = UmatKernel::none;
    umat_convention convention = {};
    CondensedForms condensed = CondensedForms::none;
    UmatKernel elastic_stiffness = UmatKernel::none;
};

/**
 * @brief The registry: every name the library serves, ordered by name.
 */
const std::map<std::string, umat_entry> &umat_registry();

/**
 * @brief The row of a name, or @c nullptr for a name the library does not serve. Never throws.
 */
const umat_entry *umat_entry_of(const std::string &umat_name);

/**
 * @brief The mean-field scheme code umat_multi() and get_L_elastic() switch on (100
 *        homogeneous, 101 Mori-Tanaka, 103 self-consistent, 104 periodic layers), 0 for any
 *        other kernel.
 */
int mean_field_scheme(const UmatKernel &kernel);

/**
 * @brief True when a kernel reads the complementarity option of a run (NCPMethod) in its
 *        cutting plane: the modular engine, the names it serves and the dedicated
 *        plasticity_iso on the mechanical faces; the modular engine and its adapters on the
 *        thermomechanical one.
 */
bool kernel_reads_ncp(const UmatKernel &kernel, const bool &thermomechanical);

/**
 * @brief True when a kernel has a closest-point update (integrator_closest_point): on the
 *        mechanical faces plasticity_iso, plasticity_chaboche, johnson_cook (its radial cutting
 *        plane is the projection by construction), the modular engine and its adapters, the SMA
 *        kernels sma_unified and sma_unified_R; on the thermomechanical face plasticity_iso,
 *        plasticity_kin_iso, the modular engine and its adapters.
 */
bool kernel_serves_closest_point(const UmatKernel &kernel, const bool &thermomechanical);

} //namespace simcoon
