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
 * @file legacy_adapters.cpp
 * @brief Implementation of the legacy-name props translators
 */

#include <simcoon/Continuum_mechanics/Umat/Modular/legacy_adapters.hpp>
#include <simcoon/Continuum_mechanics/Umat/Modular/modular_umat.hpp>
#include <cstring>
#include <map>
#include <stdexcept>
#include <vector>

using namespace arma;

namespace simcoon {

namespace {

// MODUL props stream: see ModularUMAT::configure_from_props (modular_umat.hpp)

vec translate_ELISO(const vec& p) {
    // legacy: E nu alpha
    return vec{0, 0, p(0), p(1), p(2), 0};
}

vec translate_ELIST(const vec& p) {
    // legacy: axis EL ET nuTL nuTT GLT alpha_L alpha_T  (axis FIRST)
    return vec{2, 0, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(0), 0};
}

vec translate_ELORT(const vec& p) {
    // legacy: E1 E2 E3 nu12 nu13 nu23 G12 G13 G23 alpha1 alpha2 alpha3
    return vec{3, 0, p(0), p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8),
               p(9), p(10), p(11), 0};
}

vec translate_EPICP(const vec& p) {
    // legacy: E nu alpha sigmaY k m  ->  plasticity, Mises, PowerLaw, no kin
    return vec{0, 0, p(0), p(1), p(2),
               1,
               0, 0, 2, 0, 1, 0,
               p(3), p(4), p(5)};
}

vec translate_EPCHA(const vec& p) {
    // legacy: E nu alpha sigmaY Q b C1 D1 C2 D2  ->  plasticity, Mises, Voce, 2 Armstrong-Frederick terms
    return vec{0, 0, p(0), p(1), p(2),
               1,
               0, 0, 3, 3, 1, 2,
               p(3), p(4), p(5), p(6), p(7), p(8), p(9)};
}

vec translate_PRONK(const vec& p) {
    // legacy: E0 nu0 alpha N (E_i nu_i etaB_i etaS_i) x N  ->  one viscoelastic mechanism, the branches verbatim
    const int N = static_cast<int>(p(3));
    vec out(8 + 4 * static_cast<uword>(N));
    out(0) = 0; out(1) = 0;
    out(2) = p(0); out(3) = p(1); out(4) = p(2);
    out(5) = 1;
    out(6) = 1; out(7) = N;
    for (uword i = 0; i < 4 * static_cast<uword>(N); ++i) out(8 + i) = p(4 + i);
    return out;
}

vec translate_EPKCP(const vec& p) {
    // legacy: E nu alpha sigmaY k m kX; legacy Prager X = kX*a (tensorial),
    // modular X = (2/3) C a  ->  C = 1.5 kX
    return vec{0, 0, p(0), p(1), p(2),
               1,
               0, 0, 2, 1, 1, 1,          // plasticity, Mises, PowerLaw, Prager
               p(3), p(4), p(5), 1.5 * p(6)};
}

vec translate_EPHIL(const vec& p) {
    // legacy: E nu alpha | sigmaY k m | F G H L M N  (EPTRI shares this)
    return vec{0, 0, p(0), p(1), p(2),
               1,
               0, 3, 2, 0, 1, 0,          // plasticity, Hill, PowerLaw, no kin
               p(3),
               p(6), p(7), p(8), p(9), p(10), p(11),
               p(4), p(5)};
}

// Chaboche-family head: E nu G alpha | sigmaY Q b | C1 D1 C2 D2 | criterion...
// (cubic elasticity, Voce iso, 2 AF backstresses)
vec chaboche_family(const vec& p, int yield_type, int n_crit) {
    // 14 head slots (elasticity 6 + n_mech + mech type + 5 config + sigma_Y)
    // + criterion params + Voce (2) + two AF terms (4)
    vec out(20 + static_cast<uword>(n_crit));
    out(0) = 1; out(1) = 0;                       // cubic, EnuG
    out(2) = p(0); out(3) = p(1); out(4) = p(2); out(5) = p(3);
    out(6) = 1;                                   // n_mech
    out(7) = 0;                                   // plasticity
    out(8) = yield_type; out(9) = 3; out(10) = 3; // yield, Voce, Chaboche
    out(11) = 1; out(12) = 2;                     // N_iso, N_kin
    out(13) = p(4);                               // sigma_Y
    for (int i = 0; i < n_crit; ++i) {
        out(14 + static_cast<uword>(i)) = p(11 + static_cast<uword>(i));
    }
    const uword o = 14 + static_cast<uword>(n_crit);
    out(o) = p(5); out(o + 1) = p(6);                       // Q, b
    out(o + 2) = p(7); out(o + 3) = p(8);                   // C1, D1
    out(o + 4) = p(9); out(o + 5) = p(10);                  // C2, D2
    return out;
}

vec translate_EPHAC(const vec& p) { return chaboche_family(p, 3, 6); }
vec translate_EPANI(const vec& p) { return chaboche_family(p, 5, 9); }
vec translate_EPDFA(const vec& p) { return chaboche_family(p, 4, 7); }

vec translate_EPCHG(const vec& p) {
    // legacy: E nu G alpha | sigmaY N_iso N_kin criteria | (Q,b)xN_iso |
    //         (C,D)xN_kin | criterion params
    const int N_iso = static_cast<int>(p(5));
    const int N_kin = static_cast<int>(p(6));
    const int crit_legacy = static_cast<int>(p(7));
    // legacy criteria codes: 0 = Mises, 1 = Hill, 2 = DFA, 3 = anisotropic
    static const int yield_map[4] = {0, 3, 4, 5};
    static const int n_crit_map[4] = {0, 6, 7, 9};
    if (crit_legacy < 0 || crit_legacy > 3) {
        throw std::invalid_argument("EPCHG: unknown criteria code "
                                 + std::to_string(crit_legacy));
    }
    const int n_crit = n_crit_map[crit_legacy];
    const uword legacy_crit_off = 8 + 2 * static_cast<uword>(N_iso + N_kin);

    // legacy N-term hardening shares ONE Hp: an effective Voce, b_eff = sum b_i, Q_eff = sum(b_i Q_i)/b_eff
    double b_eff = 0.0, bq = 0.0;
    for (int i = 0; i < N_iso; ++i) {
        const double Q_i = p(8 + 2 * static_cast<uword>(i));
        const double b_i = p(9 + 2 * static_cast<uword>(i));
        b_eff += b_i;
        bq += b_i * Q_i;
    }
    const double Q_eff = (b_eff > 0.0) ? bq / b_eff : 0.0;

    // 14 head slots + criterion params + effective Voce (2) + 2*N_kin
    vec out(16 + 2 * static_cast<uword>(N_kin) + static_cast<uword>(n_crit));
    out(0) = 1; out(1) = 0;
    out(2) = p(0); out(3) = p(1); out(4) = p(2); out(5) = p(3);
    out(6) = 1;
    out(7) = 0;
    out(8) = yield_map[crit_legacy];
    out(9) = 3; out(10) = 3;                       // Voce (effective), Chaboche
    out(11) = 1; out(12) = N_kin;
    out(13) = p(4);
    uword o = 14;
    for (int i = 0; i < n_crit; ++i) {
        out(o++) = p(legacy_crit_off + static_cast<uword>(i));
    }
    out(o++) = Q_eff;
    out(o++) = b_eff;
    for (int i = 0; i < 2 * N_kin; ++i) {
        out(o++) = p(8 + 2 * static_cast<uword>(N_iso) + static_cast<uword>(i));
    }
    return out;
}

vec translate_EPHIN(const vec& p) {
    // legacy: E nu alpha N | (sigmaY k m F G H L M N)xN  ->  N Hill mechanisms
    const int N = static_cast<int>(p(3));
    vec out(5 + 1 + static_cast<uword>(N) * 15);
    out(0) = 0; out(1) = 0;
    out(2) = p(0); out(3) = p(1); out(4) = p(2);
    out(5) = N;
    uword o = 6;
    for (int i = 0; i < N; ++i) {
        const uword b = 4 + static_cast<uword>(i) * 9;
        out(o++) = 0;                                        // plasticity
        out(o++) = 3; out(o++) = 2; out(o++) = 0;            // Hill, PowerLaw, no kin
        out(o++) = 1; out(o++) = 0;                          // N_iso, N_kin
        out(o++) = p(b);                                     // sigma_Y
        for (int j = 3; j < 9; ++j) {                        // F G H L M N
            out(o++) = p(b + static_cast<uword>(j));
        }
        out(o++) = p(b + 1); out(o++) = p(b + 2);            // k, m
    }
    return out;
}


const std::map<std::string, LegacyPropsTranslator>& registry() {
    static const std::map<std::string, LegacyPropsTranslator> reg = {
        {"ELISO", translate_ELISO},
        {"ELIST", translate_ELIST},
        {"ELORT", translate_ELORT},
        {"EPICP", translate_EPICP},
        {"EPCHA", translate_EPCHA},
        {"PRONK", translate_PRONK},
        {"EPKCP", translate_EPKCP},
        {"EPHIL", translate_EPHIL},
        {"EPTRI", translate_EPHIL},
        {"EPHAC", translate_EPHAC},
        {"EPANI", translate_EPANI},
        {"EPDFA", translate_EPDFA},
        {"EPCHG", translate_EPCHG},
        {"EPHIN", translate_EPHIN},
    };
    return reg;
}

}  // namespace

bool has_legacy_adapter(const std::string& umat_name) {
    return registry().count(umat_name) > 0;
}

namespace {
// per-thread memo of the last (name, props) translation of each face (plain std containers: no arma
// object outlives the call); the key is written last so a throw leaves no stale pair. The
// thermomechanical props carry [rho, c_p] in front of the mechanical legacy props.
vec translated_props(const std::string& umat_name, const vec& props, const bool thermo) {
    struct translation {
        std::string name;
        std::vector<double> props;
        std::vector<double> props_mod;
    };
    thread_local translation memo[2];
    translation& last = memo[thermo ? 1 : 0];
    const bool hit = !last.name.empty() && umat_name == last.name && props.n_elem == last.props.size()
        && (last.props.empty()
            || std::memcmp(last.props.data(), props.memptr(), last.props.size() * sizeof(double)) == 0);
    if (!hit) {
        const auto it = registry().find(umat_name);
        if (it == registry().end()) {
            throw std::invalid_argument(std::string(thermo ? "umat_legacy_modular_T" : "umat_legacy_modular")
                                     + ": no adapter registered for '" + umat_name + "'");
        }
        vec translated;
        if (thermo) {
            if (props.n_elem < 2) {
                throw std::invalid_argument("umat_legacy_modular_T: '" + umat_name + "' needs [rho, c_p, ...]");
            }
            translated = join_cols(props.subvec(0, 1), it->second(props.subvec(2, props.n_elem - 1)));
        } else {
            translated = it->second(props);
        }
        last.name.clear();
        last.props.assign(props.begin(), props.end());
        last.props_mod.assign(translated.begin(), translated.end());
        last.name = umat_name;
    }
    return vec(const_cast<double*>(last.props_mod.data()), last.props_mod.size(), false, true);   // a view
}
}  // namespace

void umat_legacy_modular(
    const std::string& umat_name,
    const arma::vec& Etot, const arma::vec& DEtot,
    arma::vec& sigma, arma::mat& Lt, arma::mat& L, const arma::mat& DR,
    const int& /*nprops*/, const arma::vec& props,
    const int& nstatev, arma::vec& statev,
    const double& T, const double& DT, const double& Time, const double& DTime,
    double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
    const int& ndi, const int& nshr, const bool& start, double& tnew_dt,
    const int& tangent_mode, const int& ncp) {

    const vec props_mod = translated_props(umat_name, props, false);   // a view of the memo
    try {
        umat_modular(umat_name, Etot, DEtot, sigma, Lt, L, DR,
                     static_cast<int>(props_mod.n_elem), props_mod,
                     nstatev, statev,
                     T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d,
                     ndi, nshr, start, tnew_dt, tangent_mode, KinematicForm::rate, ncp);
    } catch (const std::logic_error& e) {
        // a configuration error, named after the legacy law; numerical failures pass through
        // with their own type (the step-cut policy)
        throw std::invalid_argument("'" + umat_name + "' (modular adapter): "
                                 + e.what());
    }
}

void umat_legacy_modular_T(
    const std::string& umat_name,
    const arma::vec& Etot, const arma::vec& DEtot,
    arma::vec& sigma, double& r, arma::mat& dSdE, arma::vec& dSdT, arma::vec& drdE, double& drdT,
    const arma::mat& DR,
    const int& /*nprops*/, const arma::vec& props,
    const int& nstatev, arma::vec& statev,
    const double& T, const double& DT, const double& Time, const double& DTime,
    double& Wm, double& Wm_r, double& Wm_ir, double& Wm_d,
    double& Wt, double& Wt_r, double& Wt_ir,
    const int& ndi, const int& nshr, const bool& start, double& tnew_dt,
    const int& tangent_mode, const KinematicForm& kinematic_form, const int& ncp) {

    const vec props_mod = translated_props(umat_name, props, true);   // a view of the memo
    try {
        umat_modular_T(Etot, DEtot, sigma, r, dSdE, dSdT, drdE, drdT, DR,
                       static_cast<int>(props_mod.n_elem), props_mod, nstatev, statev,
                       T, DT, Time, DTime, Wm, Wm_r, Wm_ir, Wm_d, Wt, Wt_r, Wt_ir,
                       ndi, nshr, start, tnew_dt, tangent_mode, kinematic_form, ncp);
    } catch (const std::logic_error& e) {
        throw std::invalid_argument("'" + umat_name + "' (modular adapter, thermomechanical): "
                                 + e.what());
    }
}

} // namespace simcoon
