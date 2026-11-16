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

///@file criteria.cpp
///@brief Provide function for yield surfaces
///@version 1.0

#include <iostream>
#include <stdexcept>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/exception.hpp>
#include <simcoon/Continuum_mechanics/Functions/constitutive.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/transfer.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Continuum_mechanics/Functions/tensor.hpp>
#include <simcoon/Simulation/Maths/small_dense.hpp>
#include <simcoon/Simulation/Maths/small_voigt.hpp>

using namespace std;
using namespace arma;

namespace simcoon
{

    namespace {

    using vec6 = vec::fixed<6>;
    using mat6 = mat::fixed<6, 6>;
    using mat3 = mat::fixed<3, 3>;

    // the response of a stress below the floor of a criterion: zero up to the level asked
    CriterionResponse zero_response(const int &level)
    {
        CriterionResponse out;
        if (level >= 1) out.df.zeros();
        if (level >= 2) out.ddf.zeros();
        return out;
    }

    // deviator of a 3x3
    mat3 dev3(const mat3 &m)
    {
        return m - (trace(m)/3.)*mat3(fill::eye);
    }

    // strain-like Voigt of a symmetric 3x3 (doubled shears)
    vec6 strain_voigt(const mat3 &m)
    {
        vec6 x;
        voigt_from_mat(x.memptr(), m, voigt_shear_out(Tensor2Type::strain));
        return x;
    }

    // Drucker J3 correction applied to a base response (Mises or DFA)
    CriterionResponse drucker_correction(const vec &v, const CriterionResponse &base, const double &b, const double &n,
                                         const int &level)
    {
        if (!(base.f > simcoon::iota)) return zero_response(level);
        if (n >= 10.) return base;
        const mat3 s = mat_from_voigt(dev(v).memptr(), voigt_shear_in(Tensor2Type::stress));
        const double J2 = 0.5*accu(s % s);
        if (!(sqrt(3.*J2) > simcoon::iota)) return base; // DFA base: pressure part without deviator

        const double m = 1./n;
        const double J3 = det(s);
        const double J2_32 = J2*sqrt(J2);
        const double r = 1. + b*J3/J2_32;
        const double r_m = pow(r, m);
        CriterionResponse out;
        out.f = base.f*r_m;
        if (level < 1) return out;

        // gradients of J2 and J3 (strain-like) and of J3/J2^(3/2)
        const double J2_52 = J2_32*J2;
        const mat3 s2 = s*s;
        const vec6 dJ2 = strain_voigt(s);
        const vec6 dJ3 = strain_voigt(dev3(s2));
        const vec6 q = dJ3/J2_32 - (1.5*J3/J2_52)*dJ2;
        out.df = r_m*base.df + (base.f*m*b*r_m/r)*q;
        if (level < 2) return out;

        // d(dev(s s))/dsigma per Voigt component, a shear component moving both symmetric entries
        mat6 ddJ3;
        for (int c = 0; c < 6; ++c) {
            mat3 E(fill::zeros);
            E(voigt_i[c], voigt_j[c]) = 1.;
            E(voigt_j[c], voigt_i[c]) = 1.;
            const mat3 dS = dev3(E);
            ddJ3.col(c) = strain_voigt(dev3(dS*s + s*dS));
        }
        const mat6 ddJ2 = Idev2();
        const double J2_72 = J2_52*J2;
        const mat6 Q = ddJ3/J2_32 - (1.5/J2_52)*(outer6(dJ3.memptr(), dJ2.memptr()) + outer6(dJ2.memptr(), dJ3.memptr()))
                       - (1.5*J3/J2_52)*ddJ2 + (3.75*J3/J2_72)*outer6(dJ2.memptr(), dJ2.memptr());
        out.ddf = r_m*base.ddf
                + (m*b*r_m/r)*(outer6(base.df.memptr(), q.memptr()) + outer6(q.memptr(), base.df.memptr()) + base.f*Q)
                + (base.f*m*(m - 1.)*b*b*r_m/(r*r))*outer6(q.memptr(), q.memptr());
        return out;
    }

    // anisotropy operators into a caller-provided 6x6 (stack in the criteria, heap in P_Ani/P_Hill/P_DFA)
    void fill_P_Ani(const vec &params, mat &P)
    {
        assert(params.n_elem == 9); // P_11,P_22,P_33,P_12,P_13,P_23,P_44,P_55,P_66
        P.zeros();
        P(0, 0) = params(0);      // P_11
        P(1, 1) = params(1);      // P_22
        P(2, 2) = params(2);      // P_33
        P(0, 1) = params(3);      // P_12
        P(1, 0) = params(3);      // P_12
        P(0, 2) = params(4);      // P_13
        P(2, 0) = params(4);      // P_13
        P(1, 2) = params(5);      // P_23
        P(2, 1) = params(5);      // P_23
        P(3, 3) = 2. * params(6); // P_44 = P_1212
        P(4, 4) = 2. * params(7); // P_55 = P_1313
        P(5, 5) = 2. * params(8); // P_66 = P_2323
    }

    void fill_P_Hill(const vec &params, mat &P)
    {
        assert(params.n_elem == 6); // F,G,H,L,M,N
        P.zeros();
        // params = (F, G, H, L, M, N)
        P(0, 0) = params(1) + params(2); // P_11 = G+H
        P(1, 1) = params(0) + params(2); // P_22 = F+H
        P(2, 2) = params(0) + params(1); // P_33 = F+G
        P(0, 1) = -1. * params(2);       // P_12 = -H
        P(1, 0) = -1. * params(2);       // P_12 = -H
        P(0, 2) = -1. * params(1);       // P_13 = -G
        P(2, 0) = -1. * params(1);       // P_13 = -G
        P(1, 2) = -1. * params(0);       // P_23 = -F
        P(2, 1) = -1. * params(0);       // P_23 = -F
        P(3, 3) = 2. * params(5);        // P_44 = N
        P(4, 4) = 2. * params(4);        // P_55 = M
        P(5, 5) = 2. * params(3);        // P_66 = L
    }

    void fill_P_DFA(const vec &params, mat &P)
    {                               // Deshpande–Fleck–Ashby
        assert(params.n_elem == 7); // F,G,H,L,M,N,K
        P.zeros();
        // params = (F, G, H, L, M, N, K)
        P(0, 0) = params(1) + params(2) + params(6) / 9.; // P_11 = G+H+K/9
        P(1, 1) = params(0) + params(2) + params(6) / 9.; // P_22 = F+H+K/9
        P(2, 2) = params(0) + params(1) + params(6) / 9.; // P_33 = F+G+K/9
        P(0, 1) = -1. * params(2) + params(6) / 9.;       // P_12 = -H+K/9
        P(1, 0) = -1. * params(2) + params(6) / 9.;       // P_12 = -H+K/9
        P(0, 2) = -1. * params(1) + params(6) / 9.;       // P_13 = -G+K/9
        P(2, 0) = -1. * params(1) + params(6) / 9.;       // P_13 = -G+K/9
        P(1, 2) = -1. * params(0) + params(6) / 9.;       // P_23 = -F+K/9
        P(2, 1) = -1. * params(0) + params(6) / 9.;       // P_23 = -F+K/9
        P(3, 3) = 2. * params(5);                         // P_44 = N
        P(4, 4) = 2. * params(4);                         // P_55 = M
        P(5, 5) = 2. * params(3);                         // P_66 = L
    }

    // a quadratic-form criterion of the parameters (Hill, DFA, Ani)
    CriterionResponse quadratic_response(const vec &v, const vec &params, void (*fill)(const vec &, mat &), const int &level)
    {
        mat6 P;
        fill(params, P);
        return Eq_stress_P_response(v, P, level);
    }
    }

    mat P_Ani(const vec &params)
    {
        mat P(6, 6);
        fill_P_Ani(params, P);
        return P;
    }

    mat P_Hill(const vec &params)
    {
        mat P(6, 6);
        fill_P_Hill(params, P);
        return P;
    }

    mat P_DFA(const vec &params)
    {
        mat P(6, 6);
        fill_P_DFA(params, P);
        return P;
    }


    CriterionResponse Mises_stress_response(const vec &v, const int &level)
    {
        CriterionResponse out;
        double d[6];
        out.f = mises_voigt(v.memptr(), Tensor2Type::stress, d);
        if (level >= 1) out.df = flow_voigt(d, out.f, 3./2., 2.);
        if (level >= 2) mises_hessian_voigt(out.df, out.f, out.ddf, simcoon::iota);
        return out;
    }

    CriterionResponse Tresca_stress_response(const vec &v, const int &level)
    {
        assert(v.size() == 6);
        if (level >= 2) {
            throw std::invalid_argument("Tresca_stress_response: the Tresca surface has no Hessian (vertex flow)");
        }
        CriterionResponse out;
        const vec lambda = eig_sym(v2t_stress(v));   // ascending
        out.f = lambda(2) - lambda(0);
        if (level >= 1) out.df = eta_stress(v);
        return out;
    }

    CriterionResponse Eq_stress_P_response(const vec &v, const mat &H, const int &level)
    {
        assert(v.size() == 6);
        const vec6 Hv = H*v;
        const double q = dot(v, Hv);
        CriterionResponse out = zero_response(level);
        // round-off of q, bounded by |v|^T |H| |v| (Hv itself is round-off on a null direction of H)
        const double *h = H.memptr();
        double q_abs = 0.;
        for (int j = 0; j < 6; ++j) {
            double c = 0.;
            for (int i = 0; i < 6; ++i) c += std::fabs(h[6*j + i])*std::fabs(v(i));
            q_abs += c*std::fabs(v(j));
        }
        if (!(q > 8.*datum::eps*q_abs)) return out;
        out.f = sqrt(q);
        if (level < 1 || out.f < simcoon::iota) return out;
        out.df = Hv/out.f;
        if (level >= 2) out.ddf = (H - outer6(out.df.memptr(), out.df.memptr()))/out.f;
        return out;
    }

    CriterionResponse Drucker_stress_response(const vec &v, const double &b, const double &n, const int &level)
    {
        assert(b >= 0.);
        assert(n > 0.);
        return drucker_correction(v, Mises_stress_response(v, level), b, n, level);
    }

    CriterionResponse Drucker_ani_stress_response(const vec &v, const vec &params, const double &b, const double &n,
                                                  const int &level)
    {
        assert(b >= 0.);
        assert(n > 0.);
        mat6 P;
        fill_P_DFA(params, P);
        return drucker_correction(v, Eq_stress_P_response(v, P, level), b, n, level);
    }

    double Drucker_stress(const vec &v, const double &b, const double &n) { return Drucker_stress_response(v, b, n, 0).f; }
    vec dDrucker_stress(const vec &v, const double &b, const double &n) { return Drucker_stress_response(v, b, n, 1).df; }
    mat ddDrucker_stress(const vec &v, const double &b, const double &n) { return Drucker_stress_response(v, b, n, 2).ddf; }

    double Drucker_ani_stress(const vec &v, const vec &params, const double &b, const double &n) { return Drucker_ani_stress_response(v, params, b, n, 0).f; }
    vec dDrucker_ani_stress(const vec &v, const vec &params, const double &b, const double &n) { return Drucker_ani_stress_response(v, params, b, n, 1).df; }
    mat ddDrucker_ani_stress(const vec &v, const vec &params, const double &b, const double &n) { return Drucker_ani_stress_response(v, params, b, n, 2).ddf; }

    double Tresca_stress(const vec &v) { return Tresca_stress_response(v, 0).f; }
    vec dTresca_stress(const vec &v) { return Tresca_stress_response(v, 1).df; }

    double Eq_stress_P(const vec &v, const mat &H) { return Eq_stress_P_response(v, H, 0).f; }
    vec dEq_stress_P(const vec &v, const mat &H) { return Eq_stress_P_response(v, H, 1).df; }
    mat ddEq_stress_P(const vec &v, const mat &H) { return Eq_stress_P_response(v, H, 2).ddf; }

    double Hill_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Hill, 0).f; }
    vec dHill_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Hill, 1).df; }
    mat ddHill_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Hill, 2).ddf; }

    double DFA_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_DFA, 0).f; }
    vec dDFA_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_DFA, 1).df; }
    mat ddDFA_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_DFA, 2).ddf; }

    double Ani_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Ani, 0).f; }
    vec dAni_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Ani, 1).df; }
    mat ddAni_stress(const vec &v, const vec &params) { return quadratic_response(v, params, fill_P_Ani, 2).ddf; }

    double Eq_stress(const vec &v, const string &eq_type, const vec &param)
    {
        if (eq_type == "Mises")
        {
            return Mises_stress(v);
        }
        else if (eq_type == "Tresca")
        {
            return Tresca_stress(v);
        }
        else if (eq_type == "Drucker")
        {
            return Drucker_stress(v, param(0), param(1));
        }
        else if (eq_type == "Hill")
        {
            return Hill_stress(v, param);
        }
        else if (eq_type == "Ani")
        {
            return Ani_stress(v, param);
        }
        else if (eq_type == "DFA")
        {
            return DFA_stress(v, param);
        }
        else
        {
            // throw, not exit(0): an exit kills the host process (e.g. Python)
            throw std::invalid_argument("Eq_stress: unknown criterion '"
                                        + eq_type + "'");
        }
    }

    vec dEq_stress(const vec &v, const string &eq_type, const vec &param)
    {
        if (eq_type == "Mises")
        {
            return eta_stress(v);
        }
        else if (eq_type == "Tresca")
        {
            return dTresca_stress(v);
        }
        else if (eq_type == "Drucker")
        {
            return dDrucker_stress(v, param(0), param(1));
        }
        else if (eq_type == "Hill")
        {
            return dHill_stress(v, param);
        }
        else if (eq_type == "Ani")
        {
            return dAni_stress(v, param);
        }
        else if (eq_type == "DFA")
        {
            return dDFA_stress(v, param);
        }
        else
        {
            throw std::invalid_argument("dEq_stress: unknown criterion '"
                                        + eq_type + "'");
        }
    }

} // namespace simcoon
