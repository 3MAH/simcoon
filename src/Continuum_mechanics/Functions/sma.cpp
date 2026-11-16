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

///@file sma.cpp
///@brief Constitutive functions shared by the unified SMA laws

#include <cmath>
#include <stdexcept>
#include <string>
#include <armadillo>
#include <simcoon/parameter.hpp>
#include <simcoon/Continuum_mechanics/Functions/contimech.hpp>
#include <simcoon/Continuum_mechanics/Functions/criteria.hpp>
#include <simcoon/Simulation/Maths/num_solve.hpp>
#include <simcoon/Continuum_mechanics/Functions/sma.hpp>

namespace simcoon {

arma::vec sma_transformation_flow(const arma::vec &s, const double &sigmacrit, const double &Hmin,
                                  const double &Hmax, const double &k1, const bool &aniso_criteria,
                                  const arma::vec &DFA_params, const double &prager_b, const double &prager_n) {
    double sstar = Mises_stress(s) - sigmacrit;
    if (sstar < 0.) {
        sstar = 0.;
    }
    const double Hc = Hmin + (Hmax - Hmin) * (1. - std::exp(-1. * k1 * sstar));
    const CriterionResponse dr = aniso_criteria ? Drucker_ani_stress_response(s, DFA_params, prager_b, prager_n, 1)
                                                : Drucker_stress_response(s, prager_b, prager_n, 1);
    return Hc * dr.df;
}

NCPMethod sma_ncp_method(const int &ncp, const int &tangent_mode, const std::string &where) {
    const NCPMethod method = ncp_method_of(ncp, where);
    if (ncp_is_semismooth(method) && !integrates_closest_point(tangent_mode)) {
        throw std::invalid_argument(where + ": ncp='min' / 'active_set' need integrator='closest_point' "
                                    "(the cutting plane keeps the transformation strain of its iterates)");
    }
    return method;
}

SMAHardening sma_hardening(const double &xi, const double &a1, const double &a2, const double &a3,
                           const double &n1, const double &n2, const double &n3, const double &n4) {
    const double x = std::min(std::max(xi, 0.), 1.);
    const bool inside = (xi > 0.) && (xi < 1.);
    SMAHardening h;
    if ((n1 == 1.) && (n2 == 1.)) {
        h.F = a1*x + a3;
        h.dF = inside ? a1 : 0.;
    }
    else {
        h.F = 0.5*a1*(1. + std::pow(x, n1) - std::pow(1. - x, n2)) + a3;
        h.dF = inside ? 0.5*a1*(n1*std::pow(xi, n1 - 1.) + n2*std::pow(1. - xi, n2 - 1.)) : 0.;
    }
    if ((n3 == 1.) && (n4 == 1.)) {
        h.R = a2*x - a3;
        h.dR = inside ? a2 : 0.;
    }
    else {
        h.R = 0.5*a2*(1. + std::pow(x, n3) - std::pow((1. - x), n4)) - a3;
        h.dR = inside ? 0.5*a2*(n3*std::pow(xi, n3 - 1.) + n4*std::pow(1. - xi, n4 - 1.)) : 0.;
    }
    return h;
}

} // namespace simcoon
