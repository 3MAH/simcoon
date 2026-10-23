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

///@file step_meca.hpp
///@brief object that defines a mechanical step
///@version 1.0

#pragma once

#include <iostream>
#include <armadillo>
#include "step.hpp"
#include "../Phase/state_variables_M.hpp"

namespace simcoon{

/**
 * @file step_meca.hpp
 * @brief Solver functions and classes.
 */

/** @addtogroup solver
 *  @{
 */


//======================================
class step_meca : public step
//======================================
{
private:
    
protected:
    
	public :

    arma::Col<int> cBC_meca; //1 is for stress (flux), 0 if for strain (state)
    arma::vec BC_meca;
    arma::mat mecas;
    arma::mat BC_mecas;
    arma::mat BC_w;
    arma::mat BC_R;
    double BC_T;
    int cBC_T;
    arma::vec Ts;
    arma::vec BC_Ts;
    
    step_meca(); 	//default constructor
    step_meca(const unsigned int &size); 	//constructor that allocates BC_meca and cBC_meca
    step_meca(const int &number, const double &Dn_init, const double &Dn_mini, const double &Dn_inc, const int &mode, const unsigned int &control_type, const arma::Col<int>& cBC_meca, const arma::vec& BC_meca, const arma::mat& mecas, const arma::mat& BC_mecas, const double& BC_T, const int& cBC_T, const arma::vec& Ts, const arma::vec& BC_Ts, const arma::mat& BC_w, const arma::mat& BC_R); //Constructor with parameters

    step_meca(const step_meca &s); 	//Copy constructor
    virtual ~step_meca();
    
    using step::generate;
    virtual void generate(const double&, const arma::vec &, const arma::vec&, const double&);
    /**
     * @brief Generate the increment tables for a fully kinematic step (control types 5/6).
     * @param[in] mTime current (absolute) time at the start of the step
     * @param[in] mF current 3x3 deformation gradient \f$ \mathbf{F} \f$ (the step anchor)
     * @param[in] mT current temperature
     * Modes 1/2 interpolate geodesically between mF and BC_meca via the real matrix
     * logarithm/exponential; mode 3 reads the 9-component F table (absolute time column).
     */
    virtual void generate_kin(const double&, const arma::mat &, const double &);
    virtual void assess_inc(const double &, double &, const double &, phase_characteristics &, double &, const double &, const arma::mat &, const int &);
    
    virtual step_meca& operator = (const step_meca &s);

    friend  std::ostream& operator << (std::ostream& os, const step_meca &s);
};


/**
 * @brief Increments of a kinematic step (control types 5/6, modes 1/2) on the geodesic from F_prev
 * to the target: F_i = exp(w_i ln(F_target F_prev^-1)) F_{i-1}, weights w_i = inc_coef / sum(inc_coef).
 * Control type 6 passes F_target = I + grad u.
 * Shared by step_meca and step_thermomeca.
 * @param[in] F_prev deformation gradient at the start of the step
 * @param[in] BC_meca target F, 9 components row-major
 * @param[in] inc_coef increment weights (ones, or the sinusoidal mode's)
 * @param[in] BC_Time duration of the step
 * @param[out] times time increments (sized ninc)
 * @param[out] mecas generator of each increment, \f$ \Delta\mathbf{L}_i = w_i \ln(\mathbf{F}_{target}\mathbf{F}_{prev}^{-1}) \f$
 *             (9 components row-major): a fraction s of increment i is \f$ \exp(s\,\Delta\mathbf{L}_i)\,\mathbf{F} \f$
 *             (step cuts, set_driven_strain in the solver)
 * @param[out] BC_mecas F at the end of each increment
 */
void geodesic_F_increments(const arma::mat &F_prev, const arma::vec &BC_meca, const arma::vec &inc_coef, const double &BC_Time,
                           arma::vec &times, arma::mat &mecas, arma::mat &BC_mecas);

/**
 * @brief The tabulated (mode 3) kinematic step in the form geodesic_F_increments gives: BC_mecas
 * holds F at the end of each increment and mecas the generator \f$ \Delta\mathbf{L}_i =
 * \ln(\mathbf{F}_i\mathbf{F}_{i-1}^{-1}) \f$ of the increment, the row's F reached on the
 * geodesic between two table rows. Shared by step_meca and step_thermomeca.
 * @param[in] F_prev deformation gradient at the start of the step
 * @param[in] cBC_meca 2 for a component absent from the table: held at its previous value
 * @param[in] control_type 5 (the table gives F) or 6 (it gives grad u = F - I)
 * @param[out] mecas generators, one row per increment (9 components row-major)
 * @param[in,out] BC_mecas in: the table values; out: F, one row per increment
 */
void tabulated_F_increments(const arma::mat &F_prev, const arma::Col<int> &cBC_meca, const unsigned int &control_type,
                            arma::mat &mecas, arma::mat &BC_mecas);

/** @} */ // end of solver group

} //namespace simcoon
