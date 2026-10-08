/* This file is part of simcoon private.
 
 Only part of simcoon is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This file is not be distributed under the terms of the GNU GPL 3.
 It is a proprietary file, copyrighted by the authors
 */

///@file num_solve.hpp
///@brief random number generators
///@author Chemisky

#pragma once
#include <armadillo>

namespace simcoon{

/**
 * @file num_solve.hpp
 * @brief Mathematical utility functions.
 */

/** @addtogroup maths
 *  @{
 */

    
void Newton_Raphon(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_limits(const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_m(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

/**
 * @brief The convergence measure of Fischer_Burmeister_m without its Newton step.
 *
 * \f$ \sum_i |\phi_i| / |Y^{crit}_i| \f$ with the normalised Fischer-Burmeister function
 * \f$ \phi_i = \sqrt{\Phi_i^2 + \Delta p_i^{*2}} + \Phi_i - \Delta p_i^* \f$,
 * \f$ \Delta p_i^* = \Delta p_i\,|B_{ii}| \f$ — the very expression (branches included)
 * Fischer_Burmeister_m returns in its @c error argument for the entering iterate, so a
 * line search (closest_point_return_mapping) can rate a trial point without solving the
 * \f$ N \times N \f$ system. Kept in step with Fischer_Burmeister_m by construction: that
 * function's error loop is this one.
 *
 * @param[in] Phi constraints (N)
 * @param[in] Y_crit normalisations (N, non-zero)
 * @param[in] denom local Jacobian (N x N); only |diag| is used
 * @param[in] Dp multipliers (N)
 * @return the residual
 */
double Fischer_Burmeister_residual(const arma::vec &Phi, const arma::vec &Y_crit, const arma::mat &denom, const arma::vec &Dp);

void Fischer_Burmeister_m_limits(const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::vec &, double &);
    
arma::mat denom_FB_m(const arma::vec &, const arma::mat &, const arma::vec &);
    

/** @} */ // end of maths group

} //namespace simcoon
