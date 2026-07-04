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

///@brief Canonical normalised Fischer-Burmeister residual: sum_i |FB_i| / |Y_crit_i| with
///Dp scaled by |diag(denom)|. This is THE definition of the FB error — Fischer_Burmeister_m
///reports it and the closest-point return mapping (return_mapping.cpp) uses it as its
///convergence/backtracking merit; keep both on this single implementation.
double Fischer_Burmeister_residual(const arma::vec &Phi, const arma::vec &Dp, const arma::mat &denom, const arma::vec &Y_crit);

void Fischer_Burmeister_m(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_m_limits(const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::vec &, double &);
    
arma::mat denom_FB_m(const arma::vec &, const arma::mat &, const arma::vec &);
    

/** @} */ // end of maths group

} //namespace simcoon
