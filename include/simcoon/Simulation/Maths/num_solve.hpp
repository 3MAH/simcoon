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
    
void Newton_Raphon(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_limits(const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_m(const arma::vec &, const arma::vec &, const arma::mat &, arma::vec &, arma::vec &, double &);

void Fischer_Burmeister_m_limits(const arma::vec &, const arma::vec &, const arma::vec &, const arma::mat &, const arma::mat &, arma::vec &, arma::vec &, double &);
    
arma::mat denom_FB_m(const arma::vec &, const arma::mat &, const arma::vec &);
    
} //namespace simcoon
