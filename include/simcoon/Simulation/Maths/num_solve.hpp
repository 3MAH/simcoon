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
#include <string>
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

/**
 * @brief Complementarity function of a local return-mapping Newton (the loops on Fischer_Burmeister_m).
 *
 * A numerical setting of the run, not a material parameter: solver_params::ncp and
 * state_variables::ncp carry its code, set by the @c ncp option of @c solve and @c sim.umat
 * (see the UMAT catalogue). Fischer-Burmeister, the historical smoothing, is the default.
 */
enum class NCPMethod : int {
    fischer_burmeister = 0,  ///< normalised Fischer-Burmeister function (Fischer_Burmeister_m)
    min_fb_fallback = 1      ///< min-NCP semi-smooth Newton on a single surface, Fischer-Burmeister once the error grows
};

/**
 * @brief The NCP method of a code (solver_params::ncp, state_variables::ncp), checked: 0 or 1.
 * @param code the code
 * @param where the caller, for the error message
 * @throw std::invalid_argument for any other value
 */
NCPMethod ncp_method_of(const int &code, const std::string &where);

/**
 * @brief One Newton step per call of a local complementarity solve, with the chosen function.
 *
 * NCPMethod::fischer_burmeister forwards to Fischer_Burmeister_m. NCPMethod::min_fb_fallback takes
 * the semi-smooth Newton step of \f$ \min(\Delta p^*_i, -\Phi_i) = 0 \f$ (\f$ \Delta p^*_i = \Delta p_i |B_{ii}| \f$):
 * on an active row (\f$ \Phi_i > 0 \f$ or \f$ -\Phi_i < \Delta p^*_i \f$) the linearised constraint
 * \f$ \Phi_i + B_{i\cdot}\,\delta p = 0 \f$, the classic cutting plane (Ortiz and Simo 1986), on an
 * inactive one the multiplier back to zero (primal-dual active set, Hintermueller, Ito and Kunisch
 * 2002). The first step from \f$ \Delta p = 0 \f$ is the full Newton step, where the Fischer-Burmeister
 * function only takes two thirds of it: a linear hardening returns in one step instead of six.
 * Its @c error is \f$ \sum_i 2|r_i|/|Y^{crit}_i| \f$ of the entering iterate, the scale of the
 * Fischer-Burmeister measure near the solution. It applies to a single surface (one row): a
 * system of several surfaces keeps the Fischer-Burmeister step, the active set of coupled
 * surfaces making the plain semi-smooth Newton cycle. When the error grows from one call to the
 * next, or the row is singular, the solve falls back to Fischer-Burmeister for the remaining
 * calls (a nonlinear kinematic hardening can make the full step overshoot).
 *
 * One object per local solve (its state is the fallback).
 */
class ncp_step {
public:
    explicit ncp_step(const NCPMethod &method = NCPMethod::fischer_burmeister) : method_(method) {}
    /// Same arguments and effects as Fischer_Burmeister_m: Dp += dp, error of the entering iterate
    void operator()(const arma::vec &Phi, const arma::vec &Y_crit, const arma::mat &denom, arma::vec &Dp, arma::vec &dp, double &error);
    /// Whether this solve fell back to Fischer-Burmeister
    bool fell_back() const { return fallback_; }
private:
    NCPMethod method_;
    bool fallback_ = false;
    double prev_error_ = -1.;
};
    
arma::mat denom_FB_m(const arma::vec &, const arma::mat &, const arma::vec &);
    

/** @} */ // end of maths group

} //namespace simcoon
