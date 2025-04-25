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

///@file exception.hpp
///@brief exceptions of simcoon
///@version 1.0

#pragma once

#include <stdexcept>
#include <string>

namespace simcoon {

class exception_eig_sym : public std::runtime_error {
public:
    explicit exception_eig_sym(const std::string& msg)
        : std::runtime_error(msg) {}
};

class exception_det : public std::runtime_error {
    public:
        explicit exception_det(const std::string& msg)
            : std::runtime_error(msg) {}
    };

class exception_inv : public std::runtime_error {
    public:
        explicit exception_inv(const std::string& msg)
            : std::runtime_error(msg) {}
    };    

class exception_sqrtmat_sympd : public std::runtime_error {
    public:
        explicit exception_sqrtmat_sympd(const std::string& msg)
            : std::runtime_error(msg) {}
    };
    
class exception_logmat_sympd : public std::runtime_error {
    public:
        explicit exception_logmat_sympd(const std::string& msg)
            : std::runtime_error(msg) {}
    };    
    
class exception_expmat_sym : public std::runtime_error {
    public:
        explicit exception_expmat_sym(const std::string& msg)
            : std::runtime_error(msg) {}
    };

class exception_powmat : public std::runtime_error {
    public:
        explicit exception_powmat(const std::string& msg)
            : std::runtime_error(msg) {}
    };    
    
} //namespace simcoon
