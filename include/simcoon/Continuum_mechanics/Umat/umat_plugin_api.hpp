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

#pragma once
#include <string>

class umat_plugin_api {
public:
    virtual std::string name() const = 0;
    virtual void umat_(double *, double *, double *, double &, double &, double &, double &, double *, double *, double &, const double *, const double *, const double *, const double &, const double &, const double &, const double &, const double &, char *, const int &, const int &, const int &, const int &, const double *, const int &, const double &, const double *, double &, const double &, const double *, const double *, const int &, const int &, const double &, const int &, const int &, const int &) = 0;

   virtual ~umat_plugin_api() {}
};
