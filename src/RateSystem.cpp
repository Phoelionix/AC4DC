/**
 * @file RateSystem.cpp
 * @brief @copybrief RateSystem.h
 */
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

#include "RateSystem.h"
#include "Dipole.h"
#include <math.h>
#include <vector>
// #include <stringstream>
// #define NDEBUG

vector<vector<size_t>> state_type::sim_P_sizes; //= vector<vector<size_t>>(0);
size_t state_type::num_sims=0;
bool state_type::initialised=false;
bool state_type::active=false;

state_type::state_type() {
    assert(initialised == active);
    sims.resize(num_sims);
    set_P_shape(sim_P_sizes);
}

// Critical vector-space operators
state_type& state_type::operator+=(const state_type &s) {
    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]+=s.sims[i];
    }
    return *this;
}

state_type& state_type::operator*=(const double x) {
    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]*=x;
    }
    return *this;
}

state_type& state_type::operator*=(const std::vector<double> x) {
    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]*=x[i];
    }
    return *this;
}

// convenience members
state_type& state_type::operator=(const double x) {
    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]=x;
    }
    return *this;
}

state_type& state_type::operator=(const state_type &s) {
    assert(sims.size()==s.sims.size());
    for (size_t i = 0; i < sims.size(); i++) {
        sims[i]=s.sims[i];
    }
    return *this;
}

// Returns the L1 norm for each simulated volume's _c continuum
std::vector<double> state_type::norm(size_t _c) const {
    std::vector<double> return_vector;
    for(size_t i=0; i<sims.size();i++){
        double n = 0;
        for (auto& P : sims[i].atomP) {
            for (auto& p : P) {
                n += fabs(p);
            }
        }
        n += sims[i].F.norm(_c);
        return_vector.push_back(n);
    }
    return return_vector;
}