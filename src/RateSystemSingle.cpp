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
// #include <stringstream>
// #define NDEBUG

size_t spatial_index_iterator = 0;


single_state_type::single_state_type(vector<size_t> P_sizes) {
    atomP.resize(P_sizes.size());
    for (size_t i = 0; i < atomP.size(); i++) {
        atomP[i].resize(P_sizes[i]);
    }
    cumulative_photo.resize(atomP.size());
}
single_state_type::single_state_type(){}


// Critical vector-space operators
single_state_type& single_state_type::operator+=(const single_state_type &s) {
    for (size_t r = 0; r < atomP.size(); r++) {
        cumulative_photo[r] += s.cumulative_photo[r];
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] += s.atomP[r][i];
        }
    }
    F += s.F;
    bound_charge += s.bound_charge;
    return *this;
}

single_state_type& single_state_type::operator*=(const double x) {
    for (size_t r = 0; r < atomP.size(); r++) {
        cumulative_photo[r] *= x;
        for (size_t i = 0; i < atomP[r].size(); i++) {
            atomP[r][i] *= x;
        }
    }
    F *= x;
    bound_charge *=x;
    return *this;
}

// convenience members
single_state_type& single_state_type::operator=(const double x) {
    for (auto&a : cumulative_photo)
        a=x;
    for (auto& P : atomP) {
        for (auto& p : P) {
            p=x;
        }
    }
    F = x;
    bound_charge = x;
    return *this;
}

// Resizes the container to fit all of the states present in the atom ensemble
void single_state_type::set_P_shape(const vector<RateData::Atom>& atomsys) {
    vector<size_t> P_shape; 
    for(auto atom: atomsys){
        if (atom.nAtoms_in_sims[spatial_index] > 0){
            P_shape.push_back(atom.num_conf);
        } 
    }
    set_P_shape(P_shape);
}

void single_state_type::set_P_shape(const vector<size_t>& shape) {
    atomP.resize(shape.size());
    for (size_t i = 0; i < shape.size(); i++) {
        atomP[i].resize(shape[i]);
    }
    cumulative_photo.resize(shape.size());
}

// Returns the L1 norm
double single_state_type::norm(size_t _c) const {
    double n = 0;
    for (auto& P : atomP) {
        for (auto& p : P) {
            n += fabs(p);
        }
    }
    n += F.norm(_c);
    return n;
}


// Intended usage: cout<<s.atomP[a]<<endl;
ostream& operator<<(ostream& os, const bound_t& bound) {
    const double units = 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;
    for (size_t i=0; i<bound.size(); i++) {
        os << bound[i]*units << " ";
    }
    return os;
}

ostream& operator<<(ostream& os, const single_state_type& st) {
    for (size_t a=0; a<st.atomP.size(); a++) {
        os << st.atomP[a];
        if (a != st.atomP.size()-1)
            os<<"| ";
    }
    os << st.F;
    return os;
}
