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


void state_type::transform_basis_all(std::vector<double> new_knots){
    int new_basis_order = BasisSet::BSPLINE_ORDER;
    //// Get knots that have densities
    int num_new_splines = static_cast<int>(Distribution::get_trimmed_knots(new_knots).size());   // TODO replace get_trimmed_knots with get_num_funcs?. 
    std::vector<std::vector<std::vector<std::vector<double>>>> new_densities_container_all;
    for (size_t V =0; V < sims.size(); V++){
        std::vector<std::vector<std::vector<double>>> new_densities_container;
        for (size_t _c = 0; _c < Distribution::num_continuums; _c++){
            //// Compute densities for knots
            std::vector<std::vector<double>> new_densities(num_new_splines, std::vector<double>(64, 0));

            // Cackle and iterate through each new spline.
            for (size_t i=0; static_cast<int>(i)<num_new_splines; i++){
                // Use current basis to generate the density terms for gaussian integration at for each basis point.   
                // Black magic. ଘ(੭ˊᵕˋ)੭.*･｡ﾟ
                double a = new_knots[i];                  // i.e. <new_basis>.supp_min(i);
                double b = new_knots[i+new_basis_order];  // i.e. <new_basis>.supp_max(i);        
                for(size_t j=0; j < 64; j++){
                    double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;
                    new_densities[i][j] = (sims[V].F)(_c,e);  
                }
            }
            // Change distribution to empty one in new basis.    
            new_densities_container.push_back(new_densities);
        }
        for (size_t _c = 0; _c < Distribution::num_continuums; _c++){
            vector<double> new_f(num_new_splines,0);
            sims[V].F[_c] = new_f; 
        }
        new_densities_container_all.push_back(new_densities_container);
    }
    Distribution::load_knot(new_knots);
    for (size_t V =0; V < sims.size(); V++){
        for (size_t _c = 0; _c < Distribution::num_continuums; _c++){
            sims[V].F[_c].resize(Distribution::size);
            // Add densities in new basis.
            sims[V].F.add_density_distribution(_c, new_densities_container_all[V][_c]);

        }
    }
}