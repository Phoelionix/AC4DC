/**
 * @file RateSystem.h
 * @brief 
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

#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "FreeDistribution.h"
#include "RateSystemSingle.h"




/**
 * @brief Class responsible for storing the system state at a moment in time.
 * @details y[i] corresponds to t[i] throughout the plasma code, which is somewhat unfortunate.
 * y[i].F and y.atomP are everything you'd want to pop on a plot's axes for time t[i]... perhaps, even, on the y-axis against t? -S.P.
 */
class state_type
{
public:
    /// Individual rate systems
    std::vector<single_state_type> sims; 

    // Num sims is just the number of voxels or whatever volumes are being simulated, connected through SpatialConnection objects (transfer of free electrons). 
    state_type();

    // Critical vector-space devices
    state_type& operator+=(const state_type &s);
    state_type& operator*=(const double x);
    state_type& operator*=(const std::vector<double> x);
    // state_type operator+(const state_type& s2);
    // state_type operator*(double x);
    // convenience members
    state_type& operator=(const double x);
    state_type& operator=(const state_type &s);
    // state_type& operator=(const state_type& s2);

    inline single_state_type& operator[](size_t n) {
        return this->sims[n];
    }

    std::vector<double> norm(size_t _c) const;

    // Defines number and style of atomP
    // Resizes the container to fit all of the states present in the atom ensemble
    // void set_P_shape(const vector<RateData::Atom>& atomsys){
    //     //sim_P_sizes.resize(0);
    //     for (auto& sim: sims){
    //         sim.set_P_shape(atomsys);
    //         sim_P_sizes.push_back(sim.get_P_sizes());
    //     }
    // }
   void set_P_shape(const vector<vector<size_t>>& sim_shapes) {
        sims.resize(sim_shapes.size());
        for(size_t V=0; V<sim_shapes.size(); V++){
            sims[V].set_P_shape(sim_shapes[V]);
        }
    }
    void update_P_shape(){
        set_P_shape(sim_P_sizes); 
    }
    static void initialise_P_shape(const vector<RateData::Atom>& atomsys){
        assert(num_sims>0); // check num sims has been set.
        sim_P_sizes.resize(0);
        for(size_t V=0; V<num_sims; V++){
            vector<size_t> P_shape; 
            for(auto atom: atomsys){
                //if (atom.nAtoms_in_sims[V] > 0){  // Can't do this because of way precomputed Q (e.g. Q_eii) is stored.
                P_shape.push_back(atom.num_conf);
                //} 
            }
            sim_P_sizes.push_back(P_shape);
        }
    }

    static void initialise_num_sims(const size_t& num_sims){
        state_type::num_sims = num_sims; 

    }
    static size_t get_num_sims() {
        return num_sims;
    }
    static void initialise(const vector<RateData::Atom>& atom_sys, size_t num_sims){
        assert(active);
        initialise_num_sims(num_sims); // must be before initialising P shape.
        initialise_P_shape(atom_sys);
        initialised=true;
    }


    size_t P_size(size_t V, size_t a) {
        return sims[V].get_P_sizes()[a];//sim_P_sizes[V][a];
    }
    size_t num_atoms(size_t V) {
        return sims[V].get_P_sizes().size();
        //return sim_P_sizes[V].size();
    }

    
    Distribution& get_sampleF(){return sims[0].F;}
    
    void set_sample_index(const size_t& sample_index){
        state_type::sample_index = sample_index;
    }

    static size_t Num_Sims(){
        return num_sims;
    }

    static void Mark_Active(){
        active=true;
    }

private:
    static bool initialised;
    static vector<vector<size_t>> sim_P_sizes;
    static size_t num_sims;
    static size_t sample_index;
    static bool active;
};



// All f integrals have the form
// df(e)/dt = Q [f] (e)
// Expand in some basis f = \sum_k a_k(t) f_k
// \sum_k da_k(t)/dt f_k(e) = Q[f] (e)
// This does not make life easier, unless we take the f_k to have compact support.


/*         TODO: Fix these (do not currently habdle multiple P arrays implemented above.)
// Algebra definition for error-controlled steppers
state_type operator/( const state_type &s1 , const state_type &s2 ) {
    std::vector<double> tmpP = s1.P;
    std::vector<double> tmpf = s1.f;

    for (size_t i = 0; i < tmpP.size(); i++) {
        tmpP[i] /= s2.P[i];
    }
    for (size_t i = 0; i < tmpf.size(); i++) {
        tmpf[i] /= s2.f[i];
    }
    return state_type( tmpP, tmpf );
}

state_type abs( const state_type &s) {
    state_type tmp(s.P.size(), s.f.size());

    for (size_t i = 0; i < tmp.P.size(); i++) {
        tmp.P[i] = std::abs(tmp.P[i]);
    }
    for (size_t i = 0; i < tmp.f.size(); i++) {
        tmp.f[i] = std::abs(tmp.f[i]);
    }
    return state_type( tmp );
}

// also only for steppers with error control
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< state_type >
{
    typedef double result_type;
    double operator()( const state_type &p ) const
    {
        using std::abs;
        double max=0.;
        double tmp;
        for (size_t i = 0; i < p.P.size(); i++) {
            tmp = abs(p.P[i]);
            max = (max > tmp) ? max : tmp;
        }
        for (size_t i = 0; i < p.f.size(); i++) {
            tmp = abs(p.f[i]);
            max = (max > tmp) ? max : tmp;
        }
        return max;
    }
};
*/

//
// // Flag no resizing
// namespace boost { namespace numeric { namespace odeint {
// template<>
// struct is_resizeable< state_type >
// {
//     typedef boost::false_type type;
//     const static bool value = type::value;
// };
//
// } } }

// End of definitions

#endif /* end of include guard: RATESYSTEM_CXX_H */
