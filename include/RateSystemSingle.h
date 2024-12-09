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

#ifndef RATESYSTEMSINGLE_CXX_H
#define RATESYSTEMSINGLE_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "FreeDistribution.h"



/**
 * @brief Class responsible for storing the system state at a moment in time.
 * @details y[i] corresponds to t[i] throughout the plasma code, which is somewhat unfortunate.
 * y[i].F and y.atomP are everything you'd want to pop on a plot's axes for time t[i]... perhaps, even, on the y-axis against t? -S.P.
 */
class single_state_type
{
public:
    /// Probabilities of state for all atoms.
    std::vector<bound_t> atomP; 
    // Tracks sum total of photoionisation for all atoms.
    std::vector<double> cumulative_photo;     
    /// Energy distribution function
    Distribution F;   
    double bound_charge;
    

    // Since we have removed P_sizes from this class, can use this method to fetch it dynamically. 
    std::vector<size_t> get_P_sizes(){
        std::vector<size_t> P_sizes;
        P_sizes.resize(atomP.size());
        for (size_t a = 0; a < atomP.size(); a++) {
            P_sizes[a] = atomP[a].size();
        }
        return P_sizes;
    }


    // P_sizes gives the number of configs for each atom.
    single_state_type(vector<size_t> P_sizes);
    single_state_type();

    // Critical vector-space devices
    single_state_type& operator+=(const single_state_type &s);
    single_state_type& operator*=(const double x);
    // single_state_type operator+(const single_state_type& s2);
    // single_state_type operator*(double x);
    // convenience members
    single_state_type& operator=(const double x);
    // single_state_type& operator=(const single_state_type& s2);

    double norm(size_t _c) const;


    // Defines number and style of atomP
    // Resizes the container to fit all of the states present in the atom ensemble
    void set_P_shape(const vector<RateData::Atom>& atomsys);
    void set_P_shape(const vector<size_t>& shape);

    void set_spatial_index(){
        spatial_index = spatial_index_incrementor++;
    }

protected:
    size_t spatial_index; 
private:
    static size_t spatial_index_incrementor;
};

ostream& operator<<(ostream& os, const single_state_type& st);
ostream& operator<<(ostream& os, const bound_t& dist);
ostream& operator<<(ostream& os, const Distribution& dist);

// All f integrals have the form
// df(e)/dt = Q [f] (e)
// Expand in some basis f = \sum_k a_k(t) f_k
// \sum_k da_k(t)/dt f_k(e) = Q[f] (e)
// This does not make life easier, unless we take the f_k to have compact support.


/*         TODO: Fix these (do not currently habdle multiple P arrays implemented above.)
// Algebra definition for error-controlled steppers
single_state_type operator/( const single_state_type &s1 , const single_state_type &s2 ) {
    std::vector<double> tmpP = s1.P;
    std::vector<double> tmpf = s1.f;

    for (size_t i = 0; i < tmpP.size(); i++) {
        tmpP[i] /= s2.P[i];
    }
    for (size_t i = 0; i < tmpf.size(); i++) {
        tmpf[i] /= s2.f[i];
    }
    return single_state_type( tmpP, tmpf );
}

single_state_type abs( const single_state_type &s) {
    single_state_type tmp(s.P.size(), s.f.size());

    for (size_t i = 0; i < tmp.P.size(); i++) {
        tmp.P[i] = std::abs(tmp.P[i]);
    }
    for (size_t i = 0; i < tmp.f.size(); i++) {
        tmp.f[i] = std::abs(tmp.f[i]);
    }
    return single_state_type( tmp );
}

// also only for steppers with error control
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< single_state_type >
{
    typedef double result_type;
    double operator()( const single_state_type &p ) const
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
// struct is_resizeable< single_state_type >
// {
//     typedef boost::false_type type;
//     const static bool value = type::value;
// };
//
// } } }

// End of definitions

#endif /* end of include guard: RATESYSTEM_CXX_H */
