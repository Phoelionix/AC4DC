#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "AdamsIntegrator.hpp"

using namespace std;



/*
class state_type :
    boost::additive1< state_type ,
    boost::additive2< state_type , double ,
    boost::multiplicative2< state_type , double > > >
    */


// Class responsible for storing the system state.
class state_type
{
public:
    typedef std::vector<double> bound_t; // Probabilities of state

    std::vector<bound_t> atomP; // Probabilities of state for all atoms.
    std::vector<double> f; // Energy distribution function

    state_type(){
        atomP.resize(P_sizes.size());
        for (size_t i = 0; i < atomP.size(); i++) {
            atomP[i].resize(P_sizes[i]);
        }
        f.resize(f_grid.size());
    }


    // Critical vector-space devices
    state_type& operator+=(const state_type &s){
        for (size_t r = 0; r < atomP.size(); r++) {
            for (size_t i = 0; i < atomP[r].size(); i++) {
                atomP[r][i] += s.atomP[r][i];
            }
        }
        for (size_t i = 0; i < f.size(); i++) {
            f[i] += s.f[i];
        }
        return *this;
    }

    state_type& operator*=(const double x){
        for (size_t r = 0; r < atomP.size(); r++) {
            for (size_t i = 0; i < atomP[r].size(); i++) {
                atomP[r][i] *= x;
            }
        }
        for (size_t i = 0; i < f.size(); i++) {
            f[i] *= x;
        }
        return *this;
    }

    state_type operator+(const state_type& s2){
        state_type retval = *this;
        retval += s2;
        return retval;
    }

    state_type operator*(double x){
        state_type retval = *this;
        retval *= x;
        return retval;
    }

    // convenience members
    state_type& operator=(const double x){
        for (auto& P : atomP){
            for (auto& p : P) {
                p=x;
            }
        }
        for (auto& d : f){
            d=x;
        }
        return *this;
    }

    static void print_info();

    // Resizes the container to fit all of the states present in the atom ensemble
    static void set_P_shape(const vector<RateData::Atom>& atomsys) {
        P_sizes.resize(atomsys.size());
        // make the P's the right size lmao
        for (size_t a = 0; a < atomsys.size(); a++) {
            P_sizes[a] = atomsys[a].num_conf;
        }
    }

    // defines the binning of f
    static void set_elec_points(size_t n, double min_e, double max_e);
    // Defines number and style of atomP
    static void set_P_shape(const vector<size_t>& shape){
        P_sizes = shape;
    }

    static vector<double> f_grid; // f energies
    static vector<double> f_widths; // bin widths
private:
    static vector<size_t> P_sizes;
};

ostream& operator<<(ostream& os, const state_type& st);



/*         TODO: Fix these (do not currently habdle multiple P arrays implemented above.)
// Algebra definition for error-controlled steppers
state_type operator/( const state_type &s1 , const state_type &s2 ){
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

state_type abs( const state_type &s){
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
