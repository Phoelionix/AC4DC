#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
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

    // critical for definiton of vector-space algebra
    state_type(){
        // Iterates over all atomic species' P's
        for(bound_t P : atomP) {
            P.reserve(0);
        }
        f.reserve(0);
    }

    state_type(const vector<int> P_sizes, const size_t f_size){
        assert(P_sizes.size() == atomP.size());
        for (size_t i = 0; i < atomP.size(); i++) {
            atomP[i].reserve(P_sizes[i]);
        }
        f.reserve(f_size);
    }

    // Resizes the container to fit all of the states present in the atom ensemble
    state_type(const vector<RateData::Atom>& atomsys, size_t f_size) {
        atomP.reserve(atomsys.size());
        // make the P's the right size lmao
        for (size_t a = 0; a < atomsys.size(); a++) {
            atomP[a].resize(atomsys[a].num_conf);
        }
        f.reserve(f_size);
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

    // convenience members
    state_type& operator=(const double x){
        for (size_t r = 0; r < atomP.size(); r++) {
            for (size_t i = 0; i < atomP[r].size(); i++) {
                atomP[r][i] = x;
            }
        }
        for (size_t i = 0; i < f.size(); i++) {
            f[i] = x;
        }
        return *this;
    }



protected:
    // defines the binning of f
    static vector<double> f_grid; // f energies
    static vector<double> f_widths; // bin widths
};

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
