#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "AdamsIntegrator.hpp"



/*
class state_type :
    boost::additive1< state_type ,
    boost::additive2< state_type , double ,
    boost::multiplicative2< state_type , double > > >
    */

// Represents a statistical distribution of electrons.
// grid defines spacing in Hartree.
class Distribution
{
public:
    Distribution() {
        f.resize(size);
    }

    double operator[](size_t n){
        return this->f[n];
    }

    double operator()(double e){
        return this->f[i_from_e(e)];
    }

    // vector-space algebra
    Distribution& operator+=(const Distribution& d){
        for (size_t i=0; i<size; i++) {
            f[i] += d.f[i];
        }
        return *this;
    }

    Distribution& operator*=(double x){
        for (auto& fi : f){
            fi *= x;
        }
        return *this;
    }

    Distribution& operator=(const Distribution& d){
        f = d.f;
        return *this;
    }

    Distribution& operator=(double y){
        for (auto& fi : f){
            fi=y;
        }
        return *this;
    }

    template<typename T>
    T expect(T (& g) (double e)){
        T tmp = 0;
        for (size_t i = 0; i < size; i++) {
            tmp += g(grid[i])*f[i]*widths[i];
        }
        return tmp;
    }

    double eii_int (const CustomDataType::EIIdata& eii, const int i) const;
    double tbr_int (const CustomDataType::EIIdata& eii, const int i) const;
    void add_Qee(const Distribution& F);

    // N is the (dimensionless) count of the integrated energy desired
    void addDeltaLike(double e, double N);

    // Sets the object to have a MB distribution
    void set_maxwellian(double N, double T);
    // defines the f interpolation
    static void set_elec_points(size_t n, double min_e, double max_e);
    static std::string get_energies(std::string units="eV");
    vector<double> f;
    static size_t size;
private:
    static vector<double> grid;

    static vector<double> widths;
    static double e_from_i(size_t i);
    static size_t i_from_e(double e);
    static double min_e, max_e;
};

// Class responsible for storing the system state.
class state_type
{
public:
    typedef std::vector<double> bound_t; // Probabilities of state

    std::vector<bound_t> atomP; // Probabilities of state for all atoms.
    Distribution F; // Energy distribution function

    state_type();

    // Critical vector-space devices
    state_type& operator+=(const state_type &s);
    state_type& operator*=(const double x);
    // state_type operator+(const state_type& s2);
    // state_type operator*(double x);
    // convenience members
    state_type& operator=(const double x);
    // state_type& operator=(const state_type& s2);

    // Defines number and style of atomP
    // Resizes the container to fit all of the states present in the atom ensemble
    static void set_P_shape(const vector<RateData::Atom>& atomsys);
    static void set_P_shape(const vector<size_t>& shape){
        P_sizes = shape;
    }

private:
    static vector<size_t> P_sizes;
};

ostream& operator<<(ostream& os, const state_type& st);
ostream& operator<<(ostream& os, const state_type::bound_t& dist);
ostream& operator<<(ostream& os, const Distribution& dist);

// All f integrals have the form
// df(e)/dt = Q [f] (e)
// Expand in some basis f = \sum_k a_k(t) f_k
// \sum_k da_k(t)/dt f_k(e) = Q[f] (e)
// This does not make life easier, unless we take the f_k to have compact support.


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
