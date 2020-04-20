#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

#include <boost/numeric/odeint>

class state_type : boost::additive1< state_type , boost::additive2< state_type , double , boost::multiplicative2< state_type , double > > >{
    vector<double> P; // Probabilities of state
    vector<double> f; // Energy distribution function

    state_type(const std::size_t P_size, const std::size_t f_size){
        P.reserve(P_size);
        f.reserve(f_size);
    }

    state_type(const vector<double> _P&, const vector<double> _f&) : P(_P), f(_f){

    }

    state_type& operator+=(const state_type &s){
        for (size_t i = 0; i < P.size(); i++) {
            P[i] += s.P[i];
        }
        for (size_t i = 0; i < f.size(); i++) {
            f[i] += s.f[i];
        }
        return *this;
    }

    state_type& operator*=(const double x){
        for (size_t i = 0; i < P.size(); i++) {
            P[i] *= x;
        }
        for (size_t i = 0; i < f.size(); i++) {
            f[i] *= x;
        }
        return *this;
    }
};

// Algebra definition for error-controlled steppers
state_type operator/( const state_type &s1 , const state_type &s2 ){
    vector<double> tmpP(s1.P);
    vector<double> tmpf(s1.f);

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

    for (size_t i = 0; i < tmpP.size(); i++) {
        tmp.P[i] = std::abs(tmpP[i]);
    }
    for (size_t i = 0; i < tmpf.size(); i++) {
        tmp.f[i] = std::abs(tmpf[i]);
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
        using std::max;
        using std::abs;
        double max=0.;
        double tmp;
        for (size_t i = 0; i < p.P.size(); i++) {
            tmp = abs(p.P[i])
            max = (max > tmp) ? max : tmp;
        }
        for (size_t i = 0; i < p.f.size(); i++) {
            tmp = abs(p.f[i])
            max = (max > tmp) ? max : tmp;
        }
        return max;
    }
};
} } }

// End of definitions

using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

class SystemSolver
{
public:
    SystemSolver(double final_time, double _dt);
    Solve();
    vector<double> T;
    vector<state_type> Y;
private:
    //                       steps, State,   value,  Derivative, Time,   Algebra
    adams_bashforth_moulton< 5 , state_type, double, state_type, double, vector_space_algebra > abm;
    // The important updater
    int num_steps;
    double dt;
    void sys(const state_type& s, state_type& sdot, const double t);

};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
