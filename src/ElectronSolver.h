#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

/*
This file should arguably be called RateEquationSOlver, however, for historical reasons it is not.
*/
// #include <boost/numeric/odeint.hpp>
#include "AdamsIntegrator.hpp"
#include "RateSystem.h"
#include "Constant.h"
#include "MolInp.h"
#include "Input.h"
#include <iostream>
#include <fstream>



// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

class PhotonFlux
{
public:
    PhotonFlux(){};
    PhotonFlux(double fluence, double fwhm){
        set_parameters(fluence, fwhm);
    };
    void set_parameters(double, double);
    inline double operator()(double t); // Yields Photon flux in J/cm^2
private:
    double A;
    double B;
};

// W_ij for the rate equaitons
class Weight
{
public:
    Weight(size_t size);
    ~Weight();
    double* W;
    double from(size_t idx); // returns summed transition rate away from supplied index
    double to(size_t idx); // returns summed transition rate towards supplied index
    double& operator()(size_t to, size_t from){
        assert(from >= 0 && from<size);
        assert(to >= 0 && to<size);
        return W[to*size + from];
    }
protected:
    size_t size;
};


class ElectronSolver : private MolInp, private Adams_BM<state_type>
{
public:
    ElectronSolver(const char* filename, std::ofstream& log);
    void solve();
    void print(const std::string& fname);
private:
    double timespan;
    // Model parameters
    PhotonFlux pf;
    void set_flux(double fluence_in_Jcm2);
    void set_initial_conditions();
    void sys(const state_type& s, state_type& sdot, const double t);
};



#endif /* end of include guard: SYS_SOLVER_CXX_H */
