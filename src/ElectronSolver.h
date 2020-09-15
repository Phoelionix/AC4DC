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
#include <stdexcept>

// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

class PhotonFlux
{
public:
    PhotonFlux() {};
    PhotonFlux(double fluence, double fwhm) {
        set_pulse(fluence, fwhm);
    };
    void set_pulse(double, double);
    inline double operator()(double t); // Yields Photon flux in same units as A
    void save(const vector<double>& T, const std::string& file);
private:
    double A;
    double B;
};

#define MAX_T_PTS 1e6

class ElectronSolver : private ode::Adams_BM<state_type>
{
public:
    ElectronSolver(const char* filename, std::ofstream& log);
    void solve();
    void save(const std::string& folder);
    // Expose the underlying MolInp command
    void compute_cross_sections(std::ofstream& _log, bool recalc=true);
private:
    MolInp input_params;
    double timespan_au; // Atomic units
    // Model parameters
    PhotonFlux pf;

    // arrays computed at class initialisation
    vector<vector<eiiGraph> > RATE_EII;
    vector<vector<eiiGraph> > RATE_TBR;

    void get_energy_bounds(double& max, double& min);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_flux(double Jcm2_per_Haa02);
    void set_initial_conditions();
    // Components of sys that can be preallocated
    Eigen::VectorXd vec_dqdt;
    // vector<double> total_from; // accumulates total loss-per-unit-P from state [xi]
    // vector<double> total_gain; // accumulates total gain for state [xi]


    void sys(const state_type& s, state_type& sdot, const double t);
    bool hasRates = false; // flags whether Store has been populated yet.
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    void saveBound(const std::string& folder);
    state_type get_ground_state();


    bool good_state = true;
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
