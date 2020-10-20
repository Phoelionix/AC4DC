#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

/*
This file should arguably be called RateEquationSOlver, however, for historical reasons it is not.
*/
// #include <boost/numeric/odeint.hpp>
#include "HybridIntegrator.hpp"
#include "RateSystem.h"
#include "Constant.h"
#include "MolInp.h"
#include "Input.h"
#include "Pulse.h"
#include <iostream>
#include <stdexcept>

// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

#define MAX_T_PTS 1e6

class ElectronSolver : private ode::Hybrid<state_type>
{
public:
    ElectronSolver(const char* filename, ofstream& log) :
    Hybrid(3), input_params(filename, log), pf() // (order Adams method)
    {
        pf.set_shape(input_params.pulse_shape);
        pf.set_pulse(input_params.Fluence(), input_params.Width());
        timespan_au = input_params.Width()*5;
    }
    void solve();
    void save(const std::string& folder);
    void compute_cross_sections(std::ofstream& _log, bool recalc=true);
private:
    MolInp input_params;
    Pulse pf;
    double timespan_au; // Atomic units
    // Model parameters
    

    // arrays computed at class initialisation
    vector<vector<eiiGraph> > RATE_EII;
    vector<vector<eiiGraph> > RATE_TBR;

    void get_energy_bounds(double& max, double& min);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_initial_conditions();

    void sys(const state_type& s, state_type& sdot, const double t); // general dynamics (uses explicit mehtod)
    void sys2(const state_type& s, state_type& sdot, const double t); // electron-electron (uses implicit method)
    bool hasRates = false; // flags whether Store has been populated yet.
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    void saveBound(const std::string& folder);
    state_type get_ground_state();


    bool good_state = true;
    double timestep_reached = 0;
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
