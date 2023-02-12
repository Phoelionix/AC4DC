/**
 * @file ElectronRateSolver.h
 * @author Alaric Sanders 
 * @brief Defines the ElectronRateSolver class, which executes the high-level operations of the electron ODE solving.
 * @details This is the central hub of Sanders' continuum plasma extension, and is called by main.cpp. This was historically named ElectronSolver.
 * @todo The cross-sectional computations should be computed before the class is called, in main.cpp, rather than within the class, decoupling the parts of the code.
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
// (C) Alaric Sanders 2020

#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

// #include <boost/numeric/odeint.hpp>
#include "HybridIntegrator.hpp"
#include "RateSystem.h"
#include "Constant.h"
#include "MolInp.h"
#include "Input.h"
#include "Pulse.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <chrono>

// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

#define MAX_T_PTS 1e6

class ElectronRateSolver : private ode::Hybrid<state_type>
{
public:
    ElectronRateSolver(const char* filename, ofstream& log) :
    Hybrid(3), input_params(filename, log), pf() // (order Adams method)
    {
        pf.set_shape(input_params.pulse_shape);
        pf.set_pulse(input_params.Fluence(), input_params.Width());
        timespan_au = input_params.Width()*4;   // 4*FWHM, capturing effectively entire pulse.
        // Get cutoff
            truncated_timespan = timespan_au*input_params.Simulated_Fraction();
    }
    /// Solve the rate equations
    void solve(ofstream & _log);
    void save(const std::string& folder);
    /// Sets up the rate equations, which requires computing the atomic cross-sections/avg. transition rates to get the coefficients.
    void compute_cross_sections(std::ofstream& _log, bool recalc=true);
    void set_load_params(pair<string,double> name_time){load_fname = name_time.first; latest_start_time = name_time.second;}

    /// Number of secs taken for simulation to run
    long secs;

    /// time duration of get_Q_tbr and get_Q_eii
    std::chrono::duration<double, std::milli> tbr_time;
    std::chrono::duration<double, std::milli> eii_time;
    std::chrono::duration<double, std::milli> ee_time;
    std::chrono::duration<double, std::milli> apply_delta_time;
    std::chrono::duration<double, std::milli> misc_time;
private:
    MolInp input_params;  // (Note this is initialised/constructed in the above constructor)
    Pulse pf;
    double timespan_au; // Atomic units
    double truncated_timespan;
    double fraction_of_pulse_simulated;
    // Model parameters
    

    // arrays computed at class initialisation
    vector<vector<eiiGraph> > RATE_EII;
    vector<vector<eiiGraph> > RATE_TBR;

    void get_energy_bounds(double& max, double& min);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_initial_conditions();

    /////// Overrides virtual system state methods
    void sys_bound(const state_type& s, state_type& sdot, const double t); // general dynamics (uses explicit mehtod)
    void sys_ee(const state_type& s, state_type& sdot, const double t); // electron-electron (uses implicit method)
    /////// 

    bool hasRates = false; // flags whether Store has been populated yet.
    void copyInput(const std::string& src,const std::string& dir);
    /// Saves a table of free-electron dynamics to file fname
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    /// Loads the table at the given time
    void loadFreeRaw();
    /// saves a table of bound-electron dynamics , split by atom, to folder dir.
    void saveBound(const std::string& folder);
    /// Log final details pertaining to the simulation's execution to file fname (e.g. total runtime)
    void log_extra_details(ofstream & _log);

    state_type get_starting_state();
    state_type load_state();
    state_type get_ground_state();

    string load_fname = "";  // if "" don't load anything.
    double latest_start_time;  // [fs]
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
