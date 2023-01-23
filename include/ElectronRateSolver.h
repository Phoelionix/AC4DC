/**
 * @file ElectronRateSolver.h
 * @author Alaric Sanders 
 * @brief 
 * @details Core part of Sanders' continuum plasma extension. This was historically named ElectronSolver.
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
    void solve();
    void save(const std::string& folder);
    void compute_cross_sections(std::ofstream& _log, bool recalc=true);
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
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    void saveBound(const std::string& folder);
    state_type get_ground_state();


    bool good_state = true;
    double timestep_reached = 0;
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
