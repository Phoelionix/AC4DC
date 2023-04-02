/**
 * @file ElectronRateSolver.h
 * @authors Alaric Sanders & Spencer Passmore 
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
        
        if (input_params.pulse_shape ==  PulseShape::square){  //-FWHM <= t <= 3FWHM
            simulation_start_time = -timespan_au/4;
            simulation_end_time =  3*timespan_au/4; 
        } else { //-2FWHM <= t <= 2FWHM
            simulation_start_time = -timespan_au/2;
            simulation_end_time = timespan_au/2; 
        }
        if (input_params.Cutoff_Inputted()){
            simulation_end_time = input_params.Simulation_Cutoff(); // This does affect the fineness of the output
        }
        simulation_resume_time = simulation_start_time;
        set_grid_regions(input_params.elec_grid_regions);

        grid_update_period = input_params.Grid_Update_Period();

        if(input_params.Filtration_File() != ""){
            load_filtration_file();
        }
        
    }
    /// Solve the rate equations
    void solve(ofstream & _log);
    void save(const std::string& folder);
    /// Sets up the rate equations, which requires computing the atomic cross-sections/avg. transition rates to get the coefficients.
    void set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step = 0); //bool recalc=true);
    void tokenise(std::string str, std::vector<double> &out, const char delim = ' ');

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
    GridBoundaries elec_grid_regions;
    FeatureRegimes regimes;
    Cutoffs param_cutoffs; 
    Pulse pf;
    double timespan_au; // Atomic units
    double simulation_start_time;  // [Au]
    double simulation_resume_time; // [Au] same as simulation_start_time unless loading simulation state.
    double simulation_end_time;  // [Au]    
    double fraction_of_pulse_simulated;
    double grid_update_period; // time period between dynamic grid updates.

    void load_filtration_file(){};
    // Model parameters

    // arrays computed at class initialisation
    vector<vector<eiiGraph> > RATE_EII;
    vector<vector<eiiGraph> > RATE_TBR;

    // Dynamic grid
    void dirac_energy_bounds(size_t step, std::vector<double>& maximums, std::vector<double>& minimums, std::vector<double>& peaks, bool allow_shrinkage,size_t num_peaks, double peak_min_density);
    void mb_energy_bounds(size_t step, double& max, double& min, double& peak_density, bool allow_shrinkage = false);
    void transition_energy(size_t step, double& g_min);
    double approx_nearest_peak(size_t step, double start_energy,double del_energy, size_t min_sequential, double min = -1, double max =-1);  
    double approx_regime_trough(size_t step, double lower_bound, double upper_bound, double del_energy);
    double nearest_inflection(size_t step, double start_energy,double del_energy, size_t min_sequential, double min = -1, double max =-1);  
    double approx_regime_bound(size_t step, double start_energy,double del_energy, size_t min_sequential, double min_distance = 40, double min_inflection_fract = 1./4., double min = -1, double max=-1);  
    double approx_regime_peak(size_t step, double lower_bound, double upper_bound, double del_energy);  
    std::vector<double> approx_regime_peaks(size_t step, double lower_bound, double upper_bound, double del_energy, size_t num_peaks = 1, double min_density = 0);
    double approx_regime_trough(size_t step, double lower_bound, double upper_bound, double del_energy,size_t min_sequential);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_initial_conditions();

    /////// Overrides virtual system state methods
    void sys_bound(const state_type& s, state_type& sdot, state_type& s_bg, const double t); // general dynamics (uses explicit mehtod)
    void sys_ee(const state_type& s, state_type& sdot, const double t); // electron-electron (uses implicit method)
    /////// 

    bool hasRates = false; // flags whether Store has been populated yet.
    void copyInput(const std::string& src,const std::string& dir);
    /// Saves a table of free-electron dynamics to file fname
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    /// Loads the table of free-electron dynamics at the given time
    void loadFreeRaw_and_times();
    void loadBound();
    /// saves a table of bound-electron dynamics , split by atom, to folder dir.
    void saveBound(const std::string& folder);
    void saveBoundRaw(const std::string& folder);

    void set_grid_regions(GridBoundaries gb);
    void set_starting_state();
    state_type get_ground_state();

    //void high_energy_stability_check();
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
