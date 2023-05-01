/**
 * @file ElectronRateSolver.cpp
 * @author Alaric Sanders & Spencer Passmore
 * @brief @copybrief ElectronRateSolver.h
 * @details @copydetail ElectronRateSolver.h
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
// (C) Spencer Passmore 2023


#include "ElectronRateSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "SplineIntegral.h"
#include "Display.h"
#include <fstream>
#include <algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <math.h>
#include <omp.h>
#include "config.h"


void ElectronRateSolver::set_starting_state(){
    
    if (input_params.Load_Folder() != ""){ 
        cout << "[ Plasma ] loading sim state from specified files." << endl;
        loadFreeRaw_and_times();
        loadBound();
        simulation_resume_time = t.back();
    }
    else{
        cout << "[ Plasma ] Creating ground state" << endl;
        this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, IVP_step_tolerance);
    }
}

state_type ElectronRateSolver::get_ground_state() {
    state_type initial_condition;
    assert(initial_condition.atomP.size() == input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        initial_condition.atomP[a][0] = input_params.Store[a].nAtoms;
        for(size_t i=1; i<initial_condition.atomP.size(); i++) {
            initial_condition.atomP[a][i] = 0.;
        }
    }
    initial_condition.F=0;
    // std::cout<<"[ Rate Solver ] initial condition:"<<std::endl;
    // std::cout<<"[ Rate Solver ] "<<initial_condition<<std::endl;
    return initial_condition;
}

// grid initilisation call order as of 3/04/23 TODO make clearer
// set_up_grid_and_compute_cross_sections(init = true) is called first
//  if static grid, use given static grid  
//  if dynamic grid, use default start
// Then call set_starting_state() (calls in set_up_grid_... now).
//  if loading a simulation:   
//      if static grid, transform loaded data to the given static grid
//      if dynamic grid, transform grid to the loaded data's grid/basis. 
//          Detect the transition energy
// If dynamic grid, call set_up_grid_and_compute_cross_sections(init = false) 
// after each non-zero integer multiple of the update period. 
void ElectronRateSolver::set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step, bool force_update) {//, bool recalc) {
    if (!init && input_params.elec_grid_type.mode != GridSpacing::dynamic){
        return; 
    }
    //
    
    bool recalc = true;
    if (init ){
        std::cout << "[ HF ] First computation of rates, logging HF iteration excesses. Future computations will not be logged. " << std::endl;
        input_params.calc_rates(_log, recalc);
    }
    else{
        ofstream dummy_log;
        input_params.calc_rates(dummy_log,recalc); 
    }
    hasRates = true;
    

    /// Set up the grid
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        double old_trans_e = param_cutoffs.transition_e;
        if(init){
            // Empirically-based guesses for good starts.
            std::cout << "[ Dynamic Grid ] Starting with initial grid guess"<<std::endl;
            _log << "[ Dynamic Grid ] Starting with initial grid guess"<<endl;
            double e = 1/Constant::eV_per_Ha;
            param_cutoffs.transition_e = 250*e;
            regimes.mb_peak=0; regimes.mb_min=4*e; regimes.mb_max=10*e;
            // cast a wide net to ensure we capture all dirac peaks. // TODO? there has to be a much better way than this.
            double photo_peak = -1, photo_min = 1e9, photo_max = -1e9; 
            for(auto& atom : input_params.Store) {
                for(auto& r : atom.Photo) {      
                    photo_min = min(photo_min,r.energy*5/6);
                    photo_max = max(photo_max,r.energy*7/6);
                }
            }
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                regimes.dirac_minimums[i] = photo_min + i*(photo_max-photo_min)/regimes.num_dirac_peaks;
                regimes.dirac_maximums[i] = photo_min + (i+1)*(photo_max-photo_min)/regimes.num_dirac_peaks;
            }
            Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions, input_params.elec_grid_preset);
        }
        else{ //TODO? reset dirac bounds on first one of these calls (since they are initialised to be very wide) IF shrinkage is to be turned off for dirac regions.
            double old_mb_min = regimes.mb_min, old_mb_max = regimes.mb_max;
            auto old_dirac_mins = regimes.dirac_minimums, old_dirac_maxs = regimes.dirac_maximums;
            
            double last_trans_e = param_cutoffs.transition_e;
            double dirac_peak_cutoff_density = y[step].F(last_trans_e)*last_trans_e*10;
            dirac_energy_bounds(step,regimes.dirac_maximums,regimes.dirac_minimums,regimes.dirac_peaks,true,regimes.num_dirac_peaks,dirac_peak_cutoff_density);
            mb_energy_bounds(step,regimes.mb_max,regimes.mb_min,regimes.mb_peak,false);
            // Check whether we should update
            if (regimes.mb_min == old_mb_min && regimes.mb_max == old_mb_max){
                if (force_update) _log <<"Forcing grid update"<< endl; // should be unlikely w/o adaptive b spline stuff? 
                else{
                    // We must be still at the start of the simulation, it's unnecessary, even detrimental, to update.
                    _log << "SKIPPED grid update as Maxwell-Boltzmann region's grid point range was the same"<<endl;
                    regimes.dirac_minimums = old_dirac_mins; 
                    regimes.dirac_maximums = old_dirac_maxs;
                    return;
                }
            }             
            transition_energy(step, param_cutoffs.transition_e);
            // More general basis refinement, but disabled because scope.            
            if (y[step].F.seek_basis(_log, input_params.elec_grid_type, step, param_cutoffs)){
                // failed, use regular method.
                Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
                //_log << "B-spline failed to find convergent basis, using backup dynamic grid."; _log.flush();   Disabled
            }
            else{
                 _log << "B-spline basis successfully converged"; _log.flush();
            }
        } 
        if (old_trans_e != param_cutoffs.transition_e){
            std::cout <<"thermal cutoff energy updated from:\n"
            <<old_trans_e*Constant::eV_per_Ha
            << "\nto:\n"
            << param_cutoffs.transition_e*Constant::eV_per_Ha<< std::endl;
        }
        else std::cout<<"Thermal cutoff energy had NO update." <<std::endl;    

    }
    else{
        _log << "[ Manual Grid ] Setting static grid" <<endl;
    }


    
    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);


    if (init){
        // Set up the rate equations (setup called from parent Adams_B  M)
        set_starting_state(); // possibly redundant steps in here.
    }

    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ New Knots ] -------------------\n" 
            "Time: "<<t[step]*Constant::fs_per_au <<" fs; "<<"Step: "<<step<<"\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<"\n"  
            << "Grid size: "<<Distribution::size<<"\n"
            << "-----------------------------------------------------" 
            << endl;
        }        
    }

    // create the tensor of coefficients TODO these aren't constant now - change to lowercase.
    RATE_EII.resize(input_params.Store.size());
    RATE_TBR.resize(input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        size_t N = y[step].F.num_basis_funcs();
        RATE_EII[a].resize(N);
        RATE_TBR[a].resize(N*(N+1)/2);
    }
    precompute_gamma_coeffs();
    Distribution::precompute_Q_coeffs(input_params.Store);
}

void ElectronRateSolver::set_grid_regions(ManualGridBoundaries gb){
    elec_grid_regions = gb;
}

void ElectronRateSolver::solve(ofstream & _log, const std::string& tmp_data_folder) {
    assert (hasRates || "No rates found! Use ElectronRateSolver::initialise_grid_with_computed_cross_sections(log)\n");
    auto start = std::chrono::system_clock::now();

    data_backup_folder = tmp_data_folder;

    if (simulation_resume_time > simulation_end_time){cout << "\033[91m[ Error ]\033[0m Simulation is attempting to resume from loaded data with a time after that which it is set to end at." << endl;}


    
    // generate text to display during simulation run
    std::stringstream plasma_header; 
    const string banner = "================================================================================";
    plasma_header<<banner<<"\n\r";
    if (input_params.time_update_gap > 0){
        plasma_header<< "Updating display every " <<input_params.time_update_gap*Constant::fs_per_au<<" fs."<<"\n\r";
    }
    if (simulation_resume_time != simulation_start_time){
        plasma_header/*<<"\033[33m"*/<<"Loaded simulation at:  "/*<<"\033[0m"*/<<(simulation_resume_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    }
    plasma_header/*<<"\033[33m"*/<<"Final time step:  "/*<<"\033[0m"*/<<(simulation_end_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
        plasma_header <<"[ sim ] Grid update period: "<<grid_update_period * Constant::fs_per_au<<" fs"<<"\n\r";
    else 
        plasma_header << "[ sim ] Using static grid" << "\n\r";

    plasma_header<<"[ Rate Solver ] Using initial timestep size of "<<this->dt*Constant::fs_per_au<<" fs"<<"\n\r";
    plasma_header<<banner<<"\n\r";

    steps_per_grid_transform =  round(input_params.Num_Time_Steps()*(grid_update_period/timespan_au));


    std::cout << plasma_header.str()<<std::flush; // display in regular terminal, so that it is still visible after end of program

    Display::header += plasma_header.str(); // display this in ncurses screen
    //Display::deactivate();

    // Call hybrid integrator to iterate through the time steps (good state)
    good_state = true;
    this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update); // Inherited from ABM

    
    
    // 12-04 leaving this for now since I'm using it to catch really bad errors, but it now serves no more than that in practice.
    // Using finer time steps to attempt to resolve NaN encountered in ODE solving. 
    /* Redoes the ENTIRE simulation. As it seems to be intended to, it should be fixed to start from timestep_reached. 
    // But really it should start from earlier than timestep_reached since bad states tend start before they appear, in a snowball-like effect. -S.P.
    */
    double time = this->t[0];
    int retries = 0;  // int retries = 1;  <---- Turned off for now
    while (!good_state) {
        if(_log.is_open()){
            cout << " Logging."<<endl;
            _log << endl << "[ Rate Solver ] "<< "Stopped at bad state encountered at t = " << timestep_reached << " fs"  << endl;
            _log.flush();
        }        
        std::cerr<<"\033[93;1m[ Rate Solver ] Halving timestep...\033[0m"<<std::endl;  // Would be ideal to go back e.g. 0.1 fs and increase time steps, except currently a fair few bad states caused by too few time steps get through detection.
        good_state = true;
        input_params.num_time_steps *= 2;           // todo separate from input params
        if (input_params.num_time_steps > MAX_T_PTS){
            std::cerr<<"\033[31;1m[ Rate Solver ] Exceeded maximum T size. Skipping remaining iteration."<<std::endl;
            break;
        }
        if (timestep_reached - time < 0 || retries == 0){
            std::cerr<<"\033[31;1m[ Rate Solver ] Shorter timestep failed to improve convergence. Skipping remaining iteration."<<std::endl;
            break;
        } else {  // I don't think this is working properly? Seems to restart entirely -S.P.
            time = timestep_reached;
        }
        retries--;
        set_starting_state(); //TODO this would be broken by dynamic grid
        this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update); // Inherited from ABM
    }

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout << "[ Solver ] finished computation at " << ctime(&end_time) << endl;
    secs = elapsed_seconds.count();
    
    // I'm sorry for this
    std::vector<std::chrono::duration<double, std::milli>> times{
    display_time, plot_time, dyn_dt_time, backup_time, pre_ode_time,
    dyn_grid_time, user_input_time, post_ode_time,
    pre_tbr_time, eii_time, tbr_time,
    ee_time, apply_delta_time
    };
    std::vector<std::string> tags{
    "display", "live plotting", "dt updates", "data backups", "pre_ode()",
    "dynamic grid updates", "user input detection", "post_ode()",
    "misc bound processes", "get_Q_eii()", "get_Q_tbr()",
    "get_Q_ee()", "apply_delta()"
    };
    
    stringstream solver_times;
    solver_times <<"\n[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << "\n";
    solver_times << its_dinner_time(times,tags);
    std::cout<<banner<<"\033[1;32m"<<solver_times.str()<<"\033[0m\n"<<endl;
    _log <<banner<<"\n"<<solver_times.str(); _log.flush();
}

string ElectronRateSolver::its_dinner_time(std::vector<std::chrono::duration<double, std::milli>> times, std::vector<std::string> tags){
    std::stringstream out;
    for(size_t i = 0; i < times.size(); i++){
        auto m = std::chrono::duration_cast<std::chrono::minutes>(times[i]); 
        auto s = std::chrono::duration_cast<std::chrono::seconds>(times[i]);
        out <<"[ Solver ] "<<tags[i]<<" took "<<m.count()<<"m "<<s.count()%60<< "s" << "\n";
    }
    return out.str();
}


// Populates RATE_EII and RATE_TBR based on calls to Distribution::Gamma_eii/ Distribution::Gamma_tbr
// 
void ElectronRateSolver::precompute_gamma_coeffs() {
    std::cout<<"[ Gamma precalc ] Beginning coefficient computation..."<<std::endl;
    size_t N = Distribution::size;
    for (size_t a = 0; a < input_params.Store.size(); a++) {
        std::cout<<"\n[ Gamma precalc ] Atom "<<a+1<<"/"<<input_params.Store.size()<<std::endl;
        auto eiiVec = input_params.Store[a].EIIparams;
        vector<RateData::InverseEIIdata> tbrVec = RateData::inverse(eiiVec);
        size_t counter=1;
        #pragma omp parallel default(none) shared(a, N, counter, tbrVec, eiiVec, RATE_EII, RATE_TBR, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t n=0; n<N; n++) {
                #pragma omp critical
                {
                    std::cout<<"\r[ Gamma precalc ] "<<counter<<"/"<<N<<" thread "<<omp_get_thread_num()<<std::flush;
                    counter++;
                }
                #ifndef NO_EII
                Distribution::Gamma_eii(RATE_EII[a][n], eiiVec, n);
                #endif
                // Weird indexing exploits symmetry of Gamma_TBR
                // such that only half of the coefficients are stored
                #ifndef NO_TBR
                for (size_t m=n+1; m<N; m++) {
                    size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                    // // k = N... N(N+1)/2
                    Distribution::Gamma_tbr(RATE_TBR[a][k], tbrVec, n, m);
                }
                Distribution::Gamma_tbr(RATE_TBR[a][n], tbrVec, n, n);
                #endif
                
            }
        }
    }
    std::cout<<"\n[ Gamma precalc ] Done."<<std::endl;
}


// The Big One: The global ODE function, separated into sys_bound and sys_ee.
// Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j] = d/dt(average-atomic-state)
// d/dt f = Q_B[f](t)                      = d/dt(electron-energy-distribution)  // TODO what is Q_B? Q_{Bound}? Q_{ee} contributes though. 

// Non-stiff part of the system. Bound-electron dynamics.
//void ElectronRateSolver::sys_bound(const state_type& s, state_type& sdot, state_type& s_bg ,const double t) {
void ElectronRateSolver::sys_bound(const state_type& s, state_type& sdot, state_type& s_bg ,const double t) {
    const int threads = input_params.Plasma_Threads();
    
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);

    if (!good_state) return;
    auto t9 = std::chrono::high_resolution_clock::now();
    for (size_t a = 0; a < s.atomP.size(); a++) {
        const bound_t& P = s.atomP[a];
        bound_t& Pdot = sdot.atomP[a];

        // PHOTOIONISATION
        double J = pf(t); // photon flux in atomic units
        for ( auto& r : input_params.Store[a].Photo) {
            double tmp = r.val*J*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*J*P[r.from]);
            sdot.bound_charge +=  tmp;
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*J*P[r.from]);
        }

        // FLUORESCENCE
        for ( auto& r : input_params.Store[a].Fluor) {
            Pdot[r.to] += r.val*P[r.from];
            Pdot[r.from] -= r.val*P[r.from];
            // Energy from optical photon assumed lost
        }

        // AUGER
        for ( auto& r : input_params.Store[a].Auger) {
            double tmp = r.val*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*P[r.from]);
            // sdot.F.add_maxwellian(r.energy*2./3., r.val*P[r.from]);
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*P[r.from]);
            sdot.bound_charge +=  tmp;
        }

        // EII / TBR bound-state dynamics
        double Pdot_subst [Pdot.size()] = {0};    // subst = substitute.
        double sdot_bound_charge_subst = 0; 
        size_t N = Distribution::size; 
        #pragma omp parallel for num_threads(threads) reduction(+ : Pdot_subst,sdot_bound_charge_subst)     
        for (size_t n=0; n<N; n++) {
            double tmp=0; // aggregator
            
            #ifndef NO_EII
            for (size_t init=0;  init<RATE_EII[a][n].size(); init++) {
                for (auto& finPair : RATE_EII[a][n][init]) {
                    tmp = finPair.val*s.F[n]*P[init];
                    Pdot_subst[finPair.idx] += tmp;
                    Pdot_subst[init] -= tmp;
                    sdot_bound_charge_subst += tmp;   //TODO Not a minus sign like others, double check that's intentional. -S.P.
                }
            }
            #endif
            
            
            #ifndef NO_TBR
            // exploit the symmetry: strange indexing engineered to only store the upper triangular part.
            // Note that RATE_TBR has the same geometry as InverseEIIdata.
            for (size_t m=n+1; m<N; m++) {
                size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                // k = N... N(N+1)/2-1
                // W += RATE_TBR[a][k]*s.F[n]*s.F[m]*2;
                for (size_t init=0;  init<RATE_TBR[a][k].size(); init++) {
                    for (auto& finPair : RATE_TBR[a][k][init]) {
                        tmp = finPair.val*s.F[n]*s.F[m]*P[init]*2;
                        Pdot_subst[finPair.idx] += tmp;
                        Pdot_subst[init] -= tmp;
                        sdot_bound_charge_subst -= tmp;
                    }
                }
            }
            // the diagonal
            // W += RATE_TBR[a][n]*s.F[n]*s.F[n];
            for (size_t init=0;  init<RATE_TBR[a][n].size(); init++) {
                for (auto& finPair : RATE_TBR[a][n][init]) {
                    tmp = finPair.val*s.F[n]*s.F[n]*P[init];
                    Pdot_subst[finPair.idx] += tmp;
                    Pdot_subst[init] -= tmp;
                    sdot_bound_charge_subst -= tmp;
                }
            }
            #endif
        }
        // Add parallel containers to their parent containers.
        for(size_t i=0;i < Pdot.size();i++){
            Pdot[i] += Pdot_subst[i];
        }
        sdot.bound_charge += sdot_bound_charge_subst;

        auto t10 = std::chrono::high_resolution_clock::now();
        pre_tbr_time += t10 - t9;

        // Free-electron parts
        #ifdef NO_EII
        #warning No impact ionisation
        #else
        auto t1 = std::chrono::high_resolution_clock::now();
        s.F.get_Q_eii(vec_dqdt, a, P, threads);
        auto t2 = std::chrono::high_resolution_clock::now();
        eii_time += t2 - t1;
        #endif
        #ifdef NO_TBR
        #warning No three-body recombination
        #else
        auto t3 = std::chrono::high_resolution_clock::now();
        s.F.get_Q_tbr(vec_dqdt, a, P, threads);  // Serially, this is the computational bulk of the program - S.P.
        auto t4 = std::chrono::high_resolution_clock::now();
        tbr_time += t4 - t3;
        #endif
    }

    // Add change to distribution
    sdot.F.applyDeltaF(vec_dqdt);
    // if(input_params.Filtration_File() == "")
    //     // No background, use geometric-dependent loss.
    //     sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);
    // else
    //     // background handles loss, apply photoelectron filtration with background
    //     sdot.F.addFiltration(s.F, s_bg.F,input_params.loss_geometry);
    sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);
    
    if (isnan(s.norm()) || isnan(sdot.norm())) {
        cerr<<"NaN encountered in ODE iteration."<<endl;
        cerr<< "t = "<<t*Constant::fs_per_au<<"fs"<<endl;
        good_state = false;
        timestep_reached = t*Constant::fs_per_au;
    }
}


// 'badly-behaved' i.e. stiff part of the system. Electron-electron interactions.
void ElectronRateSolver::sys_ee(const state_type& s, state_type& sdot, const double t) {
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);
    
    // compute the dfdt vector
    #ifdef NO_EE
    #warning No electron-electron interactions
    #else
    // Electron-electon repulsions
    auto t5 = std::chrono::high_resolution_clock::now();
    s.F.get_Q_ee(vec_dqdt, input_params.Plasma_Threads()); 
    auto t6 = std::chrono::high_resolution_clock::now();
    ee_time += t6 - t5;
    #endif
    // 
    // Add change to distribution
    auto t7 = std::chrono::high_resolution_clock::now();
    sdot.F.applyDeltaF(vec_dqdt);
    auto t8 = std::chrono::high_resolution_clock::now();
    apply_delta_time += t8 - t7;
}




//// Mid-simulation Non-ode functions

size_t ElectronRateSolver::load_checkpoint_and_decrease_dt(ofstream &_log, size_t current_n, Checkpoint _checkpoint){    
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    
    size_t n = _checkpoint.n; // index of step to load at 
    auto saved_knots = _checkpoint.knots;
    auto saved_states = _checkpoint.last_few_states;
    //std::vector<double> knots = _checkpoint.knots;
    std::vector<double> old_knots = Distribution::get_knot_energies();

    this->regimes = _checkpoint.regimes; // needed for when we next update the regimes (it needs the previous regime for reference)

    int remaining_steps = t.size() - (n+1);
    double fact = 1.5; // factor to increase time step density by (past the checkpoint). TODO need to implement max time steps.
    const string banner = "================================================================================";
    _log <<banner<<"\nEuler iterations exceeded beyond tolerable error at t="<<t[current_n]*Constant::fs_per_au
    <<". Decreasing remaining time step size from "<<this->dt*Constant::fs_per_au<<" fs to "<< dt/fact*Constant::fs_per_au<<" fs."<<endl;
    // REDUCE time step size by factor, then add on extra steps.
    this->dt/=fact;
    assert(remaining_steps > 0 && fact > 1);
    input_params.num_time_steps = t.size() + (fact - 1)*remaining_steps; // todo separate from input params
    steps_per_time_update = max(1 , (int)(input_params.time_update_gap/(timespan_au/input_params.num_time_steps))); 
    t.resize(n+1); t.resize(input_params.num_time_steps);
    // Load distributions (necessary in case we loaded [and transformed] the same point already) and also the previous few states so that we can use runge-kutta if needed.
    y.resize(n+1-saved_states.size());
    for(state_type state:saved_states){
        y.push_back(state);
    }
    y.resize(input_params.num_time_steps);    
    
    steps_per_grid_transform = round(steps_per_grid_transform*fact);

    // Set up the t grid       
    for (size_t i=n+1; i<input_params.num_time_steps; i++){   // TODO check potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
        this->t[i] = this->t[i-1] + this->dt;
    }

    // clear knot history past loaded step
    while(Distribution::knots_history.back().step > n){
        Distribution::knots_history.pop_back();
    }    
    // with the states loaded the spline factors are untouched, so now we just need to load the appropriate knots.
    assert(saved_knots == Distribution::load_knots_from_history(n));

    _log <<"Loaded knots..."<<endl;

    // Possibly temporary, but here we are just updating basis to ensure nothing breaks.
    // Ideally would just load the grid without affecting when grid updates.
    if (old_knots != Distribution::get_knot_energies()){
        // absolutely must update, no matter what!
        update_grid(_log,n,true);
    }
    else{
        update_grid(_log,n,false);
        // (disabled) // we should always update grid now anyway with the adaptive b splines 
        // update_grid(_log,n,true);
    }

    _log <<"Updated grid..."<<endl;

    // It is standard to use a method like fourth-order Runge-Kutta after changing the step size to get us going:
    // (TODO a lack of ee dynamics currently in this step is slightly off-putting as it is more important than when initialising)
    // We get indistinguishable results skipping this entirely (for a factor of 1.5 at least).
    /*
    n += this->order;
    for (size_t m = n-this->order; m < n; m++) {
            this->step_rk4(m);
    }
    std::cout.clear();
    */

    // Remember when time step was decreased, and flag to increase time step size again
    switch (input_params.pulse_shape){

    case PulseShape::gaussian:
        // TODO impact ionisation dominates so we would need a better method to do this. Arguably it's not worth it. 
        // if (t[n] <= 0){
        //     times_to_increase_dt.push_back(-t[n]+0.05/Constant::fs_per_au);
        // }
        break;
    case PulseShape::square:
        break;  
    }   

    _log << "Loaded checkpoint - resuming."<<endl;
    
    return n;
}

/**
 * @brief Increases time step size so as to reduce the number of remaining steps.
 * @param _log 
 * @param current_n 
 * @note current implementation means if a checkpoint was loaded multiple times (as can happen when we have multiple dt decreases within 2*checkpoint period,
 * then this will be called at the same time multiple times. Not critical to make it better though.  
 * An alternative would also be to make dt proportional to the incoming intensity. What be even better would be making it proportional to the
 * net magnitude of delta(density) / delta(t)
 */
void ElectronRateSolver::increase_dt(ofstream &_log, size_t current_n){
    // /* breaks things TODO may no longer break things since I found what was wrong... retry 0.5 fs FWHM mol file after fixing excessive computation time
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    
    size_t n = current_n;

    int remaining_steps = t.size() - (n+1);
    double fact = 1.5; // factor to increase time step size dt by.
    
    _log <<"Step size increase flagged at t="<<t[current_n]*Constant::fs_per_au
    <<". Increasing remaining time step size from "<<this->dt*Constant::fs_per_au<<" fs to "<< dt*fact*Constant::fs_per_au<<" fs."<<endl;
    
    // INCREASE time step size by factor, and reduce number of remaining steps.
    this->dt*=fact;
    assert(remaining_steps > 0 && fact > 1);
    input_params.num_time_steps = input_params.num_time_steps - (1-1/fact)*remaining_steps; // todo separate from input params
    t.resize(n+1); t.resize(input_params.num_time_steps);
    y.resize(n+1); y.resize(input_params.num_time_steps);
    steps_per_grid_transform = round(steps_per_grid_transform/fact);

    // Set up the t grid       
    for (int i=n+1; i<input_params.num_time_steps; i++){   // note potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
        this->t[i] = this->t[i-1] + this->dt;
    }
    // */
}
//IOFunctions found in IOFunctions.cpp


void ElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){
    auto t_start = std::chrono::high_resolution_clock::now();
    

    ////// Display info ////// (only the regular stuff seen each step, not "popups" from popup_stream)
    auto t_start_disp = std::chrono::high_resolution_clock::now();
    if ((n-this->order)%steps_per_time_update == 0){
        Display::display_stream.str(Display::header); // clear display string
        Display::display_stream<< "\n\r"
        << "--- Press BACKSPACE/DEL to end simulation and save the data ---\n\r"   
        << "[ sim ] Next data backup in "<<(minutes_per_save - std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save)).count()<<" minute(s).\n\r"  
        << "[ sim ] Current timestep size = "<<this->dt*Constant::fs_per_au<<" fs\n\r"   
        << "[ sim ] t="
        << std::left<<std::setfill(' ')<<std::setw(6)
        << this->t[n] * Constant::fs_per_au << " fs\n\r" 
        << "[ sim ] " <<Distribution::size << " knots currently active\n\r";
        //<< Distribution::get_knot_energies() << "\n\r"; 
        // << flush; 
        Display::show(Display::display_stream);
    }        
    display_time += std::chrono::high_resolution_clock::now() - t_start_disp;  

    ////// live plotting ////// 
    auto t_start_plot = std::chrono::high_resolution_clock::now();
    if ((n-this->order+1)%100 == 0){ // TODO implement minimum time
        size_t num_pts = 4000;
        py_plotter.plot_frame(
            Distribution::get_energies_eV(num_pts),
            this->y[n].F.get_densities(num_pts,Distribution::get_knot_energies()), 
            Distribution::get_trimmed_knots(Distribution::get_knot_energies())
        );
    }        
    plot_time += std::chrono::high_resolution_clock::now() - t_start_plot;  

    ////// Dynamic time updates //////
    auto t_start_dt = std::chrono::high_resolution_clock::now();
    // store checkpoint
    if ((n-this->order+1)%checkpoint_gap == 0){  // 
        old_checkpoint = checkpoint;
        ptrdiff_t end_vect_idx = &y[n] - &y[0];  
        ptrdiff_t start_vect_idx = &y[n-order] - &y[0];  
        auto check_states = std::vector<int>(start_vect_idx, end_vect_idx+1);
        checkpoint = {n, Distribution::get_knot_energies(),this->regimes};
    }
    if (euler_exceeded || !good_state){     
        bool store_dt_increase = true;   
        if (euler_exceeded){    
            if (!increasing_dt){
                Display::popup_stream <<"\nMax euler iterations exceeded, reloading checkpoint at "<<this->t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
            }
            else{
            Display::popup_stream << "\nReloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au
            <<" fs, as increasing step size led to euler iterations being exceeded. Will attempt again after "<<2*checkpoint_gap<< "steps.\n\r";
            times_to_increase_dt.back() = t[n] + 2*checkpoint_gap*dt;  
            store_dt_increase = false;
            increasing_dt = 0;                
            }
        }
        else{
            Display::popup_stream <<"\n NaN encountered in ODE iteration. Reloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
        }
        Display::show(Display::display_stream,Display::popup_stream);
        // Reload at checkpoint's n, updating the n step, and decreasing time step length
        n = load_checkpoint_and_decrease_dt(_log,n,old_checkpoint);  // virtual function overridden by ElectronRateSolver
        if (!store_dt_increase){
            // remove the extra time added.   
            times_to_increase_dt.pop_back();           
        }
        good_state = true;
        euler_exceeded = false;
        checkpoint = old_checkpoint;
    }
    else if (!times_to_increase_dt.empty() && times_to_increase_dt.back() < t[n]){
        Display::popup_stream <<"\nAttempting to increase dt\n\r";
        Display::show(Display::display_stream,Display::popup_stream);
        if (increasing_dt > 0){
            increasing_dt--;
        
            if (increasing_dt == 0){
                //successfully increased step size
                times_to_increase_dt.pop_back();
            }
        }
        else{
            // Attempt to increase step size (i.e. decrease num remaining steps)
            this->increase_dt(_log,n);
            // Flag that we are testing to see if still converging.
            size_t num_steps_to_check = 20;
            assert(num_steps_to_check < checkpoint_gap);
            increasing_dt = num_steps_to_check;
        }
    }
    dyn_dt_time += std::chrono::high_resolution_clock::now() - t_start_dt;
    
    
    /// Save data periodically in case of crash (currently need to manually copy mol file from log folder and move folder to output if want to plot). 
    if(std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save) > minutes_per_save){
        auto t_start_backup = std::chrono::high_resolution_clock::now();
          _log << "Saving output backup" << endl;
          Display::popup_stream <<"\nSaving output backup, simulation will resume soon...\n\r";
          Display::show(Display::display_stream,Display::popup_stream);        
          time_of_last_save = std::chrono::high_resolution_clock::now();
          size_t old_size = y.size();
          y.resize(n+1); t.resize(n+1);
          save(data_backup_folder);
          y.resize(old_size); t.resize(old_size);
          // re-set the future t values       
          for (int i=n+1; i<t.size()-1; i++){   // note potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
              this->t[i] = this->t[i-1] + this->dt;
          }          
        backup_time += std::chrono::high_resolution_clock::now() - t_start_backup;
    }
    
    
    pre_ode_time += std::chrono::high_resolution_clock::now() - t_start;  
}

int ElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){
    auto t_start = std::chrono::high_resolution_clock::now();

    //////  Dynamic grid updater ////// 
    auto t_start_grid = std::chrono::high_resolution_clock::now();
    if ((n-this->order+1)%steps_per_grid_transform == 0){ // TODO would be good to have a variable that this is equal to that is modified to account for changes in time step size. If a dt decreases you push back the grid update. If you increase dt you could miss it.
        Display::popup_stream << "\nUpdating grid... \n\r"; 
        _log << "[ Dynamic Grid ] Updating grid" << endl;
        Display::show(Display::display_stream,Display::popup_stream);   
        update_grid(_log,n+1,false);       
    }   
    // move from initial grid to dynamic grid shortly after a fresh simulation's start.
    else if (n-this->order == max(2,(int)(steps_per_grid_transform/10)) && (input_params.Load_Folder() == "") && !grid_initialised){
        Display::popup_stream << "\n Performing initial grid update... \n\r"; 
        _log << "[ Dynamic Grid ] Performing initial grid update..." << endl;
        Display::show(Display::display_stream,Display::popup_stream); 
        update_grid(_log,n+1,true); 
        //////// restart simulation with better grid.
        /* bugged but I'll fix if I get time.
        n = order;
        old_checkpoint.regimes = regimes;
        old_checkpoint.knots = Distribution::get_knot_energies();
        old_checkpoint.n = n;
        checkpoint = old_checkpoint;
        while (Distribution::knots_history.size() > 0){
            Distribution::knots_history.pop_back();
        }
        Distribution::set_basis(n, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
        y.resize(n+1); y.resize(t.size());
        t.resize(n+1); t.resize(y.size());
        // Set up the t grid       
        for (int i=n+1; i<input_params.num_time_steps; i++){   // note potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
            this->t[i] = this->t[i-1] + this->dt;
        }
        // same as from update_grid()
        std::vector<double> new_energies = Distribution::get_knot_energies();
        zero_y = get_ground_state();        
        for (size_t m = n+1 - this->order; m < n+1; m++) {
            assert(this-> order < steps_per_grid_transform);
            Distribution::load_knots_from_history(n-1); // the n - 1 is correct. 
            y[m].F.transform_basis(new_energies);
        }          
        grid_initialised = true;
        */
        ////////
    }
    dyn_grid_time += std::chrono::high_resolution_clock::now() - t_start_grid;  
    
    //////  Check if user wants to end simulation early ////// 
    auto t_start_usr = std::chrono::high_resolution_clock::now();
    auto ch = wgetch(Display::win);
    if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){   
        flushinp(); // flush buffered inputs
        Display::popup_stream <<"\n\rExiting early... press backspace/del again to confirm or any other key to cancel and resume the simulation \n\r";
        Display::show(Display::display_stream,Display::popup_stream);
        nodelay(Display::win, false);
        ch = wgetch(Display::win);  // note implicitly refreshes screen
        if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){
            y.resize(n);
            t.resize(n);    
            return 1;
        }
        nodelay(Display::win, true);
        Display::show(Display::display_stream);
    }     
    
    user_input_time += std::chrono::high_resolution_clock::now() - t_start_usr;  
    
    if (Display::popup_stream.rdbuf()->in_avail() != 0){
        // clear stream
        Display::popup_stream.str(std::string());
    }   
    
    
    post_ode_time += std::chrono::high_resolution_clock::now() - t_start;    
    return 0;
}

void ElectronRateSolver::update_grid(ofstream& _log, size_t latest_step, bool force_update){
    size_t n = latest_step;
    //// Set up new grid ////        
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    // Latest step is n, so we decide our new grid based on that step, then transform N = "order" of the prior points to the new basis.
    set_up_grid_and_compute_cross_sections(_log,false,n,force_update); // virtual function overridden by ElectronRateSolver
    std::cout.clear();
    //// We need some previous points needed to perform next ode step, so transform them to new basis ////
    std::vector<double> new_energies = Distribution::get_knot_energies();
    // zero_y is used as the starting state for each step and represents the ground state, so we need to reset it so it has the right knots.
    zero_y = get_ground_state();
    cout << endl;            
    for (size_t m = n+1 - this->order; m < n+1; m++) {
        // We transform this step to the correct basis, but we also need a few steps to get us going, 
        // so we transform a few previous steps 
        // Kinda goofy but it's necessary due to the static variables. 
        assert(this-> order < steps_per_grid_transform);
        Distribution::load_knots_from_history(n-1); // the n - 1 is correct. 
        y[m].F.transform_basis(new_energies);
    }  
    // The next containers are made to have the correct size, as the initial state is set to tmp=zero_y and sdot is set to an empty state. 
} 