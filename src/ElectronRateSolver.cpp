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
    this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, 5e-3);
    if (input_params.Load_Folder() != ""){ 
        cout << "[ Plasma ] loading sim state from specified files." << endl;
        loadFreeRaw_and_times();
        loadBound();
        simulation_resume_time = t.back();
    }
    else{
        cout << "[ Plasma ] Creating ground state" << endl;
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
void ElectronRateSolver::set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step) {//, bool recalc) {
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
            std::cout << "[ Dynamic Grid ] Starting with initial grid" << std::endl;
            double e = 1/Constant::eV_per_Ha;
            param_cutoffs.transition_e = 250*e;
            regimes.mb_peak=0; regimes.mb_min=4*e; regimes.mb_max=10*e;
            // cast a wide net to ensure we capture all dirac peaks. // TODO there has to be a much better way than this.
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
        }
        else{ //TODO reset dirac bounds on first one of these calls (since they are initialised to be very wide) IF shrinkage is to be turned off for dirac regions.
            double last_trans_e = param_cutoffs.transition_e;
            double dirac_peak_cutoff_density = y[step].F(last_trans_e)*last_trans_e*10;
            dirac_energy_bounds(step,regimes.dirac_maximums,regimes.dirac_minimums,regimes.dirac_peaks,true,regimes.num_dirac_peaks,dirac_peak_cutoff_density);
            mb_energy_bounds(step,regimes.mb_max,regimes.mb_min,regimes.mb_peak);
            transition_energy(step, param_cutoffs.transition_e);
        } 
        if (old_trans_e != param_cutoffs.transition_e){
            std::cout <<"thermal cutoff energy updated from:\n"
            <<old_trans_e*Constant::eV_per_Ha
            << "\nto:\n"
            << param_cutoffs.transition_e*Constant::eV_per_Ha<< std::endl;
        }
        else std::cout<<"Thermal cutoff energy had NO update." <<std::endl;    

    }


    Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);


    if (init){
        // Set up the rate equations (setup called from parent Adams_B  M)
        set_starting_state();
    }

    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ New Knots ] -------------------\n" 
            "Time: "<<t[step]*Constant::fs_per_au <<" fs\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<""  
            << endl;
        }        
        else{}
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

void ElectronRateSolver::set_grid_regions(GridBoundaries gb){
    elec_grid_regions = gb;
}

void ElectronRateSolver::solve(ofstream & _log) {
    assert (hasRates || "No rates found! Use ElectronRateSolver::initialise_grid_with_computed_cross_sections(log)\n");
    auto start = std::chrono::system_clock::now();

    if (simulation_resume_time > simulation_end_time){cout << "\033[91m[ Error ]\033[0m Simulation is attempting to resume from loaded data with a time after that which it is set to end at." << endl;}
    
    int steps_per_time_update = max(1 , (int)(input_params.time_update_gap/(timespan_au/input_params.num_time_steps))); 


    

    // Call hybrid integrator to iterate through the time steps (good state)
    good_state = true;
    std::stringstream plasma_header; 
    const string banner = "================================================================================";
    plasma_header<<banner<<"\n\r";
    if (steps_per_time_update > 1){
        plasma_header<< "Updating time every " << steps_per_time_update << " steps." <<"\n\r";
    }
    if (simulation_resume_time != simulation_start_time){
        plasma_header<<"\033[33m"<<"Loaded simulation at:  "<<"\033[0m"<<(simulation_resume_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    }
    plasma_header<<"\033[33m"<<"Final time step:  "<<"\033[0m"<<(simulation_end_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    plasma_header<<banner<<"\n\r";
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
        plasma_header <<"[ sim ] Grid update period: "<<grid_update_period * Constant::fs_per_au<<" fs"<<"\n\r";
    else 
        plasma_header << "[ sim ] Using static grid" << "\n\r";

    plasma_header<<"[ Rate Solver ] Using timestep "<<this->dt*Constant::fs_per_au<<" fs"<<"\n\r";
    size_t steps_per_grid_transform =  round(input_params.Num_Time_Steps()*(grid_update_period/timespan_au));


    std::cout << plasma_header.str()<<std::flush; // display in regular terminal, so that it is still visible after end of program

    Display::header += plasma_header.str(); // display this in ncurses screen
    //Display::deactivate();

    this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update,steps_per_grid_transform); // Inherited from ABM

    
    
    
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
        input_params.num_time_steps *= 2;           
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
        this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update,steps_per_grid_transform); // Inherited from ABM
    }
    
    
    // Times
    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout << "[ Solver ] finished computation at " << ctime(&end_time) << endl;
    secs = elapsed_seconds.count();
    //mins
    auto eii_m = std::chrono::duration_cast<std::chrono::minutes>(eii_time); 
    auto tbr_m = std::chrono::duration_cast<std::chrono::minutes>(tbr_time);
    auto ee_m = std::chrono::duration_cast<std::chrono::minutes>(ee_time);
    auto apply_delta_m = std::chrono::duration_cast<std::chrono::minutes>(apply_delta_time);
    auto misc_m = std::chrono::duration_cast<std::chrono::minutes>(misc_time);
    //secs
    auto eii_s = std::chrono::duration_cast<std::chrono::seconds>(eii_time); 
    auto tbr_s = std::chrono::duration_cast<std::chrono::seconds>(tbr_time);
    auto ee_s = std::chrono::duration_cast<std::chrono::seconds>(ee_time);
    auto apply_delta_s = std::chrono::duration_cast<std::chrono::seconds>(apply_delta_time);
    auto misc_s = std::chrono::duration_cast<std::chrono::seconds>(misc_time);
    stringstream solver_times;
    solver_times <<"[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << "\n";
    solver_times <<"[ Solver ] get_Q_ee() took "<< ee_m.count() <<"m " << ee_s.count()%60 << "s" << "\n";
    solver_times <<"[ Solver ] get_Q_eii() took "<< eii_m.count() <<"m " << eii_s.count()%60 << "s" << "\n";
    solver_times <<"[ Solver ] get_Q_tbr() took "<< tbr_m.count() <<"m " << tbr_s.count()%60 << "s" << "\n";
    solver_times <<"[ Solver ] apply_delta() took "<< apply_delta_m.count() <<"m " << apply_delta_s.count()%60 << "s" << "\n";
    solver_times <<"[ Solver ] some misc processes took "<< misc_m.count() <<"m " << misc_s.count()%60 << "s";
    std::cout<<banner<<"\n"<<"\033[1;32m"<<solver_times.str()<<"\033[0m\n"<<endl;
    _log <<banner<<"\n"<<solver_times.str(); _log.flush();

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
        auto t9 = std::chrono::high_resolution_clock::now();
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
        misc_time += t10 - t9;

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

//IOFunctions found in IOFunctions.cpp
