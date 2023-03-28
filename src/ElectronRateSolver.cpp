/**
 * @file ElectronRateSolver.cpp
 * @author Alaric Sanders 
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

#include "ElectronRateSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "SplineIntegral.h"
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

/**
 * @brief unused
 */
void ElectronRateSolver::get_energy_bounds(double& max, double& min) {
    max = 0;
    min = 1e9;
    for(auto& atom : input_params.Store) {
        for(auto& r : atom.Photo) {
            if (r.energy > max) max=r.energy;
            if (r.energy < min) min=r.energy;
        }
        for(auto& r : atom.Auger) {
            if (r.energy > max) max=r.energy;
            if (r.energy < min) min=r.energy;
        }
    }
}

/// Finding regime boundaries. This is very rough, a gradient descent algorithm should be implemented.
// Should be cautious using in early times when wiggles are present.

/**
 * @brief Approximates midpoint between start_energy (e.g. a peak) and local density minimum.
 * @details  Finds density minimum by checking energies separated by del_energy until the energy increases for more than min_higher steps.
 * @param start_energy 
 * @param del_energy 
 * @param min don't seek past this point
 * @param max don't seek past this point
 */
double ElectronRateSolver::approx_nearest_min(size_t step, double start_energy,double del_energy, size_t min_sequential, double min, double max){
    assert(del_energy != 0);
    // Default values
    if (min < 0)
        min = 0;
    if(max < 0 || max > elec_grid_regions.bndry_E.back())
        max = elec_grid_regions.bndry_E.back();
    // Initialise
    double e = start_energy;
    double local_min = -1;
    double min_density = -1;
    double last_density = y[step].F(start_energy);
    size_t num_sequential = 0;
    if(min_sequential < 1) min_sequential = 1;
    // Seek local minimum.
    while (num_sequential < min_sequential+1){
        e += del_energy;
        if (e < min){
            local_min = min; break;}
        if (e > max){
            local_min = max; break;}

        double density = y[step].F(e);
        if(density > last_density){
            if (num_sequential == 0){
                local_min = e - del_energy;
                min_density = density;
            }
            num_sequential++;
        }
        if (density <= last_density || min_density <0){
            local_min = -1;
            num_sequential = 0;
            min_density = -1;
        }
        last_density = density; 
    }
    return local_min;
}

// num_sequential is num times the check passed.
double ElectronRateSolver::nearest_inflection(size_t step, double start_energy,double del_energy, size_t min_sequential, double min, double max){
    assert(del_energy != 0);
    // Default values
    if (min < 0)
        min = 0;
    if(max < 0 || max > elec_grid_regions.bndry_E.back())
        max = elec_grid_regions.bndry_E.back();
    // Initialise
    double e = start_energy;
    double inflection = -1;
    double min_density = -1;
    double last_density = y[step].F(start_energy);
    double last_grad = 0;
    size_t num_sequential = 0;    
    if(min_sequential < 1) min_sequential = 1;
    // Seek inflection.
    while (num_sequential < min_sequential+1){
        e += del_energy;
        if (e < min){
            inflection = min; break;}
        if (e > max){
            inflection = max; break;}
        double density = y[step].F(e);
        double mag_grad = abs((density - last_density)/del_energy);  // Should be fine unless we have extremely bad behaviour.
        if(mag_grad <= last_grad){
            if (num_sequential == 0){
                inflection = e - del_energy;
            }
            num_sequential++;
        }
        if (mag_grad >last_grad){
            inflection = -1;
            num_sequential = 0;
        }
        last_density = density; 
        last_grad = mag_grad;
    }
    return inflection;    
}

double ElectronRateSolver::approx_regime_bound(size_t step, double start_energy,double del_energy, size_t min_sequential, double min_distance, double min_inflection_fract, double _min, double _max){
    min_distance /= Constant::eV_per_Ha;
    // Find 0 of second derivative
    double inflection = nearest_inflection(step,start_energy,del_energy,min_sequential,_min,_max);
    // At min_distance, go double as 
    double A = 1/min_inflection_fract; // region between peak and inflection take up at least min_inflection_frac after min distance.
    double D = min_distance;
    int sign = (0 < del_energy) - (del_energy < 0);
    return sign*max(A*sqrt(abs(start_energy - inflection))*sqrt(D/A),D) + start_energy;

    // double local_min = approx_nearest_min(step,start_energy,del_energy,min_sequential,_min,_max);
    // // Return the midpoint as the boundary.
    // return (local_min + start_energy)/2;
}

double ElectronRateSolver::approx_regime_peak(size_t step, double lower_bound, double upper_bound, double del_energy){
    del_energy = abs(del_energy);
    assert(del_energy > 0);
    double e = lower_bound;
    double peak_density = 0;
    double peak_e = -1;
    // Seek local minimum.
    while (e < upper_bound){
        double density = y[step].F(e);
        if (density > peak_density){
            peak_density = density;
            peak_e = e;
        }
        e += del_energy; 
    }
    return peak_e;
}




/**
 * @brief Finds the energy bounds of the photoelectron region
 * @details Since it isn't obvious what function approximates the peak at a given time, 
 * we define a range that was found to experimentally give good results. 
 * @todo It may be better to have a couple of these regimes for the tallest few peaks. Wouldn't be too hard to implement.
 */
void ElectronRateSolver::dirac_energy_bounds(size_t step, double& max, double& min, double& peak, bool allow_shrinkage) {
    if(allow_shrinkage){
        max = -1;
        min = INFINITY;
    }

    double e_step_size = 50/Constant::eV_per_Ha;
    size_t num_sequential_needed = 3;

    double peak_density = -1e9;
    double min_photo_peak_considered = 3000/Constant::eV_per_Ha; // TODO instead ignore peaks that have negligible rates?
    //about halfway between the peak and the neighbouring local minima is sufficient for their mins and maxes
    for(auto& atom : input_params.Store) {
        for(auto& r : atom.Photo) {
            if (r.energy <  min_photo_peak_considered) continue;
            // Check if peak density
            double density = y[step].F(r.energy);
            if (density >= peak_density){
                peak = r.energy;
                peak_density = density; 
            }
            // Get bounds
            double lower_bound = approx_regime_bound(step,r.energy, -e_step_size, num_sequential_needed,1000,1./4.);
            double upper_bound = approx_regime_bound(step,r.energy, +e_step_size, num_sequential_needed,500,1./4.);//peak + 0.912*(peak - lower_bound);  // s.t. if lower_bound is 3/4 peak, upper_bound is 1.1*peak.
            if (upper_bound > max) max=upper_bound;
            if (lower_bound < min) min=lower_bound;
        }
    }
}

/**
 * @brief Finds the energy bounds of the MB within 2 deviations of the average energy.
 * @details 
 * @param step 
 * @param max 
 * @param min 
 * @param peak 
 */
void ElectronRateSolver::mb_energy_bounds(size_t step, double& _max, double& _min, double& peak, bool allow_shrinkage) {
    // Find e_peak = kT/2
    double min_energy = 0;
    double max_energy = 1000/Constant::eV_per_Ha;
    double e_step_size = 2/Constant::eV_per_Ha; //TODO maybe make this increase with increasing energy of previous peak.
    peak = approx_regime_peak(step,min_energy,max_energy,e_step_size);
    double kT = 2*peak;
    // CDF = Γ(3/2)γ(3/2,E/kT)
    double new_min = 0.2922*kT;  // 90% of electrons above this point
    if(_min < new_min || allow_shrinkage){
        _min = new_min;
    }
    size_t num_sequential_needed = 10; 
    double new_max = approx_regime_bound(step,peak, +e_step_size, num_sequential_needed,5,1./1.5); 
    //double new_max = 2.3208*kT; // 80% of electrons below this point (lower since not as sharp)
    if(_max < new_max || allow_shrinkage)
        _max = std::min(new_max,elec_grid_regions.bndry_E.back());
}

/**
 * @brief 
 *
 * @param _log 
 * @param init whether this is the start of the simulation and starting state must be initialised.  
 */

/**
 * @brief A basic quantifier of the ill-defined transition region
 * @details While this is a terrible approximation at later times,(it's unclear how to define the transition region as the peaks merge),
 * a) the dln(Λ)/dT < 1 as ln(Λ) propto ln(T_e^(3/2)), i.e. an accurate measure at low T is most important (kind of, at very low T transport coefficients dominate).
 * b) The coulomb logarithm is capped to 23.5 for ee collisions, which is where it is used. [and get_Q_ee limits it to approx. half of this cap though I need to determine why)
 * See https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary.
 * @param step 
 * @param g_min 
 */
void ElectronRateSolver::transition_energy(size_t step, double& g_min, bool allow_decrease){
    // We have instability in high energies at early times, so only consider high energies when there is no minimum 
    // at low energies. (This is very hacky.)
    double early_max = 1000;
    double early_min = approx_nearest_min(step,0,5,20,0,early_max);
    double new_min;
    if(early_min >= early_max){
        new_min = approx_nearest_min(step,0,10,10);
    }
    else{
        new_min = early_min;
    }
    if(allow_decrease) 
        g_min = new_min;
    else               
        g_min = max(g_min,new_min); // if the previous transition energy was higher, use that (good way to wait out instability).
}
 
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
            regimes.dirac_peak = 0;
            for(auto& atom : input_params.Store) {
                for(auto& r : atom.Photo) {     
                    if(r.energy > regimes.dirac_peak) 
                        regimes.dirac_peak = r.energy;
                }
            }
            regimes.dirac_min=regimes.dirac_peak*3/4,
            regimes.dirac_max=regimes.dirac_peak*1.14; //
        }
        else{
            dirac_energy_bounds(step,regimes.dirac_max,regimes.dirac_min,regimes.dirac_peak,true);
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

        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ New Knots ] -------------------\n" 
            "Time: "<<t[step]*Constant::fs_per_au <<" fs\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n" 
            <<"Photo [peak; range]: "<<regimes.dirac_peak*e<< "; "<< regimes.dirac_min*e<<" - "<<regimes.dirac_max*e<<"\n"
            <<"Transition energy: "<<param_cutoffs.transition_e*e<<""  
            << endl;
        }        
        else{}

    }


    Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);


    if (init){
        // Set up the rate equations (setup called from parent Adams_B  M)
        set_starting_state();
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
    const string banner = "================================================================================";
    cout<<banner<<endl;
    if (steps_per_time_update > 1){
        cout<< "Updating time every " << steps_per_time_update << " steps." <<"\n";
    }
    if (simulation_resume_time != simulation_start_time){
        cout<<"\033[33m"<<"Loaded simulation at:  "<<"\033[0m"<<(simulation_resume_time)*Constant::fs_per_au<<" fs"<<endl;
    }
    cout<<"\033[33m"<<"Final time step:  "<<"\033[0m"<<(simulation_end_time)*Constant::fs_per_au<<" fs"<<endl;
    cout<<banner<<endl;
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
        std::cout <<"[ sim ] Grid update period: "<<grid_update_period * Constant::fs_per_au<<" fs"<<std::endl;
    else 
        std::cout << "[ sim ] Using static grid" << std::endl;
    size_t steps_per_grid_transform =  round(input_params.Num_Time_Steps()*(grid_update_period/timespan_au));

    this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update,steps_per_grid_transform); // Inherited from ABM


    cout<<"[ Rate Solver ] Using timestep "<<this->dt*Constant::fs_per_au<<" fs"<<std::endl;
    
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
        std::cerr<<"\033[93;1m[ Rate Solver ] Halving timestep...\033[0m"<<std::endl;  // Unsure if this is faster than continuing the simulation with halved timesteps rather than restarting and doing so -S.P.
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
        set_starting_state();
        this->iterate(_log,simulation_start_time, simulation_end_time, simulation_resume_time, steps_per_time_update,steps_per_grid_transform); // Inherited from ABM
    }
    
    

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout << "[ Solver ] finished computation at " << ctime(&end_time) << endl;
    secs = elapsed_seconds.count();
    auto eii_m = std::chrono::duration_cast<std::chrono::minutes>(eii_time); 
    auto tbr_m = std::chrono::duration_cast<std::chrono::minutes>(tbr_time);
    auto ee_m = std::chrono::duration_cast<std::chrono::minutes>(ee_time);
    auto apply_delta_m = std::chrono::duration_cast<std::chrono::minutes>(apply_delta_time);
    auto misc_m = std::chrono::duration_cast<std::chrono::minutes>(misc_time);
    auto eii_s = std::chrono::duration_cast<std::chrono::seconds>(eii_time); 
    auto tbr_s = std::chrono::duration_cast<std::chrono::seconds>(tbr_time);
    auto ee_s = std::chrono::duration_cast<std::chrono::seconds>(ee_time);
    auto apply_delta_s = std::chrono::duration_cast<std::chrono::seconds>(apply_delta_time);
    auto misc_s = std::chrono::duration_cast<std::chrono::seconds>(misc_time);
    cout <<"[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << endl;
    cout <<"[ Solver ] get_Q_ee() took "<< ee_m.count() <<"m " << ee_s.count()%60 << "s" << endl;
    cout <<"[ Solver ] get_Q_eii() took "<< eii_m.count() <<"m " << eii_s.count()%60 << "s" << endl;
    cout <<"[ Solver ] get_Q_tbr() took "<< tbr_m.count() <<"m " << tbr_s.count()%60 << "s" << endl;
    cout <<"[ Solver ] apply_delta() took "<< apply_delta_m.count() <<"m " << apply_delta_s.count()%60 << "s" << endl;
    cout <<"[ Solver ] some misc processes took "<< misc_m.count() <<"m " << misc_s.count()%60 << "s" << endl;
    log_extra_details(_log);

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
void ElectronRateSolver::sys_bound(const state_type& s, state_type& sdot, const double t) {
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


    sdot.F.applyDelta(vec_dqdt);
    // This is loss.
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
    
    //if (!good_state) return;

    // compute the dfdt vector
    
    // for (size_t a = 0; a < s.atomP.size(); a++) {
    //     const bound_t& P = s.atomP[a];  
        
    // }
    
    

    #ifdef NO_EE
    #warning No electron-electron interactions
    #else
    auto t5 = std::chrono::high_resolution_clock::now();
    s.F.get_Q_ee(vec_dqdt, input_params.Plasma_Threads()); // Electron-electon repulsions
    auto t6 = std::chrono::high_resolution_clock::now();
    ee_time += t6 - t5;
    #endif
    auto t7 = std::chrono::high_resolution_clock::now();
    sdot.F.applyDelta(vec_dqdt);
    auto t8 = std::chrono::high_resolution_clock::now();
    apply_delta_time += t8 - t7;
}

//IOFunctions found in IOFunctions.cpp
