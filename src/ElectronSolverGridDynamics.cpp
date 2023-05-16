/**
 * @file ElectronSolverGridDynamics.cpp
 * @author Spencer Passmore
 * @brief 
 * @note density refers to electron energy density.
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
// (C) Spencer Passmore 2023

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
double ElectronRateSolver::approx_nearest_peak(size_t step, double start_energy,double del_energy, size_t min_sequential, double min, double max){
    assert(del_energy != 0);
    // Default values
    if (min < 0)
        min = 0;
    if(max < 0 || max > Distribution::get_max_E())
        max = Distribution::get_max_E();
    // Initialise
    double e = start_energy;
    double local_max = -1;
    double max_density = -1;
    double last_density = y[step].F(start_energy)*start_energy;
    size_t ascended_count = 0; // make sure we go up. We require it goes up min_sequential times to give a correlation with this param designed to address instability.
    size_t num_sequential = 0;
    if(min_sequential < 1) min_sequential = 1;
    // Seek local maximum.
    while (num_sequential < min_sequential+1){
        e += del_energy;
        if (e < min){
            local_max = min; break;}
        if (e > max){
            local_max = max; break;}

        double density = y[step].F(e)*e;  //  electron energy density.
        if(density < last_density && ascended_count >= min_sequential){
            if (num_sequential == 0){
                local_max = e - del_energy;
                max_density = density;
            }
            num_sequential++;
        }
        if (density >= last_density || max_density < 0){
            local_max = -1;
            num_sequential = 0;
            max_density = -1;
        }
        if (density >= last_density){
            ascended_count++;        
        }
        last_density = density; 
    }
    return local_max;
}

// num_sequential is num times the check passed.
// TODO instead of using min_sequential, smooth out distribution over multiple knots.
// for checks starting from a peak, could make it skip a fake inflection by checking moving towards min and min is negative, and then skipping the following inflection too.
double ElectronRateSolver::nearest_inflection(size_t step, double start_energy,double del_energy, size_t min_sequential, double min, double max){
    assert(del_energy != 0);
    // Default values
    if (min < 0)
        min = 0;
    if(max < 0 || max > Distribution::get_max_E())
        max = Distribution::get_max_E();
    // Initialise
    double e = start_energy;
    double inflection = -1;
    double last_density = y[step].F(start_energy)*start_energy;
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
        double density = y[step].F(e)*e; // electron energy density
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

// min distance in eV
double ElectronRateSolver::approx_regime_bound(size_t step, double start_energy,double del_energy, size_t min_sequential, double min_distance, double min_inflection_fract, double _min, double _max){
    min_distance /= Constant::eV_per_Ha;
    // Find 0 of second derivative
    double inflection = nearest_inflection(step,start_energy,del_energy,min_sequential,_min,_max);
    std::cout << inflection*Constant::eV_per_Ha;
    // TODO For dirac use a function that has this changed such that, from the inflection above, 
    // it finds the next minimum OR the next point at the cutoff energy, whichever comes first. 
    double A = 1/min_inflection_fract; // region between peak and inflection take up min_inflection_frac at min distance.
    double D = min_distance;
    int sign = (0 < del_energy) - (del_energy < 0);
    //return sign*max(A*sqrt(abs(start_energy - inflection))*sqrt(D/A),D) + start_energy;
    return sign*max(A*abs(start_energy - inflection),D) + start_energy;
}

double ElectronRateSolver::approx_regime_peak(size_t step, double lower_bound, double upper_bound, double del_energy){
    del_energy = abs(del_energy);
    assert(del_energy > 0);
    double e = lower_bound;
    double peak_density = 0;
    double peak_e = -1;
    // Seek maximum between low and upper bound.
    while (e < upper_bound){
        double density = y[step].F(e)*e; // electron energy density
        if (density > peak_density){
            peak_density = density;
            peak_e = e;
        }
        e += del_energy; 
    }
    return peak_e;
}
std::vector<double> ElectronRateSolver::approx_regime_peaks(size_t step, double lower_bound, double upper_bound, double del_energy, size_t num_peaks, double min_density){
    assert(del_energy > 0);
    assert(num_peaks > 0);      

    //double min_peak_separation = 400 / Constant::eV_per_Ha; 
    size_t min_sequential = 3;
    std::vector<double> peak_energies;  
    double last_peak_density = INFINITY; // for asserting expected behaviour.
    for(size_t i = 0; i < num_peaks; i++){
        // Seek maximum between low and upper bound.
        double peak_density = -1;
        double peak_e = -1;
        double e = lower_bound;
        while (e < upper_bound){
            // search for nearest peak, by looking for the nearest point that a) occurs after the density has been rising and b) is higher than the following min_sequential points separated by del_e
            e = approx_nearest_peak(step,e,del_energy,min_sequential,lower_bound,upper_bound);
            double density = y[step].F(e)*e; // electron energy density
            // separate peaks by min_peak_separation... doesnt really work since the highest peaks are at the back ah well. TODO
            if (std::find(peak_energies.begin(),peak_energies.end(),e)!= peak_energies.end()){
                //e+= min_peak_separation; //not working atm
                // continue;
                break;
            }
            if (peak_density < density){
                assert(density <= last_peak_density);
                peak_density = density;
                peak_e = e;
            }
        }
        if (peak_density < min_density) 
            peak_e = -1;
        last_peak_density = peak_density;
        peak_energies.push_back(peak_e);
    }
    return peak_energies;
}

double ElectronRateSolver::approx_regime_trough(size_t step, double lower_bound, double upper_bound, double del_energy){
    del_energy = abs(del_energy);
    assert(del_energy > 0);
    double e = lower_bound;
    double trough_density = INFINITY;
    double trough_e = -1;
    // Seek maximum between low and upper bound.
    while (e < upper_bound){
        double density = y[step].F(e)*e; // electron energy density
        if (density < trough_density){
            trough_density = density;
            trough_e = e;
        }
        e += del_energy; 
    }
    if (trough_density < 0) return -1;  // error due to wiggle - minimum is not physical.
    return trough_e;
}



/**
 * @brief Finds the energy bounds of the photoelectron region
 * @details Since it isn't obvious what function approximates the peak at a given time we can't analytically 
 * determine the boundary. Here we determine the inflection point, and choose the region's boundary 
 * to be some point past this.
 */
void ElectronRateSolver::dirac_energy_bounds(size_t step, std::vector<double>& maximums, std::vector<double>& minimums, std::vector<double>& peaks, bool allow_shrinkage,size_t num_peaks, double peak_min_density) {
    #ifdef SWITCH_OFF_DYNAMIC_BOUNDS
    return;
    #endif
    double min_photo_peak_considered = 1500/Constant::eV_per_Ha;  // An energy that is above auger energies but will catch significant peaks. //TODO replace with transition energy of last regimes?
    double peak_search_step_size = 10/Constant::eV_per_Ha;
    // Find peaks
    size_t num_sequential_needed = 3;
    peaks = approx_regime_peaks(step,min_photo_peak_considered,Distribution::get_max_E(),peak_search_step_size,num_peaks,peak_min_density);
    // Set peaks' region boundaries
    //about halfway between the peak and the neighbouring local minima is sufficient for their mins and maxes
    for(size_t i = 0; i < peaks.size(); i++) {
        double e_peak = peaks[i];
        if (e_peak < 0){
            // Found all peaks that suit our restrictions, move outside of grid so that it is ignored. 
            // Note that if allow_shrinkage == false, this removes the region permanently.
            minimums[i] = -99/Constant::eV_per_Ha;
            maximums[i] = -100/Constant::eV_per_Ha;
            continue;
        }

        assert(e_peak >= min_photo_peak_considered);
        // Allow for step size to be smaller if the *PREVIOUS* max/min  is close to the peak.
        // This won't work if the peaks change relative heights (unlikely) or comparing to the initial grid state (insignificant w/ gaussian pulse due to low intensity)
        double up_e_step = max((maximums[i] - e_peak) , 10/Constant::eV_per_Ha)/200;
        up_e_step = min(up_e_step, 10/Constant::eV_per_Ha);
        double down_e_step = min(minimums[i] - e_peak , -10/Constant::eV_per_Ha)/200;            
        down_e_step = max(down_e_step, -10/Constant::eV_per_Ha);

        // Get bounds
        double min_distance = input_params.elec_grid_preset.pulse_omega/10; // the minimum distance from the peak that the region must cover.
        std::cout << "Inflections of peak at " <<e_peak*Constant::eV_per_Ha <<" are... Lwr:";
        double lower_bound = approx_regime_bound(step,e_peak, down_e_step, num_sequential_needed,min_distance,1./4.);
        std::cout <<", Upr: ";
        double upper_bound = approx_regime_bound(step,e_peak, up_e_step, num_sequential_needed,min_distance,1./4.);//peak + 0.912*(peak - lower_bound);  // s.t. if lower_bound is 3/4 peak, upper_bound is 1.1*peak.
        std::cout << std::endl;
        if(allow_shrinkage){
            maximums[i] = -1;
            minimums[i] = INFINITY;
        }           
        if (upper_bound > maximums[i]) maximums[i] = upper_bound;
        if (lower_bound < minimums[i]) minimums[i] = lower_bound;
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
    #ifdef SWITCH_OFF_DYNAMIC_BOUNDS
    return;
    #endif
    // Find e_peak = kT/2
    double min_energy = 0;
        // a) Prevent a change if it's too big. At early times, the MB peak may not yet be taller than the auger peak.
        // b) designed to be below photoeletron peaks, TODO make this input-dependent at least.
    double max_peak_energy = max(10/Constant::eV_per_Ha,min(2*peak,2000/Constant::eV_per_Ha)); // TODO this will start to fail if updating much faster than dynamics does (i.e. low XFEL intensities) - especially since at low intensities low energy peaks become significant! meaning it will double each update. Though it works decently with abdallah neon at 0.25 fs update period.
    // Step size is hundredth of the previous peak past a peak of 1 eV. Note gets stuck on hitch at early times, but not a big deal as the minimum size of new_max lets us get through this. 
    // Alternatively could just implement local max detection for early times.
    double e_step_size = max(1./Constant::eV_per_Ha,peak)/200; 
    double new_peak = approx_regime_peak(step,min_energy,max_peak_energy,e_step_size);
    peak = new_peak;

    double kT = 2*peak;
    // CDF = Γ(3/2)γ(3/2,E/kT)
    double new_min = 0.2922*kT;  // 90% of electrons above this point
    if(_min < new_min || allow_shrinkage){
        _min = new_min;
    }
    //double min_e_seq_range  = 20/Constant::eV_per_Ha;
    //size_t num_sequential_needed = max(min_seq_needed,(int)(min_e_seq_range/e_step_size+0.5)); 
    /* Using an inflection to judge regime end may be causing issues when MB and photo merge.
    size_t num_sequential_needed = 3; // early times are safeguarded from inflections being too small by disallowing shrinkage
    double new_max = approx_regime_bound(step,peak, +e_step_size, num_sequential_needed,5,1./2.);
    */
    double new_max = 2.3208*kT; // 80% of electrons below this point (lower since not as sharp)
    if(_max < new_max || allow_shrinkage)
        _max = std::min(new_max,Distribution::get_max_E());
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
void ElectronRateSolver::transition_energy(size_t step, double& g_min){
    #ifdef SWITCH_OFF_DYNAMIC_BOUNDS
    return;
    #endif
    //TODO assert that this is called after MB and dirac regimes updated.

    // for these functions, if a trough is negative then it defaults to using g_min.
    double new_min = INFINITY;
    for(auto &dirac_peak : regimes.dirac_peaks){
        if (dirac_peak < 0) continue;
        new_min = min(new_min,approx_regime_trough(step,regimes.mb_peak,0.9*dirac_peak,2/Constant::eV_per_Ha)); 
    }
    // no peaks left, we just look for minimum between mb peak and photon energy. 
    if (new_min == INFINITY){
        new_min = min(new_min,approx_regime_trough(step,regimes.mb_peak,input_params.Omega(),2/Constant::eV_per_Ha)); 
    }
    // if(allow_decrease) 
    //     g_min = new_min;
    // else               
    g_min = max(g_min,new_min); // if the previous transition energy was higher, use that (good way to wait out instability).
}