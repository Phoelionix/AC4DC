/**
 * @file ElectronSolverGridDynamics.cpp
 * @author Spencer Passmore
 * @brief 
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
    double last_density = y[step].F(start_energy)*start_energy;
    size_t num_sequential = 0;
    if(min_sequential < 1) min_sequential = 1;
    // Seek local minimum.
    while (num_sequential < min_sequential+1){
        e += del_energy;
        if (e < min){
            local_min = min; break;}
        if (e > max){
            local_min = max; break;}

        double density = y[step].F(e)*e;  // energy density. Plot breakage seems to depend on this. TODO check why
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
        double density = y[step].F(e)*e; // energy density
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
    std::cout << inflection;
    // At min_distance, go double as 
    double A = 1/min_inflection_fract; // region between peak and inflection take up (at least if using sqrt thing) min_inflection_frac after min distance.
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
        double density = y[step].F(e)*e; // energy density
        if (density > peak_density){
            peak_density = density;
            peak_e = e;
        }
        e += del_energy; 
    }
    return peak_e;
}

double ElectronRateSolver::approx_regime_trough(size_t step, double lower_bound, double upper_bound, double del_energy){
    del_energy = abs(del_energy);
    assert(del_energy > 0);
    double e = lower_bound;
    double trough_density = INFINITY;
    double trough_e = -1;
    // Seek maximum between low and upper bound.
    while (e < upper_bound){
        double density = y[step].F(e)*e; // energy density
        if (density < trough_density){
            trough_density = density;
            trough_e = e;
        }
        e += del_energy; 
    }
    return trough_e;
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
            double density = y[step].F(r.energy)*r.energy; // energy density
            if (density >= peak_density){
                peak = r.energy;
                peak_density = density; 
            }
            // Get bounds
            std::cout << "Inflections of peak at" <<r.energy*Constant::eV_per_Ha <<" from "<<atom.name<<" are:";
            double lower_bound = approx_regime_bound(step,r.energy, -e_step_size, num_sequential_needed,1000,1./4.);
            std::cout <<", ";
            double upper_bound = approx_regime_bound(step,r.energy, +e_step_size, num_sequential_needed,600,1./5.);//peak + 0.912*(peak - lower_bound);  // s.t. if lower_bound is 3/4 peak, upper_bound is 1.1*peak.
            std::cout << std::endl;
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
    double max_energy = 2000/Constant::eV_per_Ha;
    // Step size is hundredth of the previous peak past a peak of 1 eV. Note gets stuck on hitch at early times, but not a big deal as the minimum size of new_max lets us get through this. 
    // Alternatively could just implement local max detection for early times.
    double e_step_size = max(1.,peak)/100/Constant::eV_per_Ha; 
    peak = approx_regime_peak(step,min_energy,max_energy,e_step_size);
    double kT = 2*peak;
    // CDF = Γ(3/2)γ(3/2,E/kT)
    double new_min = 0.2922*kT;  // 90% of electrons above this point
    if(_min < new_min || allow_shrinkage){
        _min = new_min;
    }
    double min_e_seq_range  = 20/Constant::eV_per_Ha;
    int min_seq_needed = 3; // -> at later times, it seeks min_seq_needed*peak/10 energy to confirm inflection.
    size_t num_sequential_needed = max(min_seq_needed,(int)(min_e_seq_range/e_step_size+0.5)); 
    std::cout << "upper inflection of MB is:";
    double new_max = approx_regime_bound(step,peak, +e_step_size, num_sequential_needed,5,1./2.); 
    std::cout << std::endl;
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
void ElectronRateSolver::transition_energy(size_t step, double& g_min){
    //TODO assert that this is called after MB and dirac regimes updated.
    double new_min = approx_regime_trough(step,regimes.mb_peak,regimes.mb_max,2); 
    // if(allow_decrease) 
    //     g_min = new_min;
    // else               
    g_min = max(g_min,new_min); // if the previous transition energy was higher, use that (good way to wait out instability).
}