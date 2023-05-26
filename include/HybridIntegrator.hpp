/** @file HybridIntegrator.hpp
 * @authors Alaric Sanders & Spencer Passmore
 * @brief Defines the Hybrid class which adds a Moulton step after each step from the inherited method of Adams_BM.
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

#ifndef HYBRID_INTEGRATE_HPP
#define HYBRID_INTEGRATE_HPP

#include "AdamsIntegrator.hpp"
#include "Constant.h"
#include "FreeDistribution.h"
#include "RateSystem.h"
#include "Display.h"
#include "Plotting.h"
#include "GridSpacing.hpp" // for checkpoint
#include <Eigen/Dense>
#include <iostream>

struct Checkpoint {
    size_t n;
    std::vector<double> knots;
    FeatureRegimes regimes;
    std::vector<state_type> last_few_states;
};

namespace ode {
template<typename T>
class Hybrid : public Adams_BM<T>{
    public:
    Hybrid<T>(unsigned int order=4, double rtol=1e-5, unsigned max_iter = 50): 
        Adams_BM<T>(order), stiff_rtol(rtol), stiff_max_iter(max_iter){};
    const static unsigned MAX_BEULER_ITER = 50;  // unused?
    FeatureRegimes regimes;
    double timestep_reached = 0;       
    private:
    virtual void sys_ee(const T& q, T& qdot) =0;
    // virtual void Jacobian2(const T& q, T& qdot, double t) =0; 
    protected:

    double stiff_rtol = 1e-4;
    unsigned stiff_max_iter = 300;//200; 
    double intolerable_stiff_divergence =0;//0.5; // allowed divergence if stiff_max_iter exceeded
    
    // stiff ode intermediate steps (i.e. steps it does without the nonstiff part)
    int msi;
    #ifndef NO_EE
    int num_stiff_ministeps = 500;
    #else
    int num_stiff_ministeps = 1;
    #endif

    // Grid/timestep dynamics    
    size_t steps_per_grid_transform;    
    Checkpoint old_checkpoint; // the checkpoint that is loaded.
    Checkpoint checkpoint; 
    int checkpoint_gap = 100; 
    size_t increasing_dt = 0;
    // TODO need to separate good_state from euler_exceeded, currently we have the case of good_state is false and 
    // euler is exceeded, which means euler exceeded. And we have the case of good state is false and euler 
    // exceeded is false, which means a nan was encountered.    
    bool good_state = true;
    bool euler_exceeded = false;  // Whether exceeded the max number of euler iterations (on prev. step)    
    // [not in use as bugged] For gaussian pulses, this stores the times that we decreased dt so that it can be increased later.
    //TODO this should be saved if want to load from a point in time, but really it's meant to replace loading.
    std::vector<double> times_to_increase_dt;  
    
    void run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update);  // TODO clean up bootstrapping. -S.P.
    void iterate(ofstream& _log, double t_initial, double t_final, const double t_resume, const int steps_per_time_update);
    void initialise_transient_y(int n); // Approximates initial stiff intermediate steps
    std::vector<T> transient_y; // stores transient/intermediate steps for stiff solver
    // More virtual funcs defined by ElectronRateSolver:
    virtual state_type get_ground_state()=0;
    virtual void pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update)=0;
    virtual int post_ode_step(ofstream& _log, size_t& n)=0;
    /// Unused
    void backward_Euler(unsigned n); 
    void step_stiff_part(unsigned n);
};

template<typename T>
// t_resume = the time to resume simulation from if loading a sim. -S.P.
// _log only used for cross-section recalcs atm.
void Hybrid<T>::iterate(ofstream& _log, double t_initial, double t_final, const double t_resume, const int steps_per_time_update) {
    if (this->dt < 1E-16) {
        std::cerr<<"WARN: step size "<<this->dt<<"is smaller than machine precision"<<std::endl;
    } else if (this->dt < 0) {
        throw std::runtime_error("Step size is negative!");
    }

    size_t npoints = (t_final - t_initial)/this->dt + 1;
    
    bool resume_sim = (t_resume == t_initial) ? false : true;

    size_t resume_idx = 0;
    if (resume_sim){
        for (size_t n=1; n<npoints; n++){
            if (this->t[n] >= t_resume){
                resume_idx = n;
                break;
            }
        }
        // The time step size does not depend on previous run's time step size. i.e. step size is same as if there was no loading.
        // TODO implement assertion that density isn't empty.
        size_t resume_idx_if_const_dt = (t_resume-t_initial)/this->dt ;
        npoints -= (resume_idx_if_const_dt + 1);
        npoints += this->t.size(); // 
        // Set checkpoint to be at the starting step
        std::vector<state_type>::const_iterator start_vect_idx = this->y.begin() - this->order + resume_idx;  
        std::vector<state_type>::const_iterator end_vect_idx = this->y.begin() + resume_idx;  
        
        std::vector<state_type> check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);
        checkpoint = {resume_idx, Distribution::get_knot_energies(),this->regimes,check_states};
    }
    old_checkpoint = checkpoint; 

    npoints = (npoints >= this->order) ? npoints : this->order;

    // Set up the containers
    this->t.resize(npoints,INFINITY);
    this->y.resize(npoints);

    // Set up the t grid       
    this->t[0] = t_initial;

    for (size_t n=1; n<npoints; n++){
        if (resume_sim && n <= resume_idx){
            continue; // Don't reset already simulated states
        }
        this->t[n] = this->t[n-1] + this->dt;
    }
    this->run_steps(_log,t_resume, steps_per_time_update);
}


// Overrides the underlying Adams method, adding a more refined but computationally expensive treatment for the stiff Q^{ee} contribution to deltaf.
/**
 * @brief  Iterates through the (t_final-t_initial = simulation_timespan)/dt timesteps.
 * .s
 * @tparam T 
 */
template<typename T>
void Hybrid<T>::run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update){
    assert(this->y.size() == this->t.size());
    assert(this->t.size() >= this->order);

    size_t n = 0; // Current step index.
    if (t_resume < this->t[this->order]){
        // 4th Order Runge-Kutta 
        // initialise enough points for multistepping to get going
        while(n < this->order) {
            this->step_rk4(n);
            n++;
        }
        // Almost guaranteed we did not load a simulation, so set first checkpoint. 
        //std::vector<state_type> check_states = std::vector<state_type>(this->y.begin(), this->y.begin()+this->order);
        std::vector<state_type>::const_iterator start_vect_idx = this->y.begin();  
        std::vector<state_type>::const_iterator end_vect_idx = this->y.begin() + this->order;  
        std::vector<state_type> check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);

        checkpoint = {this->order, Distribution::get_knot_energies(), this->regimes, check_states};
        old_checkpoint = checkpoint;
    }
    else{
        n+=this->order;
    }

    // Start with n = last step with data.
    while(this->t[n] < t_resume)
        n++;
    initialise_transient_y((int)n);

    // Set up display
    std::stringstream tol;
    tol << "[ sim ] Implicit solver uses relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n\r";
    std::cout << tol.str();  // Display in regular terminal even after ncurses screen is gone.
    Display::header += tol.str(); 
    Display::display_stream = std::stringstream(Display::header, ios_base::app | ios_base::out); // text shown that updates with frequency 1/steps_per_time_update.
    Display::popup_stream = std::stringstream(std::string(), ios_base::app | ios_base::out);  // text displayed during step of special events like grid updates
    Display::create_screen(); 

    // Run hybrid multistepping (nonstiff -> bound dynamics, stiff -> free dynamics). 
    while (n < this->t.size()-1) {
        // 1. Display general information
        // 2. Handle live plotting
        // 3. Handle dynamic updates to dt (and checkpoints, which also handle case of NaN being encountered by solver).
        this->pre_ode_step(_log, n,steps_per_time_update);
        this->step_nonstiff_part(n); 
        #ifdef DEBUG_BOUND
        for(size_t a = 0; a < this->y[n+1].atomP.size();a++)
            for(size_t i=0;i < this->y[n+1].atomP[a].size();i++){
                assert(this->y[n+1].atomP[a][i] >= 0);
            }        
        #endif    
        // this->y[n+1].from_backwards_Euler(this->dt, this->y[n], stiff_rtol, stiff_max_iter);
        this->step_stiff_part(n);

        // 1. Handle dynamic grid updates 
        // 2. Check if user wants to exit simulation early  
        if (this->post_ode_step(_log,n) == 1){
            break; // user asked to end the simulation early.
        } 
        #ifdef DEBUG_BOUND
        for(size_t a = 0; a < this->y[n+1].atomP.size();a++)
            for(size_t i=0;i < this->y[n+1].atomP[a].size();i++){
                assert(this->y[n+1].atomP[a][i] >= 0);
            }        
        #endif   
        n++;     
    }
    Display::close();
    std::cout<<"\n";
    // std::cout <<"[ sim ] Implicit solver used relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n";
    std::cout<<"[ sim ] final t = "<<this->t.back() * Constant::fs_per_au<<" fs"<< endl;  
}
/**
 * @brief Use a true implicit method to estimate change to bad part of system based solely on its own action.
 * 
 * @param n n+1 is the index of the step to predict.
 * @return template<typename T> 
 */
template<typename T>
void Hybrid<T>::step_stiff_part(unsigned n){
    if(!good_state) return;
    
    // Loop explanation:
    // msi == 'ministep' index (that was *last* updated)
    // transient_y contains the ministeps needed for next step's calculation.
    // transient_y = [1,2,3,4] for cubic splines  (order elements)
    // We use the order=3 most recent elements, and replace the oldest element with the newest
    double mini_dt = this->dt/(double)num_stiff_ministeps;
    // Hence we are targeting (msi+1)%(order) in Adams-Moulton step
    T old;
    int old_msi = msi;
    // 'rel idx' = relative index within transient y
    int last_rel_idx = msi%(this->order);
    int next_rel_idx; 
    // Store bound contribution to y dot
    old = this->y[n];
    old *= -1.;
    old += this->y[n+1];  

    #ifdef DEBUG
    assert(this->y[n].F[3] == transient_y[last_rel_idx].F[3]);
    #endif

    while (msi < old_msi + num_stiff_ministeps){  // msi = n if num_stiff_ministeps = 1. maybe call mini_n?
        
        // Adams-Moulton step I believe -S.P.
        T tmp;
        tmp *= 0;
        // tmp acts as an aggregator
        for (int i = 1; i < this->order; i++){  // work through last N=order-1 ministeps
            T ydot; // ydot stores the change this loop.
            this->sys_ee(transient_y[(1-i+msi)%(this->order)], ydot); 
            ydot *= this->b_AM[i];
            tmp += ydot;
        }
        //this->y[n+1] = this->y[n] + tmp * dt;
        tmp *= mini_dt;
        tmp += transient_y[last_rel_idx];

        // tmp now stores the additive constant:  y_n+1 = tmp + h b0 f(y_n+1, t_n+1)

        // Picard iteration
        next_rel_idx = (msi+1)%(this->order); 
        T prev;
        double diff = stiff_rtol*2;
        unsigned idx=0;
        while (diff > stiff_rtol && idx < stiff_max_iter){
            // solve by Picard iteration
            prev = transient_y[next_rel_idx];
            prev *= -1;
            T dydt;
            this->sys_ee(transient_y[next_rel_idx], dydt);
            dydt *= this->b_AM[0]*(mini_dt);
            transient_y[next_rel_idx] = tmp;
            transient_y[next_rel_idx] += dydt;
            prev += transient_y[next_rel_idx];
            diff = prev.norm()/transient_y[next_rel_idx].norm(); // Seeking convergence
            idx++;
        }
        if(idx==stiff_max_iter){
            //std::cerr<<"Max Euler iterations exceeded, err = "<<diff<<std::endl;
            if (diff > intolerable_stiff_divergence){
                //std::cerr << "Max error ("<<intolerable_stiff_divergence<<") exceeded, ending simulation early." <<std::endl; // moved to 
                this->good_state = false;
                this->timestep_reached = this->t[n+1]*Constant::fs_per_au; // t[n+1] is equiv. to t in bound !good_state case, where error condition this is modelled off is found.
                this->euler_exceeded = true;  
            }
        }
        msi++;
        last_rel_idx = msi%(this->order);
    }
    transient_y[next_rel_idx] += old;
    this->y[n+1] = transient_y[next_rel_idx];
    
}

/// Initialises intermediate steps needed for stiff solver to get going.
template<typename T>
void Hybrid<T>::initialise_transient_y(int n) {
    transient_y.resize(this->order);

    msi = this->order-1; // last index of simulation's first "loop" of transient y. 
    if (this->order/num_stiff_ministeps <= 1){
        // Approximate intermediate steps with change between last time steps.
        // We store y[n] and flag it as the most recent ministep updated. 
        T ydot;
        ydot = this->y[n-1];
        ydot*= -1.;
        ydot += this->y[n];
        ydot *= 1./(double)num_stiff_ministeps;    
        for(int i = 0; i < this->order; i++){  // Technically we don't need to fill the first idx but we may as well.
            // y[n-i] = y[n] - ydot*i
            transient_y[msi-i] = this->y[n];
            T tmp;
            tmp = ydot; tmp*=-i;
            transient_y[msi-i] += tmp;
        }
    }
    // else if <can get a whole number of steps>{

    //}
    else{
        // Not enough ministeps to remain within a single step.
        for(int i = 0; i < this->order; i++){
            transient_y[msi-i] = this->y[n-i];
        }
    }
    assert(this->y[n].F[3] == transient_y[msi].F[3]);
}

template<typename T>
void Hybrid<T>::backward_Euler(unsigned n){
    // Assumes that y_n+1 contains a guess based on sys, and estimates a solution to y_n
    unsigned idx = 0;
    //old = y[n+1] - y[n] 
    T old = this->y[n];
    old *= -1.;
    old += this->y[n+1];  // 
    //old caches the step from the regular part

    // Guess. (Naive Euler)
    sys_ee(this->y[n+1], this->y[n+1]);
    this->y[n+1] *= this->dt;
    this->y[n+1] += this->y[n]; 

    // Picard iteration TODO seek better way to do this - reduce modularisation - reduce to solving linear system?
    double diff = stiff_rtol*2;
    while (diff > stiff_rtol && idx < stiff_max_iter){ // TODO more efficient mapping to use by considering known form of equation  c.f. ELENDIF
        T tmp = this->y[n+1];
        sys_ee(tmp, this->y[n+1]);
        this->y[n+1] *= this->dt;
        this->y[n+1] += this->y[n];
        
        tmp *= -1;
        tmp += this->y[n+1];
        diff = tmp.norm()/this->y[n].norm();
        idx++;
    }
    this->y[n+1] += old;
    if(idx==stiff_max_iter) std::cerr<<"Max Euler iterations exceeded (B)"<<std::endl;
}


    





} // end of namespace 'ode'


#endif