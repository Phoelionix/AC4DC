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
};

namespace ode {
template<typename T>
class Hybrid : public Adams_BM<T>{
    public:
    Hybrid<T>(unsigned int order=4, double rtol=1e-5, unsigned max_iter = 2000): 
        Adams_BM<T>(order), stiff_rtol(rtol), stiff_max_iter(max_iter){};
    const static unsigned MAX_BEULER_ITER = 50;  // unused?
    bool good_state = true;
    bool euler_exceeded = false;  // Whether exceeded the max number of euler iterations (on prev. step)
    FeatureRegimes regimes;
    double timestep_reached = 0;       
    private:
    virtual void sys_ee(const T& q, T& qdot, double t) =0;
    // virtual void Jacobian2(const T& q, T& qdot, double t) =0; 
    protected:

    double stiff_rtol = 1e-4;
    double intolerable_stiff_err =0;//0.5;
    unsigned stiff_max_iter = 200;

    Checkpoint old_checkpoint;
    Checkpoint checkpoint;
    // For gaussian pulses, this stores the times that we decreased dt so that it can be increased later.
    //TODO this should be saved if want to load from a point in time, but really it's meant to replace loading.
    std::vector<double> times_to_increase_dt;  
    
    void run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update,const int steps_per_grid_transform);  // TODO clean up bootstrapping. -S.P.
    void iterate(ofstream& _log, double t_initial, double t_final, const double t_resume, const int steps_per_time_update,const int steps_per_grid_transform);
    /// Unused
    void backward_Euler(unsigned n); 
    void step_stiff_part(unsigned n);
    // More virtual funcs defined by ElectronRateSolver:
    virtual void set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step = 0){std::cout << "Error, attempted to use unset virtual function set_up_grid_and_compute_cross_sections." <<std::endl;} 
    virtual size_t load_checkpoint_and_decrease_dt(ofstream &log, size_t current_n, Checkpoint _checkpoint){std::cout << "Error, attempted to use unset virtual function load_checkpoint_and_decrease_dt" <<std::endl;}
    virtual void increase_dt(ofstream &log, size_t current_n){std::cout << "Error, attempted to use with unset virtual function increase_dt" <<std::endl;}
    virtual state_type get_ground_state(){}
};

template<typename T>
// t_resume = the time to resume simulation from if loading a sim. -S.P.
// _log only used for cross-section recalcs atm.
void Hybrid<T>::iterate(ofstream& _log, double t_initial, double t_final, const double t_resume, const int steps_per_time_update, const int steps_per_grid_transform) {

    if (this->dt < 1E-16) {
        std::cerr<<"WARN: step size "<<this->dt<<"is smaller than machine precision"<<std::endl;
    } else if (this->dt < 0) {
        throw std::runtime_error("Step size is negative!");
    }

    size_t npoints = (t_final - t_initial)/this->dt + 1;
    
    bool resume_sim = (t_resume == t_initial) ? false : true;

    checkpoint = {this->order, Distribution::get_knot_energies(), this->regimes};
    if (resume_sim){
        // The time step size does not depend on previous run's time step size. i.e. step size is same as if there was no loading.
        // TODO implement assertion that density isn't empty.
        size_t resume_idx = (t_resume-t_initial)/this->dt ;
        checkpoint = {resume_idx, Distribution::get_knot_energies(),this->regimes}; // THIS CURRENTLY DOESNT WORK because we don't store the knot energies.
        npoints -= (resume_idx + 1);
        npoints += this->t.size(); // 
    }
    old_checkpoint = checkpoint; 

    npoints = (npoints >= this->order) ? npoints : this->order;

    // Set up the containers
    this->t.resize(npoints,INFINITY);
    this->y.resize(npoints);

    // Set up the t grid       
    this->t[0] = t_initial;

    for (size_t n=1; n<npoints; n++){
        if (resume_sim && this->t[n] <= t_resume){
            continue; // Don't reset already simulated states
        }
        this->t[n] = this->t[n-1] + this->dt;
    }
    this->run_steps(_log,t_resume, steps_per_time_update,steps_per_grid_transform);
}


// Overrides the underlying Adams method, adding a more refined but computationally expensive treatment for the stiff Q^{ee} contribution to deltaf.
/**
 * @brief  Iterates through the (t_final-t_initial = simulation_timespan)/dt timesteps.
 * .s
 * @tparam T 
 */
template<typename T>
void Hybrid<T>::run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update,const int steps_per_grid_transform){
    assert(this->y.size() == this->t.size());
    assert(this->t.size() >= this->order);

    // activate the display
    //Display::reactivate();
    
    Plotting py_plotter;

    if (t_resume < this->t[this->order]){
        // initialise enough points for multistepping to get going
        for (size_t n = 0; n < this->order; n++) {
            this->step_rk4(n);
        }
    }
    // Run those steps 
    std::stringstream tol;
    tol << "[ sim ] Implicit solver uses relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n\r";
    std::cout << tol.str();  // Display in regular terminal even after ncurses screen is gone.
    Display::header += tol.str(); 
    std::stringstream display_stream(Display::header, ios_base::app | ios_base::out); // text shown that updates with frequency 1/steps_per_time_update.
    std::stringstream popup_stream(std::string(), ios_base::app | ios_base::out);  // text displayed during step of special events like grid updates
    Display::create_screen(); 
    size_t n = this->order - 1; // n can go back if load checkpoint.
    size_t increasing_dt = 0;
    int checkpoint_gap = 100; 
    while (n < this->t.size()-1) {
        n++;
        if (this->t[n+1] <= t_resume) continue; // Start with n = last step.

        if ((n-this->order)%steps_per_time_update == 0){
            // Display info 
            display_stream.str(Display::header); // clear display string
            display_stream<< "\n\r"
            << "--- Press BACKSPACE/DEL to end simulation and save the data ---\n\r"   
            << "[ sim ] Current timestep size ="<<this->dt*Constant::fs_per_au<<" fs\n\r"   
            << "[ sim ] t="
            << std::left<<std::setfill(' ')<<std::setw(6)
            << this->t[n] * Constant::fs_per_au << " fs\n\r" 
            << "[ sim ] " <<Distribution::size << " knots currently active\n\r"; 
            // << flush; 
            Display::show(display_stream);
        }        
        
        // live plotting
        if ((n-this->order+1)%25 == 0){ // this rate won't be significant except for samples with no heavy elements
            size_t num_pts = 4000;
            py_plotter.plot_frame(Distribution::get_energies_eV(num_pts),this->y[n].F.get_densities(num_pts,Distribution::get_knot_energies()));
        }        
        // store checkpoint
        if ((n-this->order+1)%checkpoint_gap == 0){  // 
            old_checkpoint = checkpoint;
            checkpoint = {n, Distribution::get_knot_energies(),this->regimes};
        }

        if (euler_exceeded){            
            if (!increasing_dt){
                popup_stream <<"\nMax euler iterations exceeded, reloading checkpoint at "<<this->t[old_checkpoint.n]*Constant::fs_per_au<< "and decreasing dt.\n\r";
            }
            else{
                popup_stream << "\nReloading checkpoint at "<<this->t[old_checkpoint.n]*Constant::fs_per_au
                <<", as increasing step size led to euler iterations being exceeded. Will attempt again after "<<2*checkpoint_gap<< "steps.\n\r";
                times_to_increase_dt.pop_back();  // remove the extra time added.
                times_to_increase_dt.back() = this->t[n] + 2*checkpoint_gap*this->dt;  
                increasing_dt = 0;                
            }
            Display::show(display_stream,popup_stream);
            // Reload at checkpoint's n, updating the n step, and decreasing time step length
            n = this->load_checkpoint_and_decrease_dt(_log,n,old_checkpoint);  // virtual function overridden by ElectronRateSolver
            good_state = true;
            euler_exceeded = false;
            checkpoint = old_checkpoint;
        }
        else if (!times_to_increase_dt.empty() && times_to_increase_dt.back() < this->t[n]){
            if (increasing_dt > 0){
                increasing_dt--;
                //successfully increased step size
                times_to_increase_dt.pop_back();
            }

            // // Attempt to increase step size (i.e. decrease num remaining steps)
            this->increase_dt(_log,n);
            // Flag that we are testing to see if still converging.
            size_t num_steps_to_check = 10;
            increasing_dt = num_steps_to_check;
        }
        
        this->step_nonstiff_part(n); 
        
        // this->y[n+1].from_backwards_Euler(this->dt, this->y[n], stiff_rtol, stiff_max_iter);
        this->step_stiff_part(n);

        if (euler_exceeded && increasing_dt){
            n = this->load_checkpoint_and_decrease_dt(_log,n,old_checkpoint);
            good_state = true;
            euler_exceeded = false;
        }

        // Dynamic grid updater // TODO pop in function.
        if ((n-this->order+1)%steps_per_grid_transform == 0){
            popup_stream << "\nUpdating grid... \n\r"; 
            Display::show(display_stream,popup_stream);
            std::cout.setstate(std::ios_base::failbit);  // disable character output
            // The latest step is n + 1, so we decide our new grid based on that step, then transform N = "order" of the prior points to the new basis.
            this->set_up_grid_and_compute_cross_sections(_log,false,n+1); // virtual function overridden by ElectronRateSolver
            std::cout.clear();
            // Transform enough previous points needed to get going to new basis
            std::vector<double> new_energies = Distribution::get_knot_energies();
            // zero_y is used as the starting state for each step and represents the ground state, so we need to reset it so it has the right knots.
            this->zero_y = this->get_ground_state();
            cout << endl;            
            for (size_t m = n+2 - this->order; m < n+2; m++) {
                // reload back to the old energies so we can use transform_basis(). yes it's goofy :/
                Distribution::load_knots_from_history(n);
                this->y[m].F.transform_basis(new_energies);
            }  
             // The next containers are made to have the correct size, as the initial state is set to tmp=zero_y and sdot is set to an empty state.       
        }    

        auto ch = wgetch(Display::win);
        if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){   
            flushinp(); // flush buffered inputs
            popup_stream <<"\n\rExiting early... press backspace/del again to confirm or any other key to cancel \n\r";
            Display::show(display_stream,popup_stream, false);
            nodelay(Display::win, false);
            ch = wgetch(Display::win);  // note implicitly refreshes screen
            if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){
                this->y.resize(n);
                this->t.resize(n);    
                break;
            }
            nodelay(Display::win, true);
            werase(Display::win);
        }       
        if (popup_stream.rdbuf()->in_avail() != 0){
            // clear stream
            popup_stream.str(std::string());
        }    

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
    T old = this->y[n];
    old *= -1.;
    old += this->y[n+1];

    if(!good_state) return;

    // Adams-Moulton step - (Same as one called for end of good part of system).
    // Here, tmp is the y_n+1 guessed in the preceding step, handle i=0 explicitly:
    T tmp;
    // Now tmp goes back to being an aggregator
    tmp *= 0;
    for (int i = 1; i < this->order; i++) {
        T ydot;
        this->sys_ee(this->y[1+n-i], ydot, this->t[1+n-i]);
        ydot *= this->b_AM[i];
        tmp += ydot;
    }
    
    //this->y[n+1] = this->y[n] + tmp * this->dt;
    tmp *= this->dt;
    tmp += this->y[n];

    // tmp now stores the additive constant for y_n+1 = tmp + h b0 f(y_n+1, t_n+1)

    // Picard iteration
    T prev;
    double diff = stiff_rtol*2;
    unsigned idx=0;
    while (diff > stiff_rtol && idx < stiff_max_iter){
        // solve by Picard iteration
        prev = this->y[n+1];
        prev *= -1;
        T dydt;
        this->sys_ee(this->y[n+1], dydt, this->t[n+1]);
        dydt *= this->b_AM[0]*this->dt;
        this->y[n+1] = tmp;
        this->y[n+1] += dydt;
        prev += this->y[n+1];
        diff = prev.norm()/this->y[n+1].norm();
        idx++;
    }
    if(idx==stiff_max_iter){
        std::cerr<<"Max Euler iterations exceeded, err = "<<diff<<std::endl;
        if (diff > intolerable_stiff_err){
            //std::cerr << "Max error ("<<intolerable_stiff_err<<") exceeded, ending simulation early." <<std::endl; // moved to 
            this->good_state = false;
            this->timestep_reached = this->t[n+1]*Constant::fs_per_au; // t[n+1] is equiv. to t in bound !good_state case, where error condition this is modelled off is found.
            this->euler_exceeded = true;  
        }
    }
    this->y[n+1] += old;
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
    sys_ee(this->y[n+1], this->y[n+1], this->t[n+1]);
    this->y[n+1] *= this->dt;
    this->y[n+1] += this->y[n]; 

    // Picard iteration TODO seek better way to do this - reduce modularisation - reduce to solving linear system?
    double diff = stiff_rtol*2;
    while (diff > stiff_rtol && idx < stiff_max_iter){ // TODO more efficient mapping to use by considering known form of equation  c.f. ELENDIF
        T tmp = this->y[n+1];
        sys_ee(tmp, this->y[n+1], this->t[n+1]);
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