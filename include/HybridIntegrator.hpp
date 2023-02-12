/** @file HybridIntegrator.hpp
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
#include <Eigen/Dense>

namespace ode {
template<typename T>
class Hybrid : public Adams_BM<T>{
    public:
    Hybrid<T>(unsigned int order=4, double rtol=1e-5, unsigned max_iter = 2000): 
        Adams_BM<T>(order), stiff_rtol(rtol), stiff_max_iter(max_iter){};
    const static unsigned MAX_BEULER_ITER = 50;
    bool good_state = true;
    double timestep_reached = 0;       
    private:
    virtual void sys_ee(const T& q, T& qdot, double t) =0;
    // virtual void Jacobian2(const T& q, T& qdot, double t) =0; 
    protected:

    double stiff_rtol = 1e-4;
    double intolerable_stiff_err =stiff_rtol*10;
    unsigned stiff_max_iter = 200;
     


    void run_steps();
    void iterate(double t_initial, double t_final);
    /// Unused
    void backward_Euler(unsigned n); 
    void step_stiff_part(unsigned n);
};

// template <typename T>
// void Hybrid<T>::step_backward_euler(int n){
//     T dndt;
//     // this->sys_bound(this->y[n], dndt, this->t[n]);
//     dndt *=this->dt;
//     this->y[n+1] = this->y[n];
//     this->y[n+1] += dndt;
// }

template<typename T>
void Hybrid<T>::iterate(double t_initial, double t_final) {

    if (this->dt < 1E-16) {
        std::cerr<<"WARN: step size "<<this->dt<<"is smaller than machine precision"<<std::endl;
    } else if (this->dt < 0) {
        throw std::runtime_error("Step size is negative!");
    }

    size_t npoints = (t_final - t_initial)/this->dt + 1;
    
    npoints = (npoints >= this->order) ? npoints : this->order;
    // Set up the containters
    this->t.resize(npoints);
    this->y.resize(npoints);

    // Set up the t grid
    this->t[0] = t_initial;

    for (size_t n=1; n<npoints; n++){
        this->t[n] = this->t[n-1] + this->dt;
    }
    this->run_steps();
}

// Overrides the underlying Adams method, adding a more refined but computationally expensive treatment for the stiff Q^{ee} contribution to deltaf.
/**
 * @brief  Iterates through the (t_final-t_initial = simulation_timespan)/dt timesteps.
 * 
 * @tparam T 
 */
template<typename T>
void Hybrid<T>::run_steps(){
    assert(this->y.size() == this->t.size());
    assert(this->t.size() >= this->order);

    // initialise enough points for multistepping to get going
    for (size_t n = 0; n < this->order; n++) {
        this->step_rk4(n);
    }
    // Run those steps
    std::cout << "[ sim ] Implicit solver uses relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<std::endl;
    std::cout << "[ sim ]                       ";
    for (size_t n = this->order; n < this->t.size()-1; n++) {
        std::cout << "\r[ sim ] t="
                  << std::left<<std::setfill(' ')<<std::setw(6)
                  << this->t[n] * Constant::fs_per_au << std::flush;  //TODO check if this multiplication is taxing.
     
        this->step_nonstiff_part(n); 
        
        // this->y[n+1].from_backwards_Euler(this->dt, this->y[n], stiff_rtol, stiff_max_iter);
        this->step_stiff_part(n);
    }
    std::cout<<std::endl;
}

/// Use a true implicit method to estimate change to bad part of system based solely on its own action.
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
            std::cerr << "Max error ("<<intolerable_stiff_err<<") exceeded, ending simulation early." <<std::endl;
            this->good_state = false;
            this->timestep_reached = this->t[n+1]*Constant::fs_per_au; // t[n+1] is equiv. to t in bound !good_state case
        }
    }
    this->y[n+1] += old;
}

template<typename T>
void Hybrid<T>::backward_Euler(unsigned n){
    // Assumes that y_n+1 contains a guess based on sys, and estimates a solution to y_n
    unsigned idx = 0;
    T old = this->y[n];
    old *= -1.;
    old += this->y[n+1];
    //old caches the step from the regular part

    // Guess. (Naive Euler)
    sys_ee(this->y[n+1], this->y[n+1], this->t[n+1]);
    this->y[n+1] *= this->dt;
    this->y[n+1] += this->y[n]; 

    double diff = stiff_rtol*2;
    while (diff > stiff_rtol && idx < stiff_max_iter){
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