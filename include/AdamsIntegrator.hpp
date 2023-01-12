/**
 * @file AdamsIntegrator.hpp
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

#ifndef ADAMS_CXX_H
#define ADAMS_CXX_H

#include <vector>
#include "Adams_arrays.h"
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <stdexcept>

#define ODE_MAX_MEOMRY_USE 4000000000

//
// template<typename T>
// using sysfunc_t = void (*)(const T&, T&, const double);
namespace ode{

template<typename T>
class IVPSolver{
public:
    IVPSolver();
    void setup(const T& initial_state, double _dt, double _tolerance = 1e-2);
    double dt;
    double step_tolerance;

    // void export_to(string fname);
protected:
    std::vector<T> y;
    std::vector<double> t;
    virtual void sys(const T& q, T& qdot, double t) =0;
    T zero_y;

};

template<typename T>
IVPSolver<T>::IVPSolver() {
    // TODO: GLobal refactor of y and t to have more obvious names
    y.resize(1);
    t.resize(1);
}

template<typename T>
void IVPSolver<T>::setup(const T& initial_state, double _dt, double _step_tolerance ) {
    step_tolerance = _step_tolerance;
    this->y[0] = initial_state;
    // Makes a zero std::vector in a mildly spooky way
    this->zero_y = initial_state; // do this to make the underlying structure large enough
    this->zero_y *= 0.; // set it to Z E R O
    if (_dt < 1E-16) {
        std::cerr<<"WARN: step size "<<dt<<"is smaller than machine precision"<<std::endl;
    }
    this->dt = _dt;
}


///////////////////////////////////

template<typename T>
class Adams_BM : public IVPSolver<T>{
public:
    Adams_BM(unsigned int order=4);
    void iterate(double t_initial, double t_final);

protected:
    void run_steps();
    void step(int n); // Predictor-corrector multistep method
    void step_rk4(int n); // Runge-Kutta 4th order calculator
    void step_euler(int n); // Basic explicit Euler stepper
    unsigned int order;
    const double* b_AB;
    const double* b_AM;
};

template<typename T>
Adams_BM<T>::Adams_BM(unsigned int _order):
    IVPSolver<T>()
    {
    if (_order < 2 || _order > AdamsArrays::MAX_ADAMS_ORDER) {
        std::cerr<<"ERROR: Adams order may not be greater than "<<AdamsArrays::MAX_ADAMS_ORDER;
    }
    this->order = _order;
    this->b_AB = AdamsArrays::AB_COEFF[order];
    this->b_AM = AdamsArrays::AM_COEFF[order];
}

template <typename T>
void Adams_BM<T>::step_euler(int n){
    T dndt;
    this->sys(this->y[n], dndt, this->t[n]);
    dndt *=this->dt;
    this->y[n+1] = this->y[n];
    this->y[n+1] += dndt;
}


// Computes y_n+1 based only on y_n
// For initialisng the multistep method
template<typename T>
void Adams_BM<T>::step_rk4(int n) {
    T k1, k2, k3, k4, tmp;
    this->sys(this->y[n], k1, this->t[n]);
    //tmp = this->y[n]+k1*0.5
    k1 *= this->dt;
    tmp = k1;
    tmp *= 0.5;
    tmp += this->y[n];
    this->sys(tmp, k2, this->t[n]+this->dt*0.5);
    // tmp = this->y[n]+k2*0.5
    k2 *= this->dt;
    tmp = k2;
    tmp *=0.5;
    tmp += this->y[n];
    this->sys(tmp, k3, this->t[n]+this->dt*0.5);
    // tmp = this->y[n] + k3;
    k3 *= this->dt;
    tmp = k3;
    tmp += this->y[n];
    this->sys(tmp, k4, this->t[n]+this->dt    );
    
    k4 *= this->dt;
    //y[n+1] = y[n] + (k1*(1./6) + k2*(1./3) + k3*(1./3) + k4*(1./6)) * this->dt;;
    
    
    k1 *= 1./6;
    k2 *= 1./3;
    k3 *= 1./3;
    k4 *= 1./6;
    this->y[n+1] = k1;
    this->y[n+1] += k2;
    this->y[n+1] += k3;
    this->y[n+1] += k4;
    
    this->y[n+1] *= this->dt;
    this->y[n+1] += this->y[n];
}

template<typename T>
void Adams_BM<T>::step(int n) {
    // Predicts the value y_n+1
    // Adams-Bashforth predictor routine:

    T tmp, ydot;
    tmp = this->zero_y;

    for (size_t i = 0; i < order; i++) {
        this->sys(this->y[n-i], ydot, this->t[n-i]);
        ydot *= b_AB[i];
        tmp += ydot;
    }

    // Weird syntax here is done in order to avoid calls to
    // operator+ and orparator*, which may unnecessarily create large
    // local variables.

    // y_n+1 = y_n + dt*sum_{i=0}^s-1 b_i f_(n-i)
    tmp *= this->dt;
    tmp += this->y[n];

    // Adams-Moulton corrector step
    // Here, tmp is the y_n+1 guessed in the preceding step, handle i=0 explicitly:
    T tmp2(tmp);
    this->sys(tmp2, tmp, this->t[n+1]);
    tmp *= b_AM[0];
    // Now tmp goes back to being an aggregator
    for (int i = 1; i < order; i++) {
        this->sys(this->y[n-i+1], ydot, this->t[n-i+1]);
        ydot *= b_AM[i];
        tmp += ydot;
    }
    // Store the corrected value
    //this->y[n+1] = this->y[n] + tmp * this->dt;
    tmp *= this->dt;
    tmp += this->y[n];
    this->y[n+1] = tmp;
}

template<typename T>
void Adams_BM<T>::iterate(double t_initial, double t_final) {

    if (this->dt < 1E-16) {
        std::cerr<<"WARN: step size "<<this->dt<<"is smaller than machine precision"<<std::endl;
    } else if (this->dt < 0) {
        throw std::runtime_error("Step size is negative!");
    }

    size_t npoints = (t_final - t_initial)/this->dt + 1;
    
    npoints = (npoints >= order) ? npoints : order;
    // Set up the containters
    this->t.resize(npoints);
    this->y.resize(npoints);

    // Set up the t grid
    this->t[0] = t_initial;

    for (int n=1; n<npoints; n++){
        this->t[n] = this->t[n-1] + this->dt;
    }
    this->run_steps();
}

template<typename T>
void Adams_BM<T>::run_steps(){
    assert(this->y.size() == this->t.size());
    assert(this->t.size() >= order);

    // initialise enough points for multistepping to get going
    for (int n = 0; n < order; n++) {
        //this->step_rk4(n);
	this->y[n+1] = this->y[0];
    }
    // Run those steps
    std::cout << "[ sim ]                       ";
    for (int n = this->order; n < this->t.size()-1; n++) {
        std::cout << "\r[ sim ] t="
                  << std::left<<std::setfill(' ')<<std::setw(6)
                  << this->t[n] << std::flush;
        this->step(n);
    }
    std::cout<<std::endl;

}


} // End of namespace: ode



#endif /* end of include guard: ADAMS_CXX_H */
