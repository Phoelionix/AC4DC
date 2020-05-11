#ifndef ADAMS_CXX_H
#define ADAMS_CXX_H

#include <vector>
#include "Adams_arrays.h"
#include <iostream>
#include <cstdlib>


using namespace std;

template<typename T>
class IVPSolver{
public:
    IVPSolver(T& initial_state, T (* f)(const T& q, double t), double dt);

    double dt;
    void print();
    // void export_to(string fname);
protected:
    vector<T> y;
    vector<double> t;
    T (* func)(const T& q, double t);
    T zero_y;

};

template<typename T>
class Adams_BM : public IVPSolver<T>{
public:
    Adams_BM(T& initial_state, T (* f)(const T& q, double t), double dt=1E-4, int order=5);
    long solve(double t_final, double t_initial=0); // Returns array size
    int get_order();

private:
    int order;
    const double* b_AB;
    const double* b_AM;

    void step(int n); // Predictor-corrector multistep method
    void step_rk4(int n); // Runge-Kutta 4th order calculator
};

template<typename T>
IVPSolver<T>::IVPSolver(T& initial_state, T (* f)(const T& q, double t), double dt){
    this->func = f;
    y.resize(1);
    t.resize(1);
    this->y[0] = initial_state;
    this->dt = dt;

    // Makes a zero vector in a mildly spooky way
    this->zero_y = initial_state; // do this to make the underlying std::vector large enough
    this->zero_y = 0.; // set it to Z E R O
}

template<typename T>
void IVPSolver<T>::print(){
    cout<<"t\ty\n"<<endl;
    for (size_t i = 0; i < y.size(); i++) {
        cout<<t[i]<<'\t'<<y[i]<<endl;
    }
}

///////////////////////////////////

template<typename T>
int Adams_BM<T>::get_order(){
    return this->order;
}

template<typename T>
Adams_BM<T>::Adams(T& initial_state, T (* f)(const T& q, double t), double dt, int order):
    IVPSolver<T>(initial_state, f, dt){
    if (order < 2 || order > AdamsArrays::MAX_ADAMS_ORDER){
        cerr<<"ERROR: Adams order may not be greater than "<<AdamsArrays::MAX_ADAMS_ORDER;
    }
    if (dt < 1E-16){
        cerr<<"WARN: step size "<<dt<<"is smaller than machine precision"<<endl;
    }
    this->order = order;
    this->b_AB = AdamsArrays::AB_COEFF[order];
    this->b_AM = AdamsArrays::AM_COEFF[order];
}

// Computes y_n+1 based only on y_n
// For initialisng the multistep method
template<typename T>
void Adams_BM<T>::step_rk4(int n){
    T k1 = this->func(this->y[n],this->t[n]) * this->dt;
    T k2 = this->func(this->y[n]+k1/2, this->t[n]+this->dt/2) * this->dt;
    T k3 = this->func(this->y[n]+k2/2, this->t[n]+this->dt/2) * this->dt;
    T k4 = this->func(this->y[n]+k3, this->t[n]+this->dt) * this->dt;
    this->y[n+1] = this->y[n] + k1/6 + k2/3 + k3/3 + k4/6;
}

template<typename T>
void Adams_BM<T>::step(int n){
    // Predicts the value y_n+1
    // Adams-Bashforth predictor routine:

    T tmp=this->zero_y;

    for (size_t i = 0; i < order; i++) {
        tmp += this->func(this->y[n-i], this->t[n-i]) * b_AB[i];
    }

    // Weird syntax here is done in order to avoid calls to
    // operator+ and orparator*, which may unnecessarily create large
    // local variables.

    // y_n+1 = y_n + dt*sum_{i=0}^s-1 b_i f_(n-i)
    tmp *= this->dt;
    tmp += this->y[n];

    // Adams-Moulton corrector step
    // Here, tmp is the y_n+1 guessed in the preceding step, handle i=0 explicitly:
    tmp = this->func(tmp, this->t[n+1]) * b_AM[0];
    // Now tmp goes back to bein an aggregator
    for (size_t i = 1; i < order; i++) {
        tmp += this->func(this->y[n-i+1], this->t[n-i+1]) * b_AM[i];
    }
    // Store the corrected value
    //this->y[n+1] = this->y[n] + tmp * this->dt;
    tmp *= this->dt;
    tmp += this->y[n];
    this->y[n+1] = tmp;
}

// Fixed-step method
// Maybe get fancier if it's a real problem
template<typename T>
long Adams_BM<T>::solve(double t_initial, long npoints){
    double h = this->dt;

    if (h < 1E-16){
        cerr<<"WARN: step size "<<h<<"is smaller than machine precision"<<endl;
    } else if (h < 0) {
        cerr<<"ERROR: step size is negative!"<<endl;
        return 0;
    }

    // Set up the containters
    double t_final = t_initial + npoints*this->dt;
    this->t.resize(npoints);
    this->y.resize(npoints);

    // Set up the t grid
    for (size_t j = 0; j < npoints; j++) {
        this->t[j] = t_initial + h*j;
    }

    // initialise enough points for multistepping to get going
    for (size_t n = 0; n < order; n++) {
        step_rk4(n);
    }
    // Run those steps
    for (size_t n = order - 1; n < npoints; n++) {
        step(n);
    }
    return npoints;
}



#endif /* end of include guard: ADAMS_CXX_H */
