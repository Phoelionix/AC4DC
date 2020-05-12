#ifndef ADAMS_CXX_H
#define ADAMS_CXX_H

#include <vector>
#include "Adams_arrays.h"
#include <iostream>

using namespace std;
//
// template<typename T>
// using sysfunc_t = void (*)(const T&, T&, const double);

template<typename T>
class IVPSolver{
public:
    IVPSolver();
    void setup(const T& initial_state, const double dt);
    double dt;
    void print();

    // void export_to(string fname);
protected:
    vector<T> y;
    vector<double> t;
    virtual void sys(const T& q, T& qdot, const double t) =0;
    T zero_y;

};

template<typename T>
class Adams_BM : public IVPSolver<T>{
public:
    Adams_BM(int order=5);
    double iterate(double t_initial, size_t npoints); // returns final time

private:
    int order;
    const double* b_AB;
    const double* b_AM;

    void step(int n); // Predictor-corrector multistep method
    void step_rk4(int n); // Runge-Kutta 4th order calculator
};

template<typename T>
IVPSolver<T>::IVPSolver(){
    y.resize(1);
    t.resize(1);
}

template<typename T>
void IVPSolver<T>::setup(const T& initial_state, double _dt ){
    this->y[0] = initial_state;
    // Makes a zero vector in a mildly spooky way
    this->zero_y = initial_state; // do this to make the underlying std::vector large enough
    this->zero_y = 0.; // set it to Z E R O
    if (_dt < 1E-16){
        cerr<<"WARN: step size "<<dt<<"is smaller than machine precision"<<endl;
    }
    this->dt = _dt;
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
Adams_BM<T>::Adams_BM(int _order):
    IVPSolver<T>()
    {
    if (_order < 2 || _order > AdamsArrays::MAX_ADAMS_ORDER){
        cerr<<"ERROR: Adams order may not be greater than "<<AdamsArrays::MAX_ADAMS_ORDER;
    }
    this->order = _order;
    this->b_AB = AdamsArrays::AB_COEFF[order];
    this->b_AM = AdamsArrays::AM_COEFF[order];
}


// Computes y_n+1 based only on y_n
// For initialisng the multistep method
template<typename T>
void Adams_BM<T>::step_rk4(int n){
    T k1, k2, k3, k4, tmp;
    this->sys(this->y[n], k1, this->t[n]);
    this->sys(this->y[n]+k1*0.5, k2, this->t[n]+this->dt*0.5);
    this->sys(this->y[n]+k2*0.5, k3, this->t[n]+this->dt*0.5);
    this->sys(this->y[n]+k3,     k4, this->t[n]+this->dt    );

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
void Adams_BM<T>::step(int n){
    // Predicts the value y_n+1
    // Adams-Bashforth predictor routine:

    T tmp=this->zero_y;
    T ydot;

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
    this->sys(tmp, tmp, this->t[n+1]);
    tmp *= b_AM[0];
    // Now tmp goes back to bein an aggregator
    for (size_t i = 1; i < order; i++) {
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

// Fixed-step method
// Maybe get fancier if it's a real problem
template<typename T>
double Adams_BM<T>::iterate(double t_initial, size_t npoints){
    double h = this->dt;

    if (h < 1E-16){
        cerr<<"WARN: step size "<<h<<"is smaller than machine precision"<<endl;
    } else if (h < 0) {
        cerr<<"ERROR: step size is negative!"<<endl;
        return 0;
    }

    // Set up the containters
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
    return this->t[npoints-1];
}



#endif /* end of include guard: ADAMS_CXX_H */
