#include "SystemSolver.h"

SystemSolver::SystemSolver(double final_time, double _dt){
    num_steps = (int) final_time/dt;
    dt = _dt;
    T.resize(num_steps);
    Y.resize(num_steps);
    // Set the initial T to this.
    T[0] = 0;
    abm.initialize( sys, T[0], 0, dt );
}

void SystemSolver::solve(){
    for (size_t i = 0; i < nsteps-1; i++) {
                        //in    t   out    dt
        abm.do_step( sys, Y[i], T[i], Y[i+1],  dt );
        T[i+1] = T[i] + dt;
    }
}

// The Big One: Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)
void SystemSolver::sys(const state_type& s, state_type& sdot, const double t){
    for (size_t i = 0; i < s.P.size(); i++) {
        sdot[i] = 0;
        for (size_t j = 0; j < s.P.size(); j++) {
            sdot[i] += W[i][j] * s.P[j];
        }
    }
}
