#include "NTPlasma.h"

NTPlasma::NTPlasma(const int T_POINTS, const Grid &lattice): Lattice(lattice) {
    this->n.resize(T_POINTS);
    int N = this->Lattice.size();
    for (size_t i = 0; i < n.size(); i++) {
        this->n[i].resize(N);
    }
}

// df(e,t)/dt = Q_eii[f](t) + Q_pi[f](t) + Q_a[f](t) + Q_ee[f] + Q_tbr[f](t) - \Lambda

// Electron-impact ionisation collision term at timestep n
nte_state_t NTPlasma::Q_eii(int n){

}

nte_state_t NTPlasma::Q_pi(int n){

}
// Photoionisation flux at timestep n
nte_state_t NTPlasma::Q_a(int n){

}
// Electron-electron scattering
// nte_state_t Q_ee(int n);
// Density loss due to electrons leaving the system
nte_state_t NTPlasma::Q_loss(int n){

}
