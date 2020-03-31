#include "NTPlasma.h"

NTPlasma::NTPlasma(const int T_POINTS, Grid &lattice): Lattice(lattice) {
    n.resize(T_POINTS);
    for (size_t i = 0; i < n.size(); i++) {
        n[i].resize(lattice.size());
    }
}

// Electron-Electron collision integral
double NTPlasma::EEcoll(int indx){
    return 0.;
}
