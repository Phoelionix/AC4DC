#ifndef NTPLASMA_H
#define NTPLASMA_H

#include "Grid.h"

class NTPlasma{
    NTPlasma(const int T_POINTS, Grid &lattice);
public:

private:
    vector<vector<double>> n;
    Grid Lattice;
}

#endif /* end of include guard: NTPLASMA_H */
