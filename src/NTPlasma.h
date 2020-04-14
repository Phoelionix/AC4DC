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


#ifndef NTPLASMA_H
#define NTPLASMA_H

#include <vector>
#include <fstream>
#include "Constant.h"
#include "Grid.h"
#include "Plasma.h" // for dsigmaBEB methods

typedef vector<double> nte_state_t;

class NTPlasma{
    NTPlasma(const int T_POINTS, const Grid &lattice);
public:
    // df(e,t)/dt = Q_eii[f](t) + Q_pi[f](t) + Q_a[f](t) + Q_ee[f] + Q_tbr[f](t) - \Lambda
    // Electron-impact ionisation collision term at timestep n
    nte_state_t Q_eii(int n);
    // Photoionisation flux at timestep n
    nte_state_t Q_pi(int n);
    // Auger effect
    nte_state_t Q_a(int n);
    // Electron-electron scattering
    // nte_state_t Q_ee(int n);
    // Density loss due to electrons leaving the system
    nte_state_t Q_loss(int n);


private:
    vector<nte_state_t> n;
    vector<double> e_grid; // energy values
    Grid Lattice; // Integration grid
};

#endif /* end of include guard: NTPLASMA_H */
