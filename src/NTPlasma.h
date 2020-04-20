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
#include "Numerics.h"
#include "Grid.h"
#include "Plasma.h" // for dsigmaBEB methods

typedef vector<double> nte_state_t;




class NTPlasma{
public:
    NTPlasma(const int T_POINTS, const Grid &lattice);
    // Integrated versions of the elementwise methods below go here
    // (perhaps in a federated public method)


private:
    double fintegral(double (&g)(double));
    // df(e,t)/dt = Q_eii[f](e) + Q_pi[f](e, t) + Q_a[f](e) + Q_ee[f](e) + Q_tbr[f](e) - Q_loss(e)
    // Electron-impact ionisation collision term at timestep n
    double Q_eii(nte_state_t &f, double e);
    // Photoionisation flux at timestep n
    double Q_pi(nte_state_t &f, double e);
    // Auger effect
    double Q_a(nte_state_t &f, double e);
    // Electron-electron scattering
    // double Q_ee(int n);
    // Density loss due to electrons leaving the system
    double Q_loss(nte_state_t &f, double e);

    vector<nte_state_t> e_dist;
    Grid Lattice;
};

#endif /* end of include guard: NTPLASMA_H */
