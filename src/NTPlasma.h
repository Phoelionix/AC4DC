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

class NTPlasma{
    NTPlasma(const int T_POINTS, Grid &lattice);
public:
    double EEcoll(int indx);
private:
    vector<vector<double>> n;
    Grid Lattice;
};

#endif /* end of include guard: NTPLASMA_H */
