/** @file Spatial.cpp
 * @authors Spencer Passmore
 * @brief @copybrief Spatial.hpp
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

#include "Spatial.hpp"

void Space::Clear(){
    deltaF = 0;
    F_assigned = false; 
}
void Space::CalculateOutgoingElectrons(size_t a, const LossGeometry &l, double rho){    
    // EXTREMELY CRUDE. just using charge of F and not neighbours.
    // Assumes that electrons spread out equally among neighbours (electron sinks)
    assert(F_assigned);

    Distribution F_transfer; // representing the electrons leaving the volume
    F_transfer = 0;

    F_transfer.addLoss(a,*internal_F,l,rho);
    deltaF+=F_transfer;
    
    // divvy it up equally.
    F_transfer*=(-1);
    F_transfer*=(1/electron_sinks.size());
    for (Space* source : electron_sinks){
        (*source).deltaF+=F_transfer;
        deltaF+= F_transfer;
    }
}

void Space::ApplyChanges()
{
    *internal_F+=deltaF;     
}

