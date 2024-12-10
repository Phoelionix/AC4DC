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
    F_assigned = false; 
}
void Space::ElectronTransfer(size_t a, const LossGeometry &l, double rho, single_state_type& sdot){    

    // EXTREMELY CRUDE. just using charge of F and not neighbours.
    // Assumes that electrons spread out equally among neighbours (electron sinks)
    assert(F_assigned);
    sdot.F.addLoss(a,original_F,l,rho);
    
    for (Space* source : electron_sources){        
        sdot.F.addSource(a,(*source).original_F_fraction,l,rho);
    }
}

