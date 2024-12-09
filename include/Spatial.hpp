/** @file Spatial.hpp
 * @authors Spencer Passmore
 * @brief 
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

#ifndef AC4DC_CXX_SPACE_H
#define AC4DC_CXX_SPACE_H

#include <vector>
#include <assert.h>
#include <stdexcept>
#include "FreeDistribution.h"

class Space{
public:
    //Space(Distribution& F) :F(F){}
    //Distribution& F;
    Space(Distribution* F) : internal_F(F){F_assigned = true;}
    Space(){F_assigned = false;}
    std::vector<Space *> electron_sinks;


    void CalculateOutgoingElectrons(); // Decrease deltaF and increase deltaF of electron_sinks by same amount.
    void ApplyChanges();
    void set_F(Distribution* F){
        assert(F_assigned=false);
        internal_F=F;
        F_assigned=true;
    }
    Distribution* F(){
        assert(F_assigned);
        return internal_F;
    }
private:
    Distribution deltaF;
    Distribution* internal_F;
    bool F_assigned;


};


// class SpatialConnection{
// public: 
//     Space& A, deltaA; 
//     Space& B, deltaB;
    

//     void ElectronTransfer();

//     void CalculateChanges();
//     void ApplyChanges();


// };



#endif