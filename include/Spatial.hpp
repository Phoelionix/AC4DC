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
#include "LossGeometry.hpp"
#include "RateSystem.h"

// 
struct SpatialArrangement{
    const static char concentric_shells = 0;
    const static char cube = 1;
    const static char unknown = 101;
    int mode = 101;
};

class Space{
public:
    //Space(Distribution* F) : internal_F(F){F_assigned = true;}
    Space(){F_assigned = false;}
    std::vector<Space *> electron_sinks;


    void Clear();
    virtual void CalculateOutgoingElectrons(size_t a, const LossGeometry &l, double rho); // Decrease deltaF and increase deltaF of electron_sinks by same amount.
    virtual void ApplyChanges();
    virtual void set_F(Distribution* F){
        assert(F_assigned==false);
        internal_F=F;
        F_assigned=true;

        deltaF = Distribution()=0;
    }
    // Distribution* F(){
    //     assert(F_assigned);
    //     return internal_F;
    // }


private:
    Distribution deltaF;
    Distribution* internal_F;
    bool F_assigned;


};


class Void_Space : public Space{
    public:
    void CalculateOutgoingElectrons (size_t a, const LossGeometry &l, double rho) override{}
    void ApplyChanges() override{};
    void set_F(Distribution* F) override{};
};

namespace{
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, const SpatialArrangement& sa) {
        switch (sa.mode)
        {
        case SpatialArrangement::concentric_shells:
            os << "Concentric shells";
            break;
        case SpatialArrangement::cube:
            os << "Cubes";
            break;
            os << "Unknown geometry";
            break;
        }
        return os;
    }

    [[maybe_unused]] std::istream& operator>>(std::istream& is, SpatialArrangement& sa) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No spatial arrangement specifier provided, defaulting to concentric shells"<<std::endl;
            sa.mode = SpatialArrangement::concentric_shells;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 's':
            sa.mode = SpatialArrangement::concentric_shells;
            break;
        case 'c':
            sa.mode = SpatialArrangement::cube;
            break;
        default:
            throw std::runtime_error("Unrecognised spatial arrangement type");
            break;
        }
        return is;
    }
}



// class SpatialConnection{
// public: 
//     Space& A, deltaA; 
//     Space& B, deltaB;
    

//     void ElectronTransfer();

//     void CalculateChanges();
//     void ApplyChanges();


// };



#endif