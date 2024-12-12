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
    bool confined_system=true; //TODO allow for choice in input file
};

class Space{
public:
    //Space(Distribution* F) : internal_F(F){F_assigned = true;}
    Space(){F_assigned = false;}
    std::vector<std::pair<Space *,CustomLossGeometry>> electron_sources; // second element of each tuple contains information on the ratio of the boundary surface area to the volume of this space.

    #ifdef ELECTRON_TRANSFER_DEBUG
    virtual void Clear();
    #endif
    virtual void ElectronTransfer(size_t a, double rho, single_state_type& sdot); // Modifies F by connected electron sinks original_F. 
    virtual void set_F(const Distribution* F){
        #ifndef NO_SPATIAL
        #ifdef ELECTRON_TRANSFER_DEBUG
        assert(!F_assigned);
        #endif
        internal_F=F;
        original_F = *internal_F;
        original_F_fraction = *internal_F;  // Might need to multiply this by some independent variable if doing some weird geometries but no need for now.  Good approx. for concentric shells.
        //original_F_fraction*=(1./electron_sources.size());

        F_assigned=true;
        #endif // NO_SPATIAL
    }
    void AddBoundary(Space& other_space, CustomLossGeometry boundary_geometry);
    // Distribution* F(){
    //     assert(F_assigned);
    //     return internal_F;
    // }


private:
    Distribution original_F;
    Distribution original_F_fraction;
    const Distribution* internal_F; // TODO delete
    bool F_assigned;
    


};


class Void_Space : public Space{
    public:
    void ElectronTransfer(size_t a, double rho, single_state_type& sdot) override{}
    void set_F(const Distribution* F) override{};
    #ifdef ELECTRON_TRANSFER_DEBUG
    void Clear() override{};
    #endif
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