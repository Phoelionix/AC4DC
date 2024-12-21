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
    const static char planar = 2;
    const static char experimental = 3;
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
    virtual void ElectronTransfer(const size_t& a, const double& rho, single_state_type& sdot, const double& t,const double& cross_r); // Modifies F by connected electron sinks original_F. 
    virtual void ElectronTransferV2(const size_t& a, const double& rho, single_state_type& sdot, LossGeometry& l);
    virtual void set_F(const Distribution* F){
        #ifndef NO_SPATIAL
        #ifdef ELECTRON_TRANSFER_DEBUG
        assert(!F_assigned);
        #endif
        internal_F=F;
        original_F = *internal_F;
        
        F_assigned=true;
        #endif // NO_SPATIAL
    }
    void AddBoundary(Space& other_space, CustomLossGeometry boundary_geometry);
    // Distribution* F(){
    //     assert(F_assigned);
    //     return internal_F;
    // }

    virtual void set_anchor_time(const double& time){anchor_time = time;set_last_t(anchor_time);}
    virtual void set_last_t(const double& time){last_time = time;}


private:
    Distribution original_F;
    //Distribution original_F_fraction;
    const Distribution* internal_F;
    bool F_assigned;
    double anchor_time;
    double last_time;


};


class Void_Space : public Space{
    public:
    void ElectronTransfer(const size_t& a, const double& rho, single_state_type& sdot, const double& t,const double& cross_r) override{}
    void ElectronTransferV2(const size_t& a, const double& rho, single_state_type& sdot, LossGeometry& l) override{}
    void set_F(const Distribution* F) override{};
    void set_last_t(const double& time)override {};
    void set_anchor_time(const double& time) override{};
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
        case SpatialArrangement::planar:
            os << "Planes";
            break;
        case SpatialArrangement::experimental:
            os << "Experimental";
            break;
        default:
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
        case 'p':
            sa.mode = SpatialArrangement::planar;
            break;
        case 'e':
            sa.mode = SpatialArrangement::experimental;
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