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

#ifdef ELECTRON_TRANSFER_DEBUG
void Space::Clear(){
    F_assigned = false; 
}
#endif

void Space::ElectronTransfer(const size_t& a, const double& rho, single_state_type& sdot,const double& t, const double& cross_r){    
    #ifndef NO_SPATIAL
    if (t < last_time){
        return;
    }

    #ifdef QUANTISED
    // EXTREMELY CRUDE. just using charge of F and not neighbours.
    // Assumes that electrons spread out equally among neighbours (electron sinks)
    assert(F_assigned);
    
    if (cross_r>999999){
        return;
    } 

    double min_e;
    double max_e;

    // double de = 100;
    // size_t num = (int)(Distribution::get_max_E()*Constant::eV_per_Ha+1)%(int)de;
    // de/=Constant::eV_per_Ha;

    // std::vector<std::pair<double,double>> e_t_list; 
    // for (size_t i = 1; i < num+1;i++){
    //     double max_e = de*i;
    //     double v = sqrt(2*max_e);
    //     double t_cross = cross_r/v;

    //     t_cross/=2;

    //     double next_t_cross = anchor_time;
    //     while(next_t_cross <= last_time){
    //       next_t_cross+=t_cross; 
    //     }
    //     if (next_t_cross < t){
    //         e_t_list.push_back(std::pair<double,double>{max_e,t_cross});
    //     }
    // }
    double de = 100;
    size_t num = (int)((Distribution::get_max_E()*Constant::eV_per_Ha)/de+1);

    std::vector<std::pair<double,double>> e_t_list; 
    const double fraction = 0.8; 
    assert(fraction<1);
    for (size_t i = 1; i < num+1;i++){
        double max_e = de*i;
        //double v = sqrt(2*max_e);
        double v = 5.931*sqrt(max_e); // ang/fs, with max_e in eV
        double t_cross = cross_r*Constant::Angs_per_au/v;

        t_cross*=fraction/Constant::fs_per_au;

        double next_t_cross = anchor_time;
        size_t j = 0;
        while(next_t_cross <= last_time){
          next_t_cross+=t_cross;
          j++; 
        }
        if (next_t_cross < t){
            e_t_list.push_back(std::pair<double,double>{max_e/Constant::eV_per_Ha,t_cross});
        }
        if (num == 24327){
            std::cout<<j;
        }
    }
    de/=Constant::eV_per_Ha;
    //const double dt = t-last_time;
    
    for (const auto& e_t: e_t_list){
        //TODO reimplement confinement
        const double max_e = e_t.first;
        const double t_cross = e_t.second; 
        const double min_e = max_e - de;
        //const double time_factor = t_cross/dt; //TODO handle case when time steps haven't been constant size
        const double& time_factor = t_cross;
        //const double& time_factor = t_cross/2;  // TEMPORARY
        for (const auto& source : electron_sources){        
            sdot.F.addSource(a,(*source.first).original_F,source.second,rho,min_e,max_e,time_factor);
            sdot.F.addLoss(a,original_F,source.second,rho,min_e,max_e,time_factor);
        }
    }
    #else
        const double tmp_dt = t-last_time;
        for (const auto& source : electron_sources){        
            sdot.F.addSource(a,(*source.first).original_F,source.second,rho,0,Distribution::get_max_E(),tmp_dt);
            sdot.F.addLoss(a,original_F,source.second,rho,0,Distribution::get_max_E(),tmp_dt);
        }    
    #endif // QUANTISED
    

    #endif //NO_SPATIAL
}

void Space::ElectronTransferV2(const size_t& a, const double& rho, single_state_type& sdot, LossGeometry& l){    
    #ifndef NO_SPATIAL
    // EXTREMELY CRUDE. just using charge of F and not neighbours.
    // Assumes that electrons spread out equally among neighbours (electron sinks)
    assert(F_assigned);
    
    //TODO reimplement confinement
    for (const auto& source : electron_sources){        
        //sdot.F.addSourceV2(a,(*source.first).original_F,l,source.second,rho);
        //sdot.F.addLossV2(a,original_F,l,source.second,rho);
    }

    #endif //NO_SPATIAL
}

void Space::AddBoundary(Space& other_space, CustomLossGeometry boundary_geometry){
    electron_sources.push_back(
        std::pair<Space*, CustomLossGeometry>(&other_space,boundary_geometry));
}

