#include "DynamicRegions.h"
#include <iostream>
#include <vector>
#include <variant>
#include <math.h>
#include <execution>


/**
 * @brief 
 * @details  
 * Static regions (named for the feature they are modelling when dynamic region not present):
 * - Low E divergent region
 * - Auger region
 * - Tail region
 * Dynamic regions:
 * - MB region
 * - Photoelectron region
 */

GridRegions::GridRegions(){
    // Initialise regions
    // Carbon (v. divergent example) good grid:
    // indices  1|30|55|65|85|105|145|150 
    // energies 4|10|50|200|500|4500|6500|10000   
    regions = {
        Region(1,5,10,"static"),  // low divergent
        Region(25,10,50, "static"), // low support (MB is very fine early on, grid doesn't let us make it suddenly become super fine.)
        Region(10,50,200,"static"), 
        Region(35,200,600,"static"), // auger
        Region(40,600,6500,"static"),  //TODO change 6500 to be the xray energy.
        Region(10,6500,15000,"static"), // high tail
        Region(35,-1,-1,"mb"), // Maxwell-boltzmann distribution
        Region(60,-1,-1,"dirac") // Photoelectron peak
    };

} 
/**
 * @brief Updates the regions, then transforms the grid based on the regions' updated properties.
 * @details 
 *  Static regions (named for the feature they are modelling when dynamic region not present):
 * - Low E divergent region
 * - Auger region
 * - Tail region
 * Dynamic regions:
 * - MB region
 * - Photoelectron region
 * @param _log  
 */




/**
 * @brief Update dynamic regions' energy bounds. (Using linear regions).
 * @details currently places centre of region on peak.
 */
void GridRegions::update_regions(FeatureRegimes rgm){
    std::cout << "[Dynamic Grid] Updating regions " << std::endl;
    for (size_t r = 0; r < regions.size(); r ++){
        switch(regions[r].get_type()[0]){
            case 's': // static
            break; 
            case 'd':
                regions[r].update_region(rgm.dirac_peak,rgm.dirac_min,rgm.dirac_max);
            break;
            case 'm':
                regions[r].update_region(rgm.mb_peak,rgm.mb_min,rgm.mb_max);
            break;
            default:
                std::cout <<"Error, unrecognised region type" <<regions[r].get_type() << std::endl;
        }
    }
}

void GridRegions::set_static_energies(vector<double> energy_boundaries){

}

void Region::update_region(double new_centre, double new_min, double new_max){
    assert(type != "static");    
    if (E_min != new_min || E_max != new_max){
        std::cout <<"energy range of region of type '"<<type<<"' updated from:\n" 
        <<E_min*Constant::eV_per_Ha<<" - "<<E_max*Constant::eV_per_Ha << std::endl;    
        E_min = new_min;
        E_max = new_max;
        std::cout << "to:\n" 
        <<E_min*Constant::eV_per_Ha<<" - "<<E_max*Constant::eV_per_Ha << std::endl;
    }
    else std::cout<<"Energy range of region of type '"<<type<<"' had NO update." <<std::endl;
}

// powers not implemented yet since doesn't seem necessary
/**
 * @brief Returns next knot, or -1 if it's outside the region's (energy) range. 
 * 
 * @param previous_knot 
 * @return double 
 */
double Region::get_next_knot(double previous_knot){
    if(power != 1){std::cout << "WARNING powers not implemented for dynamic grid" <<std::endl;}
    double next_knot = previous_knot + get_point_density();
    if ((E_min <= next_knot && next_knot <= E_max)){
        return next_knot;
    }
    if(E_min > previous_knot) return E_min;
    return -1;
}