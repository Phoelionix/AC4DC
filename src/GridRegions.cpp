#include "GridRegions.h"
#include "Constant.h"
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
    float e = 1/Constant::eV_per_Ha;
    regions = {
        Region(1,5*e,10*e,"static"),  // low divergent
        Region(20,10*e,600*e,"static"), // auger
        Region(10,600*e,12000*e,"static"), // high tail
        Region(35,4,10,"dirac"), // Maxwell-boltzmann distribution
        Region(40,4500,6500,"mb") // Photoelectron peak
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
void GridRegions::update_regions(double mb_peak, double mb_width, double dirac_peak, double dirac_width){
     for (size_t r = 0; r < regions.size(); r ++){
        switch(regions[r].get_type()[0]){
            case 's': // static
                return;
            break; 
            case 'd':
                regions[r].update_region(dirac_peak,dirac_width);
            break;
            case 'm':
                regions[r].update_region(mb_peak,mb_width);
            break;
            default:
                std::cout <<"Error, unrecognised region type" <<regions[r].get_type() << std::endl;
                return;
        }
    }
}

void GridRegions::set_static_energies(vector<double> energy_boundaries){

}

void Region::update_region(double new_centre, double new_width){
    if (type == "static") return;    
    E_min = new_centre - new_width/2;
    E_max = new_centre + new_width/2;
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
    return -1;
}