#include "GridRegions.h"
#include "Constant.h"
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
    Region static_low_divergent(1,5*e,10*e,"static");
    Region static_auger(20,10*e,600*e,"static");
    Region static_high_tail(10,600*e,12000*e,"static"); 
    Region dyn_mb(35,4,10,"dirac");
    Region dyn_photo(40,4500,6500,"mb");
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
 */
void GridRegions::grid_update(ofstream& _log){
    update_regions();
    

    std::vector<Region*> stat_regions = {&static_low_divergent,&static_auger,&static_high_tail};  
    std::vector<Region*> dyn_regions = {&dyn_mb,&dyn_photo};  
    Region* finer_region;
    Region* coarser_region;
    // Get each static grid's energies
    for (size_t i = 0; i < stat_regions.size(); i++){
        for(size_t j = 0; j < dyn_regions.size(); j++){
            if(stat_regions[i]->get_point_density() < dyn_regions[j]->get_point_density()){
                // If the dynamic grid is denser, consider points "filled" 
                finer_region = dyn_regions[j];
                coarser_region = stat_regions[i];
            }
            else{
                finer_region = stat_regions[i];
                coarser_region = dyn_regions[j];
            }
        }
        // Coarser region has its number of active points removed
        float prop_filled = (coarser_region->get_E_max()-coarser_region->get_E_min())/(min(finer_region->get_E_max(),coarser_region->get_E_max()) - max(finer_region.get_E_min(),coarser_region.get_E_min()));
        finer_region->set_num_points_left(finer_region->get_num_points());
        coarser_region->set_num_points_left(round(
            coarser_region->get_active_points() - prop_filled*stat_regions[i]->get_num_points()
            ));       
    }

    // Insert dynamic regions into appropriate places. We are assuming that the dynamic grids 
    std::vector<Region*> sorted_regions = stat_regions;
    sorted_regions.insert( sorted_regions.end(), dyn_regions.begin(), dyn_regions.end() );
    std::sort(sorted_regions.begin(),sorted_regions.end()); 

    // Set first index.
    std::vector<double> boundary_i = {1};
    // Set first energy to minimum of all grids present ().
    std::vector<double> boundary_E = {sorted_regions[0]->get_E_min()};
    std::vector<double> powers = {};    
    for(size_t i = 0; i < sorted_regions.size(),i++){
        powers.push(1)
        boundary_i.push(boundary_i[i]+sorted_regions[i]->get_active_points());
        if i != sorted_regions.size()-1{
            boundary_E.push(sorted_regions[i+1]->get_E_min());
        }
        else boundary_E.push(sorted_regions[i]->get_E_max());
    }

    GridBoundaries new_region = {boundary_i,boundary_E,powers}
    this->set_grid_regions(new_region);

    this->compute_cross_sections(_log,True);

}

/**
 * @brief Update individual region energy bounds. Using linear regions.
 * 
 */
void GridRegions::update_regions(){


    E_max = ;
    E_min = ;
}

void GridRegions::set_static_energies(vector<float> energy_boundaries){

}

