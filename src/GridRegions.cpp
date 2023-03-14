


#include "GridRegions.h"
#include "Constant.h"

// void ElectronRateSolver::high_energy_stability_check(){
//     // if(bad){good_state = false}
// }


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
    StaticRegion static_low_divergent(1,5*e,10*e);
    StaticRegion static_auger(20,10*e,600*e);
    StaticRegion static_high_tail(10,600*e,12000*e); 
    MBRegion dyn_mb(35,4,10);
    DiracRegion dyn_photo(40,4500,6500);
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

    auto stat_regions = {&static_low_divergent,&static_auger,&static_high_tail};  
    auto dyn_regions = {&MBrRegion,&DiracRegion};
    // Get each static grid's energies
    for(size_t i = 0; i < stat_regions.size(),i++){
        for(size_t j = 0; i < dyn_regions.size(),i++){
        if stat_regions[i].get_point_density() < dyn_regions[j].get_point_density()
            // If the dynamic grid is denser, consider points "filled" 
            finer_region = dyn_regions[j]
            coarser_region = stat_regions[i]
        }
        else{
            finer_region = stat_regions[i]
            coarser_region = dyn_regions[j]
        }
        float prop_filled = (coarser_region.get_E_max()-coarser_region.get_E_min())/(min(finer_region.get_E_max(),coarser_region.get_E_max()) - max(finer_region.get_E_min(),coarser_region.get_E_min()));
        finer_region.num_points_left = finer_region.num_points;
        coarser_region.num_points_left = round((1-prop_filled)*static_regions[i].num_points);       
    }

    // Insert dynamic regions into appropriate places, we can ignore 0 point regions
    auto * unsorted_region = new auto[state_regions.size()+dyn_regions.size()];
    std::copy(state_regions, state_regions + state_regions.size(), result);
    std::copy(dyn_regions, dyn_regions + state_regions.size(), unsorted_region + state_regions.size());
    auto sorted_regions = {};
    size_t processed = 0
    while unsorted_region.size() > 0{
        for(size_t j = 0; j < unsorted_region.size(),j++){
                if unsorted_region.num_points_left != 0{
                    if(sorted_regions.find(unsorted_region[j]) != sorted_regions.end()){
                        continue
                    }
                }
                else{
                    // skip empty regions
                    unsorted_region.remove(unsorted_region[j])
                    continue
                }
            bool smallest_E = True;
            for(size_t i = 0; i < unsorted_region.size(),i++){
                if unsorted_region[j].get_max_E() > unsorted_region[i].get_max_E(){
                    smallest_E = False;
                    break;
                }
            }
            if smallest_E{
                sorted_regions.push(unsorted_region[j]);
                unsorted_region.remove(unsorted_region[j])
            }
        }
    }
    // with the sorted grids, get the grid region params
    std::vector<double> boundary_i = {1};
    std::vector<double> boundary_E = {sorted_regions[0].get_min_E()};
    std::vector<double> powers = {};    
    for(size_t i = 0; i < sorted_regions.size(),i++){
        powers.push(1)
        boundary_i.push(boundary_i[i]+sorted_regions[i].get_active_points());
        if i != sorted_regions.size()-1{
            boundary_E.push(sorted_regions[i+1].get_min_E());
        }
        else boundary_E.push(sorted_regions[i].get_max_E());
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

