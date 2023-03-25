/**
 * @file GridRegions.h
 * @brief Handles dynamic grid mechanics
 */

#include "Constant.h"
#include "GridSpacing.hpp"
#include <vector>
#include <iostream>

using namespace std;

class Region
{
public:
    // // Dummy region
    // Region(): type{"blank"},power{1}{};
    /**
     * @brief Construct a new Region object
     * 
     * @param num_points 
     * @param E_min in eV!
     * @param E_max in eV!
     * @param type one of: "static", "dirac", "mb"
     */
    Region(int num_points, double E_min, double E_max, const string type) :
    num_points{num_points},
    E_min{E_min/Constant::eV_per_Ha},
    E_max{E_max/Constant::eV_per_Ha},
    type{type},
    power{1}   
    {
        num_points_left = num_points;
    }  
    double get_point_density(){return (E_max-E_min)/num_points;};
    double get_E_min(){return E_min;}
    double get_E_max(){return E_max;}
    const string get_type(){return type;}
    // (just used for sorting)
    bool operator < (Region& other){
        return E_max  < (other.get_E_min());
    }    
    double get_next_knot(double previous_knot);
    void update_region(double new_centre,double new_min, double new_max);

private:
    int num_points;
    int num_points_left;
    double E_min;
    double E_max;
    const string type;
    const int power;
    double get_peak();
    double dirac_peak();
    double mb_peak();
};

class GridRegions
{
public:
    GridRegions();
    void update_regions(FeatureRegimes rgm);

protected:
    void set_static_energies(vector<double> energy_boundaries);
    std::vector<Region> regions;


private:
    int NUM_STATIC_REGIONS = 3;
    int NUM_DYNAMIC_REGIONS = 2;
};

