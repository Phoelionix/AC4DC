
#include <vector>
#include <iostream>

using namespace std;

class GridRegions
{
public:
    GridRegions();
    StaticRegion static_low_divergent;
    StaticRegion static_auger;
    StaticRegion static_high_tail; 
    MBRegion dyn_mb;
    DiracRegion dyn_photo;
    void update_regions();

protected:
    void set_static_energies(vector<float> energy_boundaries);
    void grid_update(ofstream& _log);


private:
    int NUM_STATIC_REGIONS = 3;
    int NUM_DYNAMIC_REGIONS = 2;

};

class Region
{
public:
    Region(int num_points, float E_min, float E_max) :
    num_points{num_points},
    E_min{E_min},
    E_max{E_max}
    {
        num_points_left = num_points;
    }  
    float get_point_density(){return (E_max-E_min)/num_points;};
    float get_E_min(){return E_min;}
    float get_E_max(){return E_max;}
    int get_active_points(){return num_points_left;}

private:
    int num_points;
    int num_points_left;
    float E_min;
    float E_max;
};

class StaticRegion : public Region
{
public:
    StaticRegion(int num_points, float E_min, float E_max) : Region(num_points,E_min,E_max)
    {}
};

class MBRegion : public Region
{
public:
    MBRegion(int num_points, float E_min, float E_max) : Region(num_points,E_min,E_max)
    {}

private:
    float get_peak();

};

class DiracRegion : public Region
{
public:
    DiracRegion(int num_points, float E_min, float E_max) : Region(num_points,E_min,E_max)
    {}
private:
    float get_peak();
};