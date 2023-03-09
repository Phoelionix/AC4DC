/**
 * @file ElectronRateSolver.cpp
 * @author Alaric Sanders 
 * @brief @copybrief ElectronRateSolver.h
 * @details @copydetail ElectronRateSolver.h
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
// (C) Alaric Sanders 2020

#include "ElectronRateSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "SplineIntegral.h"
#include <fstream>
#include <algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <math.h>
#include <omp.h>
#include "config.h"

// IO functions
void ElectronRateSolver::save(const std::string& _dir) {
    string dir = _dir; // make a copy of the const value
    dir = (dir.back() == '/') ? dir : dir + "/";

    saveFree(dir+"freeDist.csv");
    saveFreeRaw(dir+"freeDistRaw.csv");
    saveBound(dir);
    saveBoundRaw(dir);

    std::vector<double> fake_t; // TODO double check why I called this fake_t, probably doesn't make sense now. -S.P.
    int num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size(); // Fineness of output is only limited by num time steps.
    float t_fineness = (simulation_start_time - simulation_end_time) / num_t_points;
    float previous_t = t[0];
    int i = -1;
    while (i <  static_cast<int>(t.size())-1){  //TODO make this some constructed function or something -S.P. 
        i++;
        if(t[i] < previous_t + t_fineness){
            continue;
        }        
        fake_t.push_back(t[i]);
        previous_t = t[i];
    }
    pf.save(fake_t,dir+"intensity.csv");

    // fake_t.resize(0);
    // num_t_points = num_t_points = t.size();
    // t_idx_step = t.size() / num_t_points;
    // for (size_t i=0; i<num_t_points; i++) {
    //     fake_t.push_back(t[i*t_idx_step]);
    // }
    // pf.save(fake_t,dir+"intensityRaw.csv");    

}

void ElectronRateSolver::log_extra_details(ofstream & _log){
    // Saves details pertaining to the simulation's execution to file fname
    if(_log.is_open()){
        cout << "[ Details ] Logging run-specific details..."<<endl;
        _log << endl << "[ Rate Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << endl;
        _log.flush();
    }
}

void ElectronRateSolver::saveFree(const std::string& fname) {
    // Saves a table of free-electron dynamics to file fname
    ofstream f;
    cout << "[ Free ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Density @ energy (eV):" <<endl;
    f << "#           | "<<Distribution::output_energies_eV(this->input_params.Out_F_size())<<endl;

    assert(y.size() == t.size());
    int num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size(); // Fineness of output is only limited by num time steps.
    float t_fineness = (simulation_start_time - simulation_end_time)  / num_t_points;
    float previous_t = t[0];
    int i = -1;
    while (i <  static_cast<int>(t.size())-1){
        i++;
        if(t[i] < previous_t + t_fineness){
            continue;
        }
        f<<t[i]*Constant::fs_per_au<<" "<<y[i].F.output_densities(this->input_params.Out_F_size())<<endl;
        previous_t = t[i];
    }
    f.close();
}

/**
 * @brief Saves each time and corresponding B-spline coefficients.
 * 
 * @param fname 
 */
void ElectronRateSolver::saveFreeRaw(const std::string& fname) {
    ofstream f;
    cout << "[ Free ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Energy Knot: "<< Distribution::output_knots_eV() << endl;
    f << "# Time (fs) | Expansion Coeffs (not density)"  << endl;

    assert(y.size() == t.size());
    
    for (size_t i=0; i<t.size(); i++) {
        f<<t[i]*Constant::fs_per_au<<" "<<y[i].F<<endl;  // Note that the << operator divides the factors by Constant::eV_per_Ha. -S.P.
    }
    f.close();
}


void ElectronRateSolver::saveBound(const std::string& dir) {
    // saves a table of bound-electron dynamics , split by atom, to folder dir.
    assert(y.size() == t.size());
    // Iterate over atom types
    for (size_t a=0; a<input_params.Store.size(); a++) {
        ofstream f;
        string fname = dir+"dist_"+input_params.Store[a].name+".csv";
        cout << "[ Atom ] Saving to file "<<fname<<"..."<<endl;
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | ";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : input_params.Store[a].index_names) {
            f << cfgname << " ";
        }
        f<<endl;
        // Iterate over time.
        int num_t_points = input_params.Out_T_size();
        if ( num_t_points >  t.size() ) num_t_points = t.size();
        float t_fineness = (simulation_start_time - simulation_end_time)  / num_t_points;    
        float previous_t = t[0];
        int i = -1;
        while (i <  static_cast<int>(t.size())-1){
            i++;
            if(t[i] < previous_t + t_fineness){
                continue;
            }            
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(input_params.Store.size() == y[i].atomP.size());
            
            f<<t[i]*Constant::fs_per_au << ' ' << y[i].atomP[a]<<endl;   // prob. multiplied by 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au
            previous_t = t[i];
        }
        f.close();
    }
}



void ElectronRateSolver::saveBoundRaw(const std::string& dir) {
    for (size_t a=0; a<input_params.Store.size(); a++) {
        string fname = dir+"dist_"+input_params.Store[a].name+"_Raw.csv";
        ofstream f;
        cout << "[ Bound ] Saving to file "<<fname<<"..."<<endl;
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | ";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : input_params.Store[a].index_names) {
            f << cfgname << " ";
        }
        f<<endl;        
        for (size_t i=0; i<t.size(); i++) {
            assert(input_params.Store.size() == y[i].atomP.size());
            f<<t[i]*Constant::fs_per_au << ' ' << y[i].atomP[a]<<endl;
        }
        f.close();
    }
}

/**
 * @brief 
 * @details destructive
 */


void ElectronRateSolver::tokenise(std::string str, std::vector<double> &out, const char delim)
{
    out.resize(0); // clear 
    std::istringstream ss(str);

    std::string s;
    while (std::getline(ss, s, delim)) {
        if (s.size() == 0 || (s.size() == 1 && s[0] == '\r')){continue;}  // This or part is just to deal with excel being a pain.
        double d = stod(s); 
        out.push_back(d);
    }
}

/**
 * @brief Loads the times and corresponding spline factors from a simulation's raw output into this simulation's y[i].F. 
 * The final time step is made as close to the load_time specified by input_params, without being over it or obviously divergent.   
 * @details This code loads the spline scale factors and associated times from the *raw* distribution output file. 
 * The code parasitically overrides the distribution object with the original spline basis coefficients, 
 * after which it transforms to the grid point basis submitted for *this* run via transform_basis().
 * It should be noted that transform_basis() is an *approximation* of the non-polynomial behaviourpyt, albeit it is very good due to using 
 * order 64 gaussian quadrature.
 */
void ElectronRateSolver::loadFreeRaw_and_times() {
    vector<string> time_and_BS_factors;
    const std::string& fname = input_params.Load_Folder() + "freeDistRaw.csv";

    cout << "\n[ Free ] Applying free distribution from file with order 64 gaussian quadrature: "<<fname<<"..."<<endl;
    cout << "[ Caution ] Ensure same atomic input files are used!"<<endl;
    
    ifstream infile(fname);
	if (infile.good())
        std::cout<<"Opened successfully!"<<endl;
    else {
        std::cerr<<"Opening failed."<<endl;
		exit(EXIT_FAILURE); // Quit with a huff.
        return;
    }    
    
    // READ  
    string comment = "#";
    string grid_point_flag = "# Energy Knot:";
        
     // Get indices of lines to load
    std::string line;
    float previous_t;
    float t_fineness = 0.01;  // Fineness should be kept to at least this, so that raw can be kept fine throughout loadings. 
    vector<int> step_indices;
    int i = -4;
    while (std::getline(infile, line)){
        i++;
        if(i >= 0){
            std::istringstream s(line);
            string str_time;
            s >> str_time;    
            float t = stod(str_time);          
            if(i >= 1 && t < previous_t + t_fineness){
                continue;
            }        
            step_indices.push_back(i);
            previous_t =  t;
        }
    }     
    int num_steps = step_indices.size();    
    // get each raw line
    infile.clear();  
    infile.seekg(0, std::ios::beg);
    i = -4;
    std::vector<double> saved_knots;
	while (!infile.eof())
	{    
        i++;
        string line;
        getline(infile, line);
         // KNOTS
        if (!line.compare(0,grid_point_flag.length(),grid_point_flag)){
            // Remove boilerplate
            line.erase(0,grid_point_flag.size()+1);
            // Get grid points (knots)
            this->tokenise(line,saved_knots);
            // Convert to Atomic units
            for (size_t k = 0; k < saved_knots.size();k++){
                saved_knots[k]/=Constant::eV_per_Ha; 
            }
            continue;
        }
        else if (!line.compare(0, 1, comment)){
            continue;
        }
        if (!line.compare(0, 1, "")) continue;
        // continue if i not in step_indices
        if (std::find(step_indices.begin(), step_indices.end(), i) == step_indices.end()) continue;
        time_and_BS_factors.push_back(line);
    }

    // Populate solver with distributions at each time step
    t.resize(num_steps,0);  
    y.resize(num_steps,y[0]);    
    bound_t saved_time(num_steps,0);  // this is a mere stand-in for t at best.
    i = 0;
    std::vector<double> saved_f;  // Defined outside of scope so that saved_f will be from the last time step)
    for(const string elem : time_and_BS_factors){
        // TIME
        std::istringstream s(elem);
        string str_time;
        s >> str_time;
        saved_time[i] = stod(str_time);
        // Convert to right units (based on saveFreeRaw)
        saved_time[i] /= Constant::fs_per_au;

        if(saved_time[i] > input_params.Load_Time_Max() || i >= y.size()){
            // time is past the maximum
            y.resize(i); // (Resized later by integrator for full sim.)
            t.resize(i);
            break;
        }
        // SPLINE FACTORS
        this->tokenise(elem,saved_f);
        saved_f.erase(saved_f.begin());    
        // Convert to old scale.
        for(size_t j = 0; j < saved_f.size();j++){
            saved_f[j] *= Constant::eV_per_Ha;
        }
        // FREE DISTRIBUTION
        std::vector<double> new_knots =  y[0].F.get_knot_energies();      
        y[i].F.set_distribution(saved_knots,saved_f);
        // To ensure compatibility, "translate" old distribution to new grid points.    

        y[i].F.transform_basis(new_knots);         
        t[i] = saved_time[i];
        i++;
    }

    // Translate time to match input of THIS run. 
    if (simulation_start_time != saved_time[0]){
        for(size_t i = 0; i < saved_time.size(); i++){
            t[i] += simulation_start_time - saved_time[0];
        }
    }
    // Shave off end until get a non-obviously-divergent starting point.
    auto too_large = [](vector<state_type>& y){return y.end()[-1].F(0) > y.end()[-2].F(0);};
    auto too_small = [](vector<state_type>& y){return y.end()[-1].F(0) > 0;};
    while (too_large(y) || too_small(y)){
            y.resize(y.size()-1);
            t.resize(t.size()-1);
    }
    // Clear obviously divergent F
    for(size_t i = 0;i++;i < y.size()){
        if (too_large(y) || too_small(y)){
            y[i].F = 0;
        }
    }  
    vector<double> new_knots = y.back().F.get_knot_energies();
    string col = "\033[95m"; string clrline = "\033[0m\n";
    cout <<  col + "Loaded simulation. Param comparisons are as follows:" << clrline
    << "Grid points at simulation resume time (" + col << t.back()*Constant::fs_per_au <<"\033[0m):" << clrline
    << "gp i        | (energy , spline factor)  " << clrline;
    vector<size_t> out_idxs = {0,1, 10,100};
    for(size_t j : out_idxs){
        if (j >= y.size()) continue;  //use setw
        cout << "Source gp "<<j <<" | " + col
        << "(" << saved_knots[j] << " , " << saved_f[j] << ")" << clrline
        << "New gp "<<j <<"   | " + col                                              
        << "(" << new_knots[j] << " , " << y.back().F[j] << ")" 
        << clrline << "------------------" << clrline;
    }
    cout << endl;
}   

/**
 * @brief Loads all times and free e densities from previous simulation's raw output, and uses that to populate y[i].F, the free distribution.
 * @attention Currently assumes same input states
 * @todo need to change load_bound_path to a map from atom names to the specific bound file (currently just takes 1 bound file)
 */
void ElectronRateSolver::loadBound() {
    cout << "Loading atoms' bound states"<<endl; 
    cout << "[ Caution ] Ensure same atoms and corresponding .inp files are used!"<<endl; 


    for (size_t a=0; a<input_params.Store.size(); a++) {
        // (unimplemented) select atom's bound file  
        const std::string& fname = input_params.Load_Folder() + "dist_" + input_params.Store[a].name + "_Raw.csv";

        cout << "[ Free ] Loading bound states from file. "<<fname<<"..."<<endl;
         
        ifstream infile(fname);
        // READ
        // get each raw line
        string comment = "#";
        vector<string> saved_occupancies;
        while (!infile.eof())
        {    
            string line;
            getline(infile, line);
            if (!line.compare(0, 1, comment)) continue;
            if (!line.compare(0, 1, "")) continue;
            saved_occupancies.push_back(line);
        }        
        

        int num_steps = y.size();
        if (y.size()==0){
            cout << y.size() << " <- y.size()" << endl;
            cout << "[[Dev warning]] It seems loadBound was run before loadFreeRaw_and_times, but this means loadBound won't know what times to use." << endl;
        }
        
        //  initialise - fill with initial state
        for(size_t count = 1; count < num_steps; count++){
            this->y[count].atomP[a] = this->y[0].atomP[a];
        }
        
        
        // Iterate through and find each time that matches. 
        int matching_idx;
        for(string elem : saved_occupancies){
            // TIME
            std::stringstream s(elem);
            double elem_time;
            string str_time;
            s >> str_time;            
            elem_time = stod(str_time);
            // Convert to right units (based on saveFreeRaw)
            elem_time /= Constant::fs_per_au;
            if(elem_time > t.back()){
                break;
            }            
            matching_idx = find(t.begin(),t.end(),elem_time) - t.begin(); 
            if (matching_idx >= t.size()){
                continue;
            }
            else{
                std::vector<double> occ_density;
                tokenise(elem,occ_density);
                occ_density.erase(occ_density.begin()); // remove time element
                // Convert to correct units
                const double units = 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;  
                for(size_t k = 0; k < occ_density.size();k++){
                    occ_density[k] /= units;
                }                 
                y[matching_idx].atomP[a] = occ_density;
            }
        }
        // // Shave time and state containers to last matching state.
        //y.resize(matching_idx + 1);
        //t.resize(matching_idx + 1);
        if(t.size() != matching_idx + 1){
            throw std::runtime_error("No bound state found for the final loaded step (assuming program lded free distribution)."); 
        }
    }
}
