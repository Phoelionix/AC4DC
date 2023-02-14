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

    std::vector<double> fake_t;
    size_t num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size();
    size_t t_idx_step = t.size() / num_t_points;
    for (size_t i=0; i<num_t_points; i++) {
        fake_t.push_back(t[i*t_idx_step]);
    }
    pf.save(fake_t,dir+"intensity.csv");

}

void ElectronRateSolver::log_extra_details(ofstream & _log){
    // Saves details pertaining to the simulation's execution to file fname
    if(_log.is_open()){
        cout << "[ Details ] Logging run-specific details..."<<endl;
        _log << endl << "[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << endl;
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
    size_t num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size(); // Fineness of output is only limited by num time steps.
    size_t t_idx_step = t.size() / num_t_points;
    for (size_t i=0; i<num_t_points; i++) {
        f<<t[i*t_idx_step]*Constant::fs_per_au<<" "<<y[i*t_idx_step].F.output_densities(this->input_params.Out_F_size())<<endl;
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
        if (s.size() == 0){continue;}
        double d = stod(s);
        out.push_back(d);
    }
}

/**
 * @brief Loads all times and free e densities from previous simulation's raw output, and uses that to populate y[i].F, the free distribution.
 * @details This code loads the spline scale factors and associated times from the *raw* distribution output file. 
 * The code parasitically overrides the distribution object with the original spline basis coefficients, 
 * after which it transforms to the grid point basis submitted for *this* run via transform_to_new_basis().
 * It should be noted that transform_to_new_basis() is an *approximation*, albeit it is very good due to using 
 * order 64 gaussian quadrature. (Though it seems excessive)
 *  
 * I had no idea what I was doing so it is possible it doesn't need to be an approximation.
 * 
 * @todo this should probably be moved to a new methods file.
 */
void ElectronRateSolver::loadFreeRaw_and_times() {
    vector<string> time_and_BS_factors;
    const std::string& fname = load_free_fname;
    
    
    cout << "[ Free ] Applying free distribution from file with order 10 gaussian quadrature: "<<fname<<"..."<<endl;
    cout << "[ Caution ] Ensure same input files are used!"<<endl;
    
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
        
     // Count num lines
    int num_steps = -3;
    std::string line;
    while (std::getline(infile, line))
        ++num_steps;

    int step_skip_size = 1;
    int max_num_loaded_steps = 5000;
    while(num_steps/step_skip_size > max_num_loaded_steps){
        step_skip_size*=2;
    }            
    num_steps = num_steps/step_skip_size;           
    // get each raw line
    infile.clear();  
    infile.seekg(0, std::ios::beg);
    int count = -4;
    std::vector<double> saved_knots;
	while (!infile.eof())
	{    
        count++;
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
        if (count%(step_skip_size) != 0 || count < 0) continue;
        time_and_BS_factors.push_back(line);
    }

    // Populate solver with distributions at each time step
    t.resize(num_steps,0);  
    y.resize(num_steps,y[0]);    
    bound_t saved_time(num_steps,0);  // this is a mere stand-in for t at best.
    int i = 0;
    std::vector<double> saved_f;  // Defined outside of scope so that saved_f will be from the last time step)
    for(const string elem : time_and_BS_factors){
        // TIME
        std::istringstream s(elem);
        string str_time;
        s >> str_time;
        saved_time[i] = stod(str_time);
        // Convert to right units (based on saveFreeRaw)
        saved_time[i] /= Constant::fs_per_au;

        if(saved_time[i] > loaded_data_time_boundary || i >= y.size()){
            // time is past the maximum
            y.resize(i); // (Resized later by integrator for full sim.)
            t.resize(i);
            vector<double> new_knots = y.back().F.get_knot_energies();
            string col = "\033[95m"; string clrline = "\033[0m\n";
            cout <<  col + "Loaded simulation. Param comparisons are as follows:" << clrline
            << "Gird points at simulation resume time (" + col << t.back()*Constant::fs_per_au <<"\033[0m):" << clrline
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

        y[i].F.transform_to_new_basis(new_knots);         
        t[i] = saved_time[i];
        i++;
    }

    // Translate time to match input of THIS run. 
    if (simulation_start_time != saved_time[0]){
        for(size_t i = 0; i < saved_time.size(); i++){
            t[i] += simulation_start_time - saved_time[0];
        }
    }
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
        size_t num_t_points = input_params.Out_T_size();
        if ( num_t_points >  t.size() ) num_t_points = t.size();
        size_t t_idx_step = t.size() / num_t_points;        
        for (size_t i=0; i<num_t_points; i++) {
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(input_params.Store.size() == y[i].atomP.size());
            
            f<<t[i*t_idx_step]*Constant::fs_per_au << ' ' << y[i*t_idx_step].atomP[a]<<endl;
        }
        f.close();
    }

}

/**
 * @brief Loads all times and free e densities from previous simulation's raw output, and uses that to populate y[i].F, the free distribution.
 * @attention Currently assumes same input states
 * @todo need to change load_bound_fname to a map from atom names to the specific bound file (currently just takes 1 bound file)
 */
void ElectronRateSolver::loadBound() {
    cout << "Loading atoms' bound states"<<endl; 
    cout << "[ Caution ] Ensure same atoms and corresponding .inp files are used!"<<endl; 


    
    for (size_t a=0; a<input_params.Store.size(); a++) {
        // (unimplemented) select atom's bound file  
        const std::string& fname = load_bound_fname; 
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
        
        // Until saving raw files, fill with initial state and just make last element correct.
        for(size_t count = 1; count < num_steps; count++){
            this->y[count].atomP[a] = this->y[0].atomP[a];
        }
        
        std::vector<double> last_occ_density;
        // Iterate backwards until reach a time that matches.
        reverse(saved_occupancies.begin(),saved_occupancies.end());      
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
                continue;
            }            
            matching_idx = find(t.begin(),t.end(),elem_time) - t.begin(); 
            if (matching_idx >= t.size()){
                continue;
            }
            else{
                last_occ_density.resize(0);
                tokenise(elem,last_occ_density);
                last_occ_density.erase(last_occ_density.begin()); // remove time element
                break;
            }
        }
        // Shave time and state containers to the time that matches with the bound state.
        y.resize(matching_idx + 1);
        t.resize(matching_idx + 1);        
        this->y.back().atomP[a] = last_occ_density;

    }
}
