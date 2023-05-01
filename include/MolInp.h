/**
 * @file MolInp.h
 * @brief Molecular input for coupled atom/electron plasma calculations.
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

#ifndef AC4DC_CXX_MOLINP_H
#define AC4DC_CXX_MOLINP_H

#include "Input.h"
#include "Constant.h"
#include "SplineBasis.h"
#include "GridSpacing.hpp"
#include "LossGeometry.hpp"
#include "Pulse.h"

class MolInp
{
	// Molecular input for coupled atom/electron plasma calculations.
public:
	MolInp(const char* filename, ofstream & _log);
	~MolInp() {}

	/// Vector of atomic input objects
	vector<Input> Atomic; 
	/// Stores all atomic parameters: EII, photoionisation, fluorescence, 
	vector<RateData::Atom> Store; 

	vector<Potential> Pots;
	vector<vector<RadialWF>> Orbits;
	vector<Grid> Latts;
	vector<vector<vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Cutoff_Inputted() {return cutoff_flag;}
	double Simulation_Cutoff() {return simulation_cutoff_time;}
	double Fluence() {return fluence;}
	int Num_Time_Steps() {return num_time_steps;}
	double dropl_R() {return radius;}

  	// void Set_Fluence(double new_fluence) {fluence = new_fluence;}
	bool Write_Charges() {return write_charges; }
	bool Write_Intensity() {return write_intensity; }
	bool Write_MD_data() {return write_md_data; }

	int Out_T_size() {return out_T_size; }
	int Out_F_size() {return out_F_size; }

	int Plasma_Threads(){return omp_threads;}

	string Filtration_File(){return filtration_file;}

	string Load_Folder(){return load_folder;}
	double Load_Time_Max(){return simulation_resume_time_max/Constant::fs_per_au;}
	bool Using_Input_Timestep(){return loading_uses_input_timestep;}

	double Grid_Update_Period(){return grid_update_period;}

	string name = "";

	// Scans for atomic process rates in the folder output/[atom name]/Xsections and recalculates if absent.
	// 'Recalculate' flag skips scanning and forces recomputation.
	// The result is available as Store.
    void calc_rates(ofstream &_log, bool recalc=true);

	GridSpacing elec_grid_type;
	DynamicGridPreset elec_grid_preset;
	Cutoffs param_cutoffs;
	ManualGridBoundaries elec_grid_regions;
	LossGeometry loss_geometry;
	PulseShape pulse_shape = PulseShape::none;

	int num_time_steps = -1; // Number of time steps for time dynamics

	double time_update_gap = 0; // interval between each cout'd time.	 
protected:

	bool validate_inputs();

	double omega = -1;// XFEL photon energy, au.
	double width = -1; // XFEL pulse width in au. Gaussian profile hardcoded.

	bool use_fluence = false;
	bool use_count = false;
	bool use_intensity = false;
	double fluence = -1; // XFEL pulse fluence, au.
	double photon_count = -1;
	double max_intensity = -1;
	
	// Simulation end time, inputted as fs. Note that simulation (w/ rate equations) temporal width is 4*FWHM.
	double simulation_cutoff_time = 1; 
	// Whether the cutoff was given
	bool cutoff_flag = false;

	int out_T_size = 0; // Unlike atomic input, causes to output all points.
	int out_F_size = 100; // Number of F grid points to use when outputting.
	double radius = -1; // Droplet radius
	int omp_threads = 1;

	// Flags for outputting
	bool write_charges = false;
	bool write_intensity = false;
	bool write_md_data = true;

	// unit volume.
	double unit_V = -1.;

	// Simulation loading parameters
	double simulation_resume_time_max; // will attempt to load closest to this time but not after.
	string load_folder = ""; // If "" don't load anything. 
	string filtration_file = ""; // If "" don't load anything.
	bool loading_uses_input_timestep = false; // 

	// Dynamic grid
	double grid_update_period; // time period between dynamic grid updates, fs.
};


#endif
