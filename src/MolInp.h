#ifndef AC4DC_CXX_MOLINP_H
#define AC4DC_CXX_MOLINP_H

#include "Input.h"
#include "Constant.h"
#include "SplineBasis.h"

class MolInp
{
	// Molecular input for coupled atom/electron plasma calcualtions.
public:
	MolInp(const char* filename, ofstream & log);
	~MolInp() {}

	
	vector<Input> Atomic; // Vector of atomic input objects
	vector<RateData::Atom> Store; // Stores all atomic parameters: EII, photoionisation, fluorescence, 

	vector<Potential> Pots;
	vector<vector<RadialWF>> Orbits;
	vector<Grid> Latts;
	vector<vector<vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Fluence() {return fluence;}
	int Num_Time_Steps() {return num_time_steps;}
	double dropl_R() {return radius;}

  	// void Set_Fluence(double new_fluence) {fluence = new_fluence;}
	bool Write_Charges() {return write_charges; }
	bool Write_Intensity() {return write_intensity; }
	bool Write_MD_data() {return write_md_data; }

	int Out_T_size() {return out_T_size; }
	int Out_F_size() {return out_F_size; }

	double Max_Elec_E() {return max_elec_e;}
	double Min_Elec_E() {return min_elec_e;}
	size_t Num_Elec_Points() {return num_elec_points;}


	string name = "";

	// Scans for atomic process rates in the folder output/[atom name]/Xsections and recalculates if absent.
	// 'Recalculate' flag skips scanning and forces recomputation.
	// The result is available as Store.
    void calc_rates(ofstream &_log, bool recalc=true);

	// Possible values: "linear", "quadratic", "exponential"
	GridSpacing elec_grid_type;

	int num_time_steps = -1; // Number of time steps for time dynamics

protected:

	bool validate_inputs();

	double omega = -1;// XFEL photon energy, au.
	double width = -1; // XFEL pulse width in au. Gaussian profile hardcoded.
	double fluence = -1; // XFEL pulse fluence, au.
	
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

	// Electron grid style
	double min_elec_e = -1;
	double max_elec_e = -1;
	size_t num_elec_points = -1; // Number of cells in the free-electron distribution expansion

};


#endif
