#ifndef AC4DC_CXX_MOLINP_H
#define AC4DC_CXX_MOLINP_H

#include "Input.h"

class MolInp
{
	// Molecular input for coupled atom/electron plasma calcualtions.
public:
	MolInp(const char* filename, ofstream & log);
	~MolInp() {}

	vector<Input> Atomic;
	vector<RateData::Atom> Store;

	vector<Potential> Pots;
	vector<vector<RadialWF>> Orbits;
	vector<Grid> Latts;
	vector<vector<vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Fluence() {return fluence;} 
	int ini_T_size() {return num_time_steps;}
	double dropl_R() {return radius;}

  	void Set_Fluence(double new_fluence) {fluence = new_fluence;}
	bool Write_Charges() {return write_charges; }
	bool Write_Intensity() {return write_intensity; }
	bool Write_MD_data() {return write_md_data; }

	int Out_T_size() {return out_T_size; }

	string name = "";

    void calc_rates(ofstream &_log, bool recalc=true);

protected:

	double omega = 5000;// XFEL photon energy, au.
	double width = 5; // XFEL pulse width in au. Gaussian profile hardcoded.
	double fluence = 0; // XFEL pulse fluence, 10^4 J/cm^2.
	int num_time_steps = 1000; // Guess number of time steps for time dynamics.
	int out_T_size = 0; // Unlike atomic input, causes to output all points.
	double radius = 1000.;
	int omp_threads = 1;

	bool write_charges = false;
	bool write_intensity = false;
	bool write_md_data = true;

	// unit volume.
	double unit_V = 1.;

	// ONLY USED BY solver.cpp NONTHERMAL SIMULATION
	// NOT READ BY AC4DC SIMULATION
	double min_elec_e = 100;
	double max_elec_e = 5000;
	size_t num_elec_points = 200; // Number of cells in the free-electron distribution expansion

};


#endif
