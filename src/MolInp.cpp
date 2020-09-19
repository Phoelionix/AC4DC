#include "MolInp.h"
#include "Constant.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>
#include "HartreeFock.h"
#include "ComputeRateParam.h"


MolInp::MolInp(const char* filename, ofstream & log)
{
	// Input file for molecular ionization calculation.
	map<string, vector<string>> FileContent;

	std::cout<<"Opening molecular file "<<filename<<"... ";
	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	ifstream infile(filename);
	if (infile.good())
        std::cout<<"Success!"<<endl;
    else {
        std::cerr<<"Failed."<<endl;
		exit(EXIT_FAILURE); // chuck a hissy fit and quit.
        return;
    }
	string comment = "//";
	string curr_key = "";

	while (!infile.eof())
	{
		string line;
		getline(infile, line);
		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 1, "")) continue;
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = vector<string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	int num_atoms = FileContent["#ATOMS"].size();

	Orbits.clear();
	Orbits.resize(num_atoms);
	Latts.clear();
	Latts.resize(num_atoms, Grid(0));
	Pots.clear();
	Pots.resize(num_atoms);
	Atomic.clear();
	Store.clear();
	Store.resize(num_atoms);
	Index.clear();
	Index.resize(num_atoms);

	for (int n = 0; n < FileContent["#VOLUME"].size(); n++) {
		stringstream stream(FileContent["#VOLUME"][n]);

		if (n == 0) stream >> unit_V;
		if (n == 1) stream >> radius;
	}

	for (int n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		stringstream stream(FileContent["#OUTPUT"][n]);
		char tmp;

		if (n == 0) stream >> out_T_size;
		if (n == 1) stream >> out_F_size;
		if (n == 2) {
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 3) {
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
		if (n == 4) {
			stream >> tmp;
			if (tmp != 'Y') write_md_data = false;
		}
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
	}

	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		stringstream stream(FileContent["#NUMERICAL"][n]);

		if (n == 0) stream >> num_time_steps;
		if (n == 1) stream >> omp_threads;
		if (n == 2) stream >> min_elec_e;
		if (n == 3) stream >> max_elec_e;
		if (n == 4) stream >> num_elec_points;
		if (n == 5) stream >> elec_grid_type;
		if (n == 6) stream >> elec_grid_type.transition_E;
		if (n == 7) stream >> elec_grid_type.num_exp;

	}
	const string bc = "\033[97;1m"; // begin colour escape code
	const string clr = "\033[0m"; // clear escape code
	const string banner = "\033[34m[ ================================ ]\033[0m";
	cout<<banner<<endl;
	cout<<bc<<"Unit cell size: "<<clr<<unit_V<<"A^3"<<endl;
	cout<<bc<<"Droplet radius: "<<clr<<radius<<"A"<<endl<<endl;

	cout<<bc<<"Photon energy: "<<clr<<omega<<" eV"<<endl;
	cout<<bc<<"Pulse fluence: "<<clr<<fluence*10000<<"J/cm^2"<<endl;
	cout<<bc<<"Pulse FWHM: "<<clr<<width<<"fs"<<endl<<endl;

	cout<<bc<<"Electron grid: "<<clr<<min_elec_e<<" ... "<<max_elec_e<<" eV"<<endl;
	cout<<"               "<<num_elec_points<<" points"<<endl;
	cout<<bc<<"Grid type: "<<clr<<elec_grid_type<<endl;
	cout<<endl;

	cout<<bc<<"ODE Iteration: "<<clr<<num_time_steps<<" timesteps"<<endl<<endl;

	cout<<bc<<"Output parameters: "<<clr<<endl;
	cout<<bc<<"Timesteps: "<<clr<<out_T_size<<endl;
	cout<<bc<<"Energy points: "<<clr<<out_F_size<<endl;
	cout<<banner<<endl;

	// Convert to number of photon flux.
	omega /= Constant::eV_per_Ha;
	fluence *= 10000/Constant::Jcm2_per_Haa02/omega;

	elec_grid_type.transition_E /= Constant::eV_per_Ha;

	// Convert to atomic units.
	width /= Constant::fs_per_au;
	radius /= Constant::Angs_per_au;
	unit_V /= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;

	min_elec_e /= Constant::eV_per_Ha;
	max_elec_e /= Constant::eV_per_Ha;

	// Reads the very top of the file, expecting input of the form
	// H 2
	// O 1
	// Then scans for atomic files of the form 
	// input/H.inp
	// input/O.inp
	// Store is then populated with the atomic data read in below.
	for (int i = 0; i < num_atoms; i++) {
		string at_name;
		double at_num;

		stringstream stream(FileContent["#ATOMS"][i]);
		stream >> at_name >> at_num;

		Store[i].nAtoms = at_num/unit_V;
		Store[i].name = at_name;
		Store[i].R = radius;

		at_name = "input/" + at_name + ".inp";

		Atomic.push_back(Input((char*)at_name.c_str(), Orbits[i], Latts[i], log));
		// Overrides pulses found in .inp files
		Atomic.back().Set_Pulse(omega, fluence, width);
		Atomic.back().Set_Num_Threads(omp_threads);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}

	if (!validate_inputs()) {
		cerr<<endl<<endl<<endl<<"Exiting..."<<endl;
		throw runtime_error(".mol input file is invalid");
	}
}

bool MolInp::validate_inputs() {
	bool is_valid=true;
	cerr<<"\033[31;1m";
	if (omega <= 0 ) { cerr<<"ERROR: pulse omega must be positive"; is_valid=false; }
	if (width <= 0 ) { cerr<<"ERROR: pulse width must be positive"; is_valid=false; }
	if (fluence <= 0 ) { cerr<<"ERROR: pulse fluence must be positive"; is_valid=false; }
	if (num_time_steps <= 0 ) { cerr<<"ERROR: got negative number of timesteps"; is_valid=false; }
	if (out_T_size <= 0) { cerr<<"ERROR: system set to output zero timesteps"; is_valid=false; }
	if (out_F_size <= 0) { cerr<<"ERROR: system set to output zero energy grid points"; is_valid=false; }
	if (radius <= 0) { cerr<<"ERROR: radius must be positive"; is_valid=false; }
	if (omp_threads <= 0) { omp_threads = 4; cerr<<"Defaulting number of OMP threads to 4"; }

	if (elec_grid_type.mode == GridSpacing::unknown){cerr<<"ERROR: Grid spacing (must start with (l)inear, (q)uadratic, (e)xponential, (h)ybrid) not recognised"; is_valid=false;}
	if (elec_grid_type.mode == GridSpacing::hybrid){
		if (elec_grid_type.num_exp <= 0 || elec_grid_type.num_exp >= num_elec_points) { 
			cerr<<"Defaulting number of exponential points to "<<num_elec_points/2;
			elec_grid_type.num_exp = num_elec_points/2;
		}
		if (elec_grid_type.transition_E <= 0 || elec_grid_type.transition_E > max_elec_e) {
			 cerr<<"Defaulting crossover energy to "<<max_elec_e/5; 
		}
	}

	// unit cell volume.
	if (unit_V <= 0) { cerr<<"ERROR: unit xell volume must be positive"; is_valid=false; }

	// Electron grid style
	if(min_elec_e < 0 || max_elec_e < 0 || max_elec_e <= min_elec_e) { cerr<<"ERROR: Electron grid specification invalid"; is_valid=false; }
	if (num_time_steps <= 0 ) { cerr<<"ERROR: got negative number of energy steps"; is_valid=false; }
	cerr<<"\033[0m";
	return is_valid;
}

void MolInp::calc_rates(ofstream &_log, bool recalc) {
	// Loop through atomic species.
	for (int a = 0; a < Atomic.size(); a++) {
		HartreeFock HF(Latts[a], Orbits[a], Pots[a], Atomic[a], _log);

		// This Computes the parameters for the rate equations to use, loading them into Init.
		ComputeRateParam Dynamics(Latts[a], Orbits[a], Pots[a], Atomic[a], recalc);
		vector<int> final_occ(Orbits[a].size(), 0);
		vector<int> max_occ(Orbits[a].size(), 0);
		for (int i = 0; i < max_occ.size(); i++) {
			if (fabs(Orbits[a][i].Energy) > Omega()) final_occ[i] = Orbits[a][i].occupancy();
			max_occ[i] = Orbits[a][i].occupancy();
		}

		string name = Store[a].name;
		double nAtoms = Store[a].nAtoms;

		Store[a] = Dynamics.SolvePlasmaBEB(max_occ, final_occ, _log);
		Store[a].name = name;
		Store[a].nAtoms = nAtoms;
		Store[a].R = dropl_R();
		Index[a] = Dynamics.Get_Indexes();
	}
}
