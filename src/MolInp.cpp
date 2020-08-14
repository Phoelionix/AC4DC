#include "MolInp.h"
#include "Constant.h"
#include <fstream>
#include <sstream>
#include <iostream>
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
		if (n == 1) {
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 2) {
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
		if (n == 3) {
			stream >> tmp;
			if (tmp != 'Y') write_md_data = false;
		}
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
		if (n == 3) stream >> num_time_steps;
	}

	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		stringstream stream(FileContent["#NUMERICAL"][n]);

		if (n == 0) stream >> num_time_steps;
		if (n == 1) stream >> omp_threads;
		if (n == 2) stream >> min_elec_e;
		if (n == 3) stream >> max_elec_e;
		if (n == 4) stream >> num_elec_points;

	}

	// Convert to number of photon flux.
	omega /= Constant::eV_per_Ha;
	fluence *= 10000/Constant::Jcm2_per_Haa02/omega;

	// Convert to atomic units.
	width /= Constant::fs_per_au;
	radius /= Constant::Angs_per_au;
	unit_V = 1./Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;

	min_elec_e /= Constant::eV_per_Ha;
	max_elec_e /= Constant::eV_per_Ha;

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
		Atomic.back().Set_Pulse(omega, fluence, width);
		Atomic.back().Set_Num_Threads(omp_threads);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}

}

void MolInp::calc_rates(ofstream &_log, bool recalc){
	// Loop through atomic species.
	for (int a = 0; a < Atomic.size(); a++) {
		HartreeFock HF(Latts[a], Orbits[a], Pots[a], Atomic[a], _log);

		// This Computes the parameters for the rate equations to use, loading them into Init.
		ComputeRateParam Dynamics(Latts[a], Orbits[a], Pots[a], Atomic[a], false);
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
