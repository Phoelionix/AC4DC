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

/*
 *  author: Alexander Kozlov <alexx.kozloff@gmail.com>
 *
 *  First posted: 20-01-2020
 *  Modified: 24-02-2020 Alaric Sanders <alaric.sanders@gmail.com>
 */

#include "stdafx.h"
#include "Input.h"
#include "Grid.h"
#include "Potential.h"
#include "RadialWF.h"
#include "DecayRates.h"
#include <fstream>
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "Constant.h"
#include <sys/stat.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include "DecayRates.h"
#include "Plasma.h"
#include <assert.h>

using namespace std;

/*
USAGE: ac4dc input.inp
*/


void export_psi(string file_name, RadialWF& Psi, Grid& latt){
	cout<<"Attempting to save orbital "<<file_name<<endl;
	const int N = latt.size();
	assert(Psi.G.size() == N && Psi.F.size() == N);

	std::ofstream os;
	os.open(file_name.c_str());
	os<<"# n="<< Psi.N() << ", l="<< Psi.L() <<std::endl;
	os<<"# r F G"<<std::endl;
	for (size_t i = 0; i < N; i++) {
		os << latt.R(i) <<' '<< Psi.F[i] <<' '<< Psi.G[i] << endl;
	}
	os.close();
}


int main(int argc, char *argv[])
{
	string logname = "./output/log_";
	string filename;
	string tail;

	bool export_orbs = false;

	if (argc > 1) {
		filename = argv[1];
		size_t lastdot = filename.find_last_of(".");
		if (lastdot != std::string::npos) {
			tail = filename.substr(lastdot, filename.size());
			filename = filename.substr(0, lastdot);
		}
		size_t lastslash = filename.find_last_of("/");
		if (lastslash != std::string::npos) filename = filename.substr(lastslash+1);

		logname = logname + filename + ".txt";
	} else {
		std::cout << "Could not find file" << filename << "Exiting..." <<endl;
		return 1;
	}
	if (argc>2 && strcmp(argv[2], "--export_orbs") == 0) export_orbs=true;
	cout<<argv[2]<<endl;

	ofstream log(logname);

	int c_start = clock();

	if (tail == ".txt") {

		// Molecular input.

		MolInp Molecule(argv[1], log);
		Molecule.calc_rates(log);

		// Solve a coupled system of equations for atoms and electron plasma.
		ComputeRateParam Dynamics(Molecule.Latts[0], Molecule.Orbits[0], Molecule.Pots[0], Molecule.Atomic[0]);

		Dynamics.SetupAndSolve(Molecule, log);

	} else {

		// Atomic input.

		Grid Lattice(0);//dummy grid. Will be modified by configuration class
		vector<RadialWF> Orbitals;
		Input Init(argv[1], Orbitals, Lattice, log);

		std::cout << "Nuclear charge: " << Init.Nuclear_Z() << endl;
		std::cout << "Nuclear potential: pointlike Coulomb" << endl << endl;

		Potential U(&Lattice, Init.Nuclear_Z(), Init.Pot_Model());
		HartreeFock HF(Lattice, Orbitals, U, Init, log);


		// Solve the system of equations for atomic charge state dynamics.
		if (Init.TimePts() != 0) {
			ComputeRateParam Dynamics(Lattice, Orbitals, U, Init);
			vector<int> final_occ(Orbitals.size(), 0);
			vector<int> max_occ(Orbitals.size(), 0);

			for (int i = 0; i < max_occ.size(); i++) {
				if (fabs(Orbitals[i].Energy) > Init.Omega()) final_occ[i] = Orbitals[i].occupancy();
				max_occ[i] = Orbitals[i].occupancy();
			}

			Dynamics.SolveFrozen(max_occ, final_occ, log);
			Dynamics.SetupAndSolve(log);
		}

		if (export_orbs){
			string orbfile = "./output/";
			orbfile+=filename+"/Orbital";
			mkdir(orbfile.c_str(), ACCESSPERMS);
			std::cout << "Exporting orbitals..." << endl;
			for ( auto& psi : Orbitals) {
				std::stringstream ss;
				ss<<orbfile;
				ss<<"/N"<<psi.N()<<"L"<<psi.L()<<".csv";
				export_psi(ss.str(), psi, Lattice);
			}
		}
	}


	int c_stop = clock();
	log << "====================================================" << endl
		<<"Total execution time " << float(c_stop - c_start)/1000000 << " sec" << endl;
	log.close();

	return 0;
}
