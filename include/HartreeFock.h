/**
 * @file HartreeFock.h
 * @brief Entrance to a terrifying rabbit hole.
 * 
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
#pragma once

#include "Input.h"
#include "Grid.h"
#include "RadialWF.h"
#include "Potential.h"
#include <vector>
#include <fstream>
#include "Numerics.h"
#include "Constant.h"


using namespace std;

//Class that performs Hartree-Fock calculations
class HartreeFock
{
public:
//	HamMod = 0 - Hartree-Fock; 1 - LDA.
	HartreeFock(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp, ofstream &log);

	int Get_Virtual(vector<RadialWF> &Virtual, vector<RadialWF> &Orbitals, Potential &U, ofstream &log);
	int LDA_Get_Virtual(vector<RadialWF> &Virtual, vector<RadialWF> &Orbitals, Potential &U, ofstream &log);
	int Master(Grid * Lattice, RadialWF * Current, Potential * U, double Energy_tolerance, ofstream &log);
	// Total configuration energy.
	double Conf_En(vector<RadialWF> &Orbitals, Potential &U);
	// Some electrons occupy virtual orbitals.
	double Conf_En(vector<RadialWF> &Orbitals, vector<RadialWF> &Virtual, Potential &U);
	// Primitive hybridisation induced by uniform field E_at_nuc. Same idea as CI.
	CustomDataType::polarize Hybrid(vector<RadialWF> &Orbitals, vector<RadialWF> &Virtual, Potential &U);


	~HartreeFock();
private:

	double Master_tolerance = pow(10, -10);
	double No_exchange_tolerance = pow(10, -3);
	double HF_tolerance = pow(10, -6);
	int max_HF_iterations = 1000;
	int max_Virt_iterations = 70;
	Grid * lattice;
	double OrthogonalityTest(vector<RadialWF> &Orbitals);
	void MixOldNew(RadialWF * New_Orbital, RadialWF * Old_Orbital);
};

class GreensMethod : Adams
{
public:
	GreensMethod(Grid * Lattice, RadialWF * Current, Potential * U);

	vector<double> GreenOrigin(RadialWF * Psi_O);
	vector<double> GreenInfinity(RadialWF * Psi_Inf);
private:
	Grid * lattice;
	RadialWF * psi;
	Potential * u;
};
