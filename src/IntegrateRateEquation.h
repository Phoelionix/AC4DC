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
#include <vector>
#include <list>
#include <fstream>
#include "Constant.h"
#include "Plasma.h"

#include <assert.h>

using namespace std;

//Class for solving system of linear coupled diffecrential equations of the form
//
// dP(i)/dt = Sum{j!=i} (R(j->i)P(j) - R(i->j)P(i))
//
//where the rates R(j->i) represent electronic transitions from configuration j to
//configuration i (Photo-ionization, fluorescense, and Auger). Photoionization im
//external beam with a given flux I(t) is assumed. Assumption - electron in the
//continuum orbital is out of the system.

class IntegrateRateEquation
{
	/*
	This class takes photoionisation, Auger and EII cross sections computed by some
	means and solves a system of rate equations for the system dynamics.

	*/
	vector<vector<double>> dpdt;
	vector<vector<double>> p;
	vector<vector<double>> p_storage;

	vector<double> time_storage;
	vector<double>& t;
	vector<double>& dt;
	vector<double> f;
	RateData::Atom& store;

	int adams_n = 5;
public:
	// Rate equations for a single atom.
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, RateData::Atom& Store, vector<double> InitCond, const vector<double>& Intensity = vector<double>());
	int Solve(double P_min = 0, double P_max = 1, int storage_time_pts = 500);
	// Rate equations for single chemical element + electron plasma.
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, RateData::Atom& Store, Plasma & Electrons, vector<double> InitCond, const vector<double>& Intensity = vector<double>());
	int Solve(Plasma & Electrons, double P_min = 0, double P_max = 1, int storage_time_pts = 500);
	// Rate equations for molecule + electron plasma.
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, vector<RateData::Atom> & Store, Plasma & Electrons, const vector<double>& Intensity = vector<double>());
	int Solve(Plasma & Electrons, vector<RateData::Atom> & Store, int storage_time_pts = 500);



	// // Rate equations for single chemical element + nonthermal plasma.
	// IntegrateRateEquation(vector<double> &dT, vector<double> &T, RateData::Atom & Store, NTPlasma & Elecs, vector<double> InitCond, const vector<double>& Intensity = vector<double>());
	// int Solve(NTPlasma & Elecs, double P_min, double P_max, int storage_time_pts = 500);
	// // Rate equations for molecule + non-thermal electron plasma.
	// IntegrateRateEquation(vector<double>& dT, vector<double>& T, vector<RateData::Atom> & Store, NTPlasma & Electrons, const vector<double>& Intensity = vector<double>());
	// int Solve(NTPlasma & Electrons, vector<RateData::Atom> & Store, int storage_time_pts = 500);



	int WriteCharge(vector<RateData::Atom> & Store);
	// P_min and P_max are upper and lower bound for P.
	// storage_time_pts is the number of time points to store. Sometimes calculations might require
	// too many time points, so it is cheaper to store some and interpolate later if needed.

	int Write(ofstream & charge);

	vector<vector<double>> GetP() { return p_storage; }
	vector<double> GetT() { return time_storage; }

	~IntegrateRateEquation();
};






class MonteCarlo
{
	// Monte Carlo implementation of rate equation solver.
};
