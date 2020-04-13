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
#include <fstream>
#include "Constant.h"

using namespace std;

struct elec_state_t {
	double N;
	double E;
	double Np;
	double Ep;

	elec_state_t operator+(const elec_state_t &a){
		return {a.N+N, a.E+E, a.Np+Np, a.Ep+Ep};
	}

	// Scalar multiplication
	elec_state_t operator*( double c){
		return {c*N, c*E, c*Np, c*Ep};
	}

	// entrywise multiplication
	elec_state_t operator*(const elec_state_t a){
		return {a.N*N, a.E*E, a.Np*Np, a.Ep*Ep};
	}

	elec_state_t operator=(double c){
		return {c,c,c,c};
	}

	elec_state_t operator=(int x){
		double c = x;
		return {c,c,c,c};
	}
	void operator+=(const elec_state_t &a){
		*this = *this + a;
	}
	void operator*=(const double c){
		*this = *this * c;
	}
	void operator*=(const elec_state_t &a){
		*this = *this * a;
	}
	bool operator==(const elec_state_t &a){
		return N==a.N && E==a.E && Np==a.Np && Ep == a.Ep;
	}

};


class Plasma
{
//======================================================================
//	Electron impact ionization (EII) and three body recombination (TBR)
//	cross-sections from binary encounter models (BEB and BED).
//  Also stores plasma paramameters, such as photo- and secondary electron number densities and energy densities.
//  BEB is much cheaper since it doesn't require
//  differential ocillatory strength at every point.
//  For binary encounter model see Yong-Ki Kim, M. Eugene Rudd, PRA 50(5), 3954 (1994).
//  For sum rules used in derivation see A. Dalgarno, N. Lynn, Proc. Phys. Soc. A 70, 802 (1957).
//  For oscillator strengths see H. A. Bethe, E. E. Salpeter book, 265 (1977).
//======================================================================
public:
	Plasma(int size);

	~Plasma() {}

	vector<elec_state_t> state; // State vector
	vector<elec_state_t> delta; // First derivatives

	void resize(int size);

	// p - impactor electron momentum
	// p_s - secondary electron momentum
	// B - ionization potential
	// occ - orbitals occupancy
	double DsigmaBEB(double T, double W, double B, double t, int occ);
	double sigmaBEB(double T, double B, double u, int occ);
	// Int_0^(p*p/2 - B) dW W d(sigmaBED)/dW - first moment of a secondary electron energy.
	double sigmaBEBw1(double T, double B, double u, int occ);

	void SetMaxwellPF(double Temperature); // Set temperature and norm for Maxwellian PF.
	double MaxwellPF(double W);
	double MaxwellEII(double B, double u, int occ);
	double TESTMaxwellEII(double B, double u, int occ);


	// Adams-Bashforth and Adams-Moulton single-step methods
	void update_AB(int m, vector<double>& dt);
	void update_AM(int m, vector<double>& dt);

	double BettaInt(double y);// Auxillary function for photo-secondary electrons energy exchange.

	// void set_last(int m); // Sets variables[i]  = variables[i-1]

	//void setup_EII(vector<RadialWF> &Virtual, double k_min, double k_max);
private:
	double MaxwellT = 1;
	double MaxwellNorm = 1;

	int adams_n = 5;

	//void Get_Ni(Grid & Lattice, vector<RadialWF> & Orbitals, vector<RadialWF> &Virtuals);
	//void Get_Qi(Grid & Lattice, vector<RadialWF> & Orbitals, vector<RadialWF> &Virtuals );
	// Calculate "N_i" (BED) and "Q_i" (BEB), and orbitals electron kinetic energies parameters EII model.
	//vector<double> N;
	//vector<double> Q;
};
