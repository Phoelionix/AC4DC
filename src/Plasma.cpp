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
#include "Constant.h"
#include <cmath>
#include "Plasma.h"
#include "Dipole.h"

static const elec_state_t zero_state = {0,0,0,0};

Plasma::Plasma(int size)
{
	resize(size);
}

void Plasma::resize(int size)
{
	state.clear();
	state.resize(size, zero_state); // Number of photoelectrons.

	delta.clear();
	delta.resize(size, zero_state); // Number of photoelectrons.

	MaxwellT = 1;
	MaxwellNorm = 1;
}

void Plasma::SetMaxwellPF(double Temperature)
{
	MaxwellT = Temperature;
	MaxwellNorm = pow(0.5/Temperature/Constant::Pi, 1.5);
}

double Plasma::MaxwellPF(double W)
{
	// Multiply by MAxwellian Distribution Norm after the rate is calculated.
	return MaxwellNorm*exp(-W/MaxwellT);
}


double Plasma::BettaInt(double y)
{
	return (-505.214774*y + 384.199825)/(-3033.508326*y*y + 2756.504848*y +
	863.334618) ;
}

// void Plasma::set_last(int m)
// {
// 	state[m] = state[m-1];
// }

void Plasma::update_AB(int m, vector<double>& dt)
{
	elec_state_t st = state[m-1];

	for (int j = 0; j < adams_n; j++) {
		st += delta[m-j-1] * Bashforth_5[j] * dt[m - j - 1];
	}
	state[m] = st;
}

void Plasma::update_AM(int m, vector<double>& dt)
{
	elec_state_t st = state[m-1];
	for (int j = 0; j < adams_n; j++) {
		st += delta[m-j-1] * Moulton_5[j] *  dt[m-j];
	}
	state[m] = st;
}

// TODO: VERIFY IF NORMALISATION IS ACCURATE.
double Plasma::MaxwellEII(double B, double u, int occ)
{
	// Secondary electron impact ionization cross-section.
	// Integral : N(T)*dp*p^3*exp(-p^2/2/T)*sigmaBEB(p, B, u, occ).
	// Find maximum value of p^3*exp(-p^2/2/T).
	double W_min = B;
	double W_max = min(12*MaxwellT + W_min, 300*W_min);// peak to cutoff ration of exp(-12) ~10^(-6).

	double W = 0, Result = 0, k = 0.5*(W_max - W_min), l = 0.5*(W_max + W_min);
	for (int i = 0; i < 13; i++) {
		W = k*gaussX_13[i] + l;
		Result += gaussW_13[i]*W*exp(-W/MaxwellT)*Dipole::sigmaBEB(W, B, u, occ);
	}
	Result *= (W_max - W_min)*MaxwellNorm;// Extra 2 comes from p^2 = 2W in the integrandintegrand.
	return Result;
}
