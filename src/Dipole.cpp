#include "Dipole.h"
#include <cmath>
#include "Constant.h"

double Dipole::sigmaBEB(double T, double B, double u, int occ)
{
	// Electron impact ionization cross-section
	// in Binary Encounter Bethe approximation : Kim, Rudd, PRA 50(5), 3954 (1994)
	double t = T/B;
	if (t < 1) return 0;
	double lnt = log(t);
	double S = Constant::Pi*occ/B/B;

	//double Result = S/(t + u + 1)*(0.5*(1 - 1./t/t)*lnt + (1 - 1/t - lnt/(1+t)));
	double Result = S/(t + u + 1)*(0.5*(1 - 1./t/t)*lnt + (1 - 1/t - lnt/(1+t)));
	return Result;
}


double Dipole::sigmaBEBw1(double T, double B, double u, int occ)
{
	// Calculates the integral \int_0^(p*p/2 - B) dW W d(sigma)/dW using gaussian quadrature.
	double Result = 0;
	double Wmax = 0.5*(T - B);
	// Critically important to capture the interval where d(sigma)/dW substantially != 0.
	if (Wmax <= 0) return 0;
	double W = 0, k = 0.5*Wmax;
	// -1 .. 1 -> 0 .. Wmax [ W = 0.5*Wmax*(x + 1) ]
	for (int i = 0; i < 13; i++) {
		W = k*(gaussX_13[i] + 1);
		Result += gaussW_13[i]*W*Dipole::DsigmaBEB(T, W, B, u, occ);
	}

	Result *= 0.5*Wmax;
	return Result;
}


double Dipole::DsigmaBEB(double T, double W, double B, double u, int occ)
{
	// T - impactor electron energy.
	// W - ejected electron energy.
	// B - binding energy.
	// i - orbital index from which W is ejected.
	// Q[i] are set to 1. See commented functions if need to be calculated.
	// see eq. (52)
	double t = T/B;
	double w = W/B;
	if (t < 1 + w) return 0;

	double F = 1.;//F2 first, Q = 1
	double x = 1./(w + 1);
	double y = 1./(t - w);
	double Result = -1*(x + y)/(t+1) + (x*x + y*y) + log(t)*(x*x*x + y*y*y);

	Result *= Constant::Pi*occ/B/B/(t + u + 1);
	return Result;
}
