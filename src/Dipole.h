#ifndef AC4DC_CXX_DIPOLE_H
#define AC4DC_CXX_DIPOLE_H


namespace Dipole
{
	// p - impactor electron momentum
	// p_s - secondary electron momentum
	// B - ionization potential
	// occ - orbitals occupancy
	double DsigmaBEB(double T, double W, double B, double t, int occ);
	double sigmaBEB(double T, double B, double u, int occ);
	// Int_0^(p*p/2 - B) dW W d(sigmaBED)/dW - first moment of a secondary electron energy.
	double sigmaBEBw1(double T, double B, double u, int occ);
};


#endif /* end of include guard: AC4DC_CXX_DIPOLE_H */
