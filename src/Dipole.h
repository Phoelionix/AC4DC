#ifndef AC4DC_CXX_DIPOLE_H
#define AC4DC_CXX_DIPOLE_H


namespace Dipole
{
	// p - impactor electron momentum
	// p_s - secondary electron momentum
	// B - ionization potential
	// occ - orbitals occupancy
	double DsigmaBEB(double T, double W, double B, double u, int occ);
	double sigmaBEB(double T, double B, double u, int occ);
	// Int_0^(T - B) dW W d(sigmaBED)/dW - total energy absorbance
	double sigmaBEBw1(double T, double B, double u, int occ);
}


#endif /* end of include guard: AC4DC_CXX_DIPOLE_H */
