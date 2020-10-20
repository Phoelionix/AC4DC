#include "Pulse.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include "Constant.h"


void Pulse::save(const std::vector<double>& Tvec, const std::string& fname) {
    std::ofstream f;
    std::cout << "[ Flux ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Time (fs) | Intensity (pht/cm2/fs)" <<endl;
    for (auto& t : Tvec) {
        f << t*Constant::fs_per_au << " ";
        double intensity = (*this)(t);
        intensity /= Constant::fs_per_au;
        intensity /= Constant::cm_per_au*Constant::cm_per_au;
        f << intensity <<std::endl;
    }
    f.close();
}

void GaussianPulse::set_pulse(double fluence, double fwhm) {
    // The photon flux model
    // Gaussian model A e^{-t^2/B}
    std::cout<<"[ Flux ] fluence="<<fluence<<", fwhm="<<fwhm<<endl;
    B = fwhm*fwhm/(4*0.6931471806); // f^2/4ln(2)
    A = fluence/pow(Constant::Pi*B,0.5);
}

inline double GaussianPulse::operator()(double t) {
    // Returns flux at time t (same units as fluence)
    return A*exp(-t*t/B);
}

void SquarePulse::set_pulse(double fluence, double fwhm) {
    // The photon flux model
    // Gaussian model A e^{-t^2/B}
    std::cout<<"[ Flux ] fluence="<<fluence<<", width="<<fwhm<<endl;
    B = fwhm; // f^2/4ln(2)
    A = fluence/fwhm;
}

inline double SquarePulse::operator()(double t) {
    // Returns flux at time t (same units as fluence)
    if (t< -B || t >0) {
        return 0;
    } else {
        return A;
    }
}

