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

void Pulse::set_pulse(double fluence, double _fwhm) {
    // The photon flux model
    // Gaussian model A e^{-t^2/B}
    std::cout<<"[ Flux ] fluence="<<fluence<<", fwhm="<<fwhm<<endl;
    I0 = fluence/_fwhm;
    fwhm = _fwhm;

}

double Pulse::operator()(double t) {
    // Returns flux at time t (same units as fluence)
    const double norm = sqrt(Constant::Pi/4/log(2));
    switch (shape)
    {
    case PulseShape::gaussian:
        return I0/fwhm/norm*pow(2,-t*t*4/fwhm*fwhm);
        break;
    case PulseShape::square:
        if (t< -fwhm || t >0) {
            return 0;
        } else {
            return I0;
        }
        break;
    default:
        throw runtime_error("Pulse shape has not been set.");
        break;
    }   
}
