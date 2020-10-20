#ifndef AC4DC_PULSE_CXX_H
#define AC4DC_PULSE_CXX_H

#include <vector>
#include <string>
#include <iostream>

enum class PulseShape {
    square, gaussian, none
};


class Pulse
{
public:
    Pulse() {};
    Pulse(double fluence, double _fwhm, PulseShape ps) {
        this->set_pulse(fluence, _fwhm);
        shape=ps;
    };
    double operator()(double t); // Yields Photon flux in units of A/fwhm
    void save(const std::vector<double>& T, const std::string& file);
    void set_pulse(double fluence, double width);
    void set_shape(PulseShape ps){
        if (ps != PulseShape::none){
            shape = ps;
        } else {
            std::cerr<<"[ Pulse ] Defaulting pulse shape to Gaussian"<<std::endl;
            shape=PulseShape::gaussian;
        }  
    };
private:
    PulseShape shape;
    double I0;
    double fwhm;
};
/*
class GaussianPulse : public Pulse
{
public:
    // GaussianPulse() {};
    GaussianPulse(double fluence, double fwhm) {
        set_pulse(fluence, fwhm);
    };
    void set_pulse(double fluence, double width);
    inline double operator()(double t); // Yields Photon flux in same units as A
    // void save(const vector<double>& T, const std::string& file);
private:
    double A;
    double B;
};

class SquarePulse : public Pulse
{
public:
    SquarePulse(double fluence, double fwhm) {
        set_pulse(fluence, fwhm);
    };
    void set_pulse(double fluence, double width);
    inline double operator()(double t); // Yields Photon flux in same units as A
    // void save(const vector<double>& T, const std::string& file);
private:
    double A;
    double B;
};
*/

namespace {
    
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, PulseShape ps) {
        switch (ps)
        {
        case PulseShape::gaussian:
            os << "Gaussian";
            break;
        case PulseShape::square:
            os << "rectangular";
            break;
        default:
            os << "Unknown pulse shape";
            break;
        }
        return os;
    }

    [[maybe_unused]] std::istream& operator>>(std::istream& is, PulseShape& ps) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No pulse shape provided, defaulting to Gaussian..."<<std::endl;
            ps = PulseShape::gaussian;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'g':
            ps = PulseShape::gaussian;
            break;
        case 's':
        case 'r':
            ps = PulseShape::square;
            break;
        default:
            std::cerr<<"Unrecognised pulse shape \""<<tmp<<"\""<<std::endl;
            ps = PulseShape::gaussian;
            break;
        }
        return is;
    }
}

#endif