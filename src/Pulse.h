#ifndef AC4DC_PULSE_CXX_H
#define AC4DC_PULSE_CXX_H

#include <vector>
#include <string>

class Pulse
{
public:
    Pulse() {};
    Pulse(double fluence, double w) {
        this->set_pulse(fluence, w);
    };
    virtual void set_pulse(double fluence, double width);
    virtual double operator()(double t); // Yields Photon flux in units of A/fwhm
    void save(const std::vector<double>& T, const std::string& file);
};

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


#endif