#include "src/AdamsIntegrator.hpp"
#include "src/RateSystem.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// INTERNAL UNITS:
// eV, fs
class EEcoll : ode::Adams_BM<Distribution>{
public:
    EEcoll(double Temp, size_t num_points, double _dt, int _num_ptser) : Adams_BM<Distribution>(_num_ptser) {
        //                            num, min, max
        GridSpacing g(GridSpacing::linear);
        Distribution::set_elec_points(100, 0, 6*Temp, g);
        Distribution init;
        // 100 electrons
        cerr<<"Initial conditions: N=100 Maxwellian T="<<Temp*Constant::eV_per_Ha<<" eV, ";
        cerr<<"curve with N=15 delta at "<<4*Temp*Constant::eV_per_Ha<<"eV"<<endl;
        init.set_maxwellian(Temp, 100);
        init.addDeltaSpike(4*Temp, 15);
        // HACK: give it an empty Store
        vector<RateData::Atom> dummy(0);
        Distribution::precompute_Q_coeffs(dummy);
        this->setup(init, _dt);
        this->iterate(0, num_points*_dt); // returns final time
        // this->y and this->t are now populated
    }

    void print(string fname) {
        std::ofstream os(fname);
        os<<"#t | E (eV)"<<endl;
        os<<"# |  "<<Distribution::output_energies_eV(Distribution::size)<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            os<<t[i]*Constant::fs_per_au<<' '<<y[i]<<endl;
        }
    }
protected:
    void sys(const Distribution& q, Distribution& qdot, const double t) {
        qdot=0;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(Distribution::size);;
        q.get_Q_ee(v);
        qdot.applyDelta(v);
        if (isnan(qdot.norm())) throw runtime_error("NaN encountered in sdot");
    }
};

int main(int argc, char const *argv[]) {

    if (argc < 4) {
        cerr<<"Usage: "<<argv[0]<<" [fname] [T (Ha)] [fin_time] [num_points]"<<endl;
        return 1;
    }
    string fname(argv[1]);
    double temperature = atof(argv[2]);
    double fin_time = atof(argv[3]);
    int num_pts = atoi(argv[4]);
    fin_time /= Constant::fs_per_au;
    double step = fin_time/num_pts;
    cerr<<"[electron-electron] "<<num_pts<<"electron points, h="<<step<<endl;
    EEcoll I(temperature, num_pts, step, 3);
    I.print(fname);
    return 0;
}
