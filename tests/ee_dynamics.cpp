#include "src/AdamsIntegrator.hpp"
#include "src/RateSystem.h"
#include <iostream>
#include <cmath>

using namespace std;

// INTERNAL UNITS:
// eV, fs
class EEcoll : ode::Adams_BM<Distribution>{
public:
    EEcoll(double T, size_t num_points, double _dt, int _num_ptser) : Adams_BM<Distribution>(_num_ptser) {
        //                            num, min, max
        GridSpacing g(GridSpacing::linear);
        Distribution::set_elec_points(100, 0, 50*T*Constant::kb_Ha, g);
        Distribution init;
        // 100 electrons
        cerr<<"Initial conditions: N=100 Maxwellian T="<<T*Constant::kb_eV<<" eV, ";
        cerr<<"curve with N=100 delta at "<<4*T*Constant::kb_eV<<"eV"<<endl;
        init.set_maxwellian(100, T);
        init.addDeltaSpike(4*T*Constant::kb_Ha, 100);
        this->setup(init, _dt);
        this->iterate(0, num_points); // returns final time
        // this->y and this->t are now populated
    }

    void print() {
        cout<<"#t | E (eV)"<<endl;
        cout<<"# |  "<<Distribution::output_energies_eV(Distribution::size)<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            cout<<t[i]<<' '<<y[i]<<endl;
        }
    }
protected:
    void sys(const Distribution& q, Distribution& qdot, const double t) {
        qdot=0;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(Distribution::size);;
        q.get_Q_ee(v);
        qdot.applyDelta(v);
    }
};

int main(int argc, char const *argv[]) {

    if (argc <3) {
        cerr<<"Usage: "<<argv[0]<<" [T] [step] [num_points]"<<endl;
        return 1;
    }
    double temperature = atof(argv[1]);
    double step = atof(argv[2]);
    int num_pts = atoi(argv[3]);
    cerr<<"[electron-electron] "<<num_pts<<"electron points, h="<<step<<endl;
    EEcoll I(temperature, 100, step, num_pts);
    I.print();
    return 0;
}
