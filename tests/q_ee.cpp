#include "AdamsIntegrator.hpp"
#include "RateSystem.h"
#include <iostream>
#include <cmath>

using namespace std;

// INTERNAL UNITS:
// eV, fs
class EEcoll : ode::Adams_BM<Distribution>{
public:
    EEcoll(double T, size_t num_points, double _dt, int _order) : Adams_BM<Distribution>(_order){
        //                            num, min, max
        Distribution::set_elec_points(100, 0, 50*T*Constant::kb_Ha);
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

    void print(){
        cout<<"#t | E (eV)"<<endl;
        cout<<"# |  "<<Distribution::output_energies_eV(Distribution::size)<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            cout<<t[i]<<' '<<y[i]<<endl;
        }
    }
protected:
    void sys(const Distribution& q, Distribution& qdot, const double t){
        qdot=0;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(Distribution::size);;
        // q.get_Q_ee(v);
        qdot.applyDelta(v);
    }
};

int main(int argc, char const *argv[]) {

    if (argc <3){
        cerr<<"Usage: integral_verif [T] [step] [order]"<<endl;
        return 1;
    }
    double temperature = atof(argv[1]);
    double step = atof(argv[2]);
    int ord = atoi(argv[3]);
    // cerr<<"[electron-electron] "<<ord<<"th order method, h="<<step<<endl;
    // EEcoll I(temperature, 200, step, ord);
    // I.print();
    return 0;
}
