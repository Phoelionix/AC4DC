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
    EEcoll(double max_e, size_t num_e_points) : Adams_BM<Distribution>(3) {
        //                            num, min, max
        GridSpacing g(GridSpacing::linear);
        g.zero_degree = 1;
        Distribution::set_elec_points(num_e_points, 0, max_e, g);
        
        // HACK: give it an empty Store
        vector<RateData::Atom> dummy(0);
        Distribution::precompute_Q_coeffs(dummy);
        
    }

    void set_initial_condition(double E, double density, double spike_density, double _dt) {
        Distribution init;
        cerr<<"Initial conditions: n="<<density<<" Maxwellian T="<<E*Constant::eV_per_Ha<<" eV, ";
        cerr<<"curve with n="<<spike_density<<" delta at "<<4*E*Constant::eV_per_Ha<<"eV"<<endl;
        init.set_maxwellian(E, density);
        init.addDeltaSpike(4*E, spike_density);
        this->setup(init, _dt);
        // this->y and this->t are now sized appropriately
    }

    void run_sim(double final_time){
        this->iterate(0, final_time);
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

    if (argc < 7) {
        cerr<<"Usage: "<<argv[0]<<" [fname] [T (eV)] [fin_time (fs)] [num_t_points] [num_e_points] [density (A^-3)]"<<endl;
        return 1;
    }
    string fname(argv[1]);
    double temperature = atof(argv[2]);
    double fin_time = atof(argv[3]);
    int num_t_pts = atoi(argv[4]);
    int num_e_pts = atoi(argv[5]);
    double density = atof(argv[6]);

    cerr<<"Density ="<<density<<" elec per Angstrom3"<<endl;
    cerr<<"Final Time: "<<fin_time<<" fs in "<<num_t_pts<<" timesteps"<<endl;

    fin_time /= Constant::fs_per_au;
    temperature /= Constant::eV_per_Ha;
    density *= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;

    double step = fin_time/num_t_pts;
    EEcoll I(6*temperature, num_e_pts);
    cerr<<"final time: "<<fin_time<<" au"<<endl;
    I.set_initial_condition(temperature, density, density*0.2, step);
    I.run_sim(fin_time);
    I.print(fname);
    return 0;
}
