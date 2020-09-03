#include "FreeDistribution.h"
#include "Constant.h"
#include <iostream>

using namespace std;


RateData::EIIdata get_fake_eii(){
    RateData::EIIdata tmp;
    tmp.init = 10;
    tmp.fin.push_back(4);
    tmp.occ.push_back(2);
    tmp.ionB.push_back(11.3);
    tmp.kin.push_back(1.9);
}

int main(int argc, char const *argv[]) {
    if (argc < 5) cerr << "usage: "<<argv[0]<<" [density (au^-3)] [temp (eV)] [max_e (eV)] [num_basis]";
    double density = atof(argv[1]);
    double T = atof(argv[2])*Constant::kb_Ha;
    double max_e = atof(argv[3])/Constant::eV_per_Ha;
    int N = atoi(argv[4]);

    Distribution::set_elec_points(N, 0, max_e);
    Distribution F;
    F.set_maxwellian(density, T);
    RateData::EIIdata eii_process = get_fake_eii();
    eiiGraph eg;
    Distribution::basis.Gamma_eii(eg, eii_process, J);
    cout<<Distribution::output_energies_eV(N*10)<<endl;
    cout<<F.output_densities(N*10)<<endl;
    return 0;
}
