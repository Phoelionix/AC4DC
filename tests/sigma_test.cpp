#include "Dipole.h"
#include "Constant.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    if(argc < 3){
        cerr<<"Usage: "<<argv[0]<<" T B"<<endl;
        return 1;
    }
    double T = atof(argv[1]);
    double B = atof(argv[2]);
    double U = B/7;
    double occ = 2;

    cout<<"Total cross-section = "<<Dipole::sigmaBEB(T, B, U, occ)<<endl;
     
    double max_e = T-B;
    double tmp =0;
    for (size_t i = 0; i < 13; i++)
    {
        double W = gaussX_13[i]*max_e/2 + max_e/2;
        tmp += Dipole::DsigmaBEB(T, W, B, U, occ)*gaussW_13[i];
    }
    cout<<"Integrated differential cross-section over [0, T-B] = "<<tmp<<endl;
    
    
    
    return 0;
}
