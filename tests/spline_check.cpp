#include "src/SplineBasis.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc < 5) {
        cout<<"Usage: "<<argv[0]<<" [fstem] [num_funcs] [min_e (Ha)] [max_e (Ha)] (linear|quadratic|exponential) "<<endl;
    }
    std::string dummy(argv[1]);
    ofstream fout(dummy+"_vals.csv");
    ofstream dfout(dummy+"_derivative.csv");
    int num_funcs = atoi(argv[2]);
    double min = atof(argv[3]);
    double max = atof(argv[4]);
    istringstream is(argv[5]);
    GridSpacing gs;
    is >> gs;

    BasisSet basis;
    basis.set_parameters(num_funcs, min, max, gs, 0);
    double de = max/(num_funcs*20);
    double e=0;
    for (size_t i=0; i<num_funcs*20; i++) {
        fout<<e<<" ";
        dfout<<e<<" ";
        for (size_t j=0; j<num_funcs; j++) {
            fout<<basis(j, e)<<" ";
            dfout<<basis.D(j, e)<<" ";
        }
        fout<<endl;
        dfout<<endl;
        e += de;
    }
    fout.close();
    dfout.close();
    return 0;
}
