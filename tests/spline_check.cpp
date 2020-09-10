#include "SplineBasis.h"
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc < 6){
        cout<<"Usage: "<<argv[0]<<" [num_funcs] [min_e (Ha)] [max_e (Ha)] [zero_regularity] (linear|quadratic|exponential "<<endl;
    }
    int num_funcs = atoi(argv[1]);
    double min = atof(argv[2]);
    double max = atof(argv[3]);
    int zero_degree = atoi(argv[4]);
    istringstream is(argv[5]);
    GridSpacing gs;
    is >> gs;

    BasisSet basis;
    basis.set_parameters(num_funcs, min, max, zero_deg, gs);
    double de = max/(num_funcs*20);
    double e=0;
    for (size_t i=0; i<num_funcs*20; i++){
        cout<<e<<" ";
        for (size_t j=0; j<num_funcs; j++){
            cout<<basis(j, e)<<" ";
        }
        cout<<endl;
        e += de;
    }
    return 0;
}
