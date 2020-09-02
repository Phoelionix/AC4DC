#include "SplineBasis.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[]) {
    int num_funcs = atoi(argv[1]);
    double max = atof(argv[2]);
    BasisSet basis;
    basis.set_parameters(num_funcs, 0, max);
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
