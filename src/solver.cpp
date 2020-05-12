/*
Solver executable. Reads output of the bin/ratecalculator executable.



*/

#include "RateEquationSolver.h"
#include "ElectronSolver.h"
#include "Input.h"
#include "Constant.h"
#include <iostream>

using namespace std;

// Rate system solver.
// Uses precomputed rates from AC4DC for all atomic cross-section data.
// KEEP IN MIND:
// - For every atom X listed in the .mol file, AC4DC must be run for the file X.inp
// - AC4DC has input parameters for pulse width, energy and fluence.
// - Only photon energy affects the rate calculations.
// Let scripts/run.py handle all of these details.

int get_file_names(const char* infile_, string &name, string &logfile) {
    string infile = string(infile_);
    size_t namestart = infile.rfind('/');
    size_t nameend = infile.rfind('.');
    namestart = (namestart==string::npos) ? 0 : namestart + 1;// Exclude leading slash
    nameend = (nameend==string::npos) ? infile.size() : nameend;
    name = infile.substr(namestart, nameend-namestart);
    logfile = "output/run_" + name + ".log";
    // check correct format
    string extension = infile.substr(nameend);
    if (extension != ".mol") {
        cerr<<"This file is for coupled calculations. Please provide a .mol file similar to Lysozyme.mol"<<endl;
        return 1;
    }
    return 0;
}



int main(int argc, const char *argv[]) {

    bool recalc = false;
    if (argc < 2){
        cout << "No input file supplied. Exiting...";
        return 1;
    }
    if (argc > 2){
        // flags
        if (argv[2][0] != '-') goto YES_THIS_IS_A_GOTO_SUE_ME;
        char c;
        int i = 0;
        while ((c = argv[2][++i]) != '\0'){
            switch (c) {
                case 'r':
                    // recalculate
                    recalc = true;
                    break;
                default:
                    cout<<"Flag '"<<c<<"' is not a recognised flag."<<endl;
            }
        }
    }
    YES_THIS_IS_A_GOTO_SUE_ME:

    string in_name, logname;

    if (get_file_names(argv[1], in_name, logname) == 1) return 1;

    cout<<"Running simulation for target "<<in_name<<endl;
    cout << "logfile name: " << logname <<endl;
    ofstream log(logname);
    cout << "Initialising... " <<endl;
    ElectronSolver S(argv[1], log); // Contains all of the collision parameters.

    S.compute_cross_sections(log, recalc);
    S.solve();



    return 0;
}
