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



// Note. this code is NOT compatible with Windows.
// Rewriting with boost::filesystem is a good idea if this is required.

void try_mkdir(const std::string& fname){
    if (mkdir(fname.c_str(), ACCESSPERMS) == -1) {
        if (errno != EEXIST)
            cerr<<"mkdir error attempting to create "<< fname << ":" << errno;
    }
}

int get_file_names(const char* infile_, string &tag, string &logfile, string&outdir) {
    // Takes infile of the form "DIR/Lysozyme.mol"
    // Stores "Lysozyme" in tag, "output/log/run_Lysozyme" in logfile
    string infile = string(infile_);
    size_t tagstart = infile.rfind('/');
    size_t tagend = infile.rfind('.');
    tagstart = (tagstart==string::npos) ? 0 : tagstart + 1;// Exclude leading slash
    tagend = (tagend==string::npos) ? infile.size() : tagend;
    tag = infile.substr(tagstart, tagend-tagstart);
    // guarantee the existence of a folder structure
    try_mkdir("output");
    try_mkdir("output/log");
    try_mkdir("output/__Molecular");
    logfile = "output/log/run_" + tag + ".log";
    // check correct format
    outdir = "output/__Molecular";
    string extension = infile.substr(tagend);
    if (extension != ".mol") {
        cerr<<"This file is for coupled calculations. Please provide a .mol file similar to Lysozyme.mol"<<endl;
        return 1;
    }
    return 0;
}

class CmdParser{
public:
    CmdParser(int argc, const char *argv[]){
        if (argc < 2){
            cout << "Usage: solver path/to/molecular/in.mol [-rh]";
            valid_input = false;
        }
        for (int a=2; a<argc; a++){
            if (argv[a][0] != '-')
                continue;

            int i=0;
            while (argv[a][i++]!='\0'){
                switch (argv[a][i]) {
                    case 'r':
                        // recalculate
                        recalc = true;
                        break;
                    case 'h':
                        // Usage help.
                        cout<<"This is physics code, were you really expecting documentation?"<<endl;
                        break;
                    default:
                        cout<<"Flag '"<<argv[a][i]<<"' is not a recognised flag."<<endl;
                }
            }
        }
    }
    bool recalc = false;
    bool valid_input = true;
};


int main(int argc, const char *argv[]) {
    CmdParser runsettings(argc, argv);
    if (!runsettings.valid_input){
        return 1;
    }

    string name, logname, outdir;

    if (get_file_names(argv[1], name, logname, outdir) == 1)
        return 1;

    cout<<"Running simulation for target "<<name<<endl;
    cout << "logfile name: " << logname <<endl;
    ofstream log(logname);
    cout << "Initialising... " <<endl;
    ElectronSolver S(argv[1], log); // Contains all of the collision parameters.

    S.compute_cross_sections(log, runsettings.recalc);
    S.solve();
    S.save(outdir);


    return 0;
}
