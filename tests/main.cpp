/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/

// (C) Alaric Sanders 2020

#include "ComputeRateParam.hpp"
#include "ElectronSolver.hpp"
#include "Input.hpp"
#include "Constant.hpp"
#include <iostream>

using namespace std;

// Rate system solver.
// Uses precomputed rates from AC4DC for all atomic cross-section data.
// KEEP IN MIND:
// - For every atom X listed in the .mol file, AC4DC must be run for the file X.inp
// - AC4DC has input parameters for pulse width, energy and fluence.
// - Only photon energy affects the rate calculations.
// Let scripts/run.py handle all of these details.


void print_banner(const char* fname){
    std::ifstream ifs(fname, ifstream::in);

    char c = ifs.get();
    while (ifs.good()) {
        std::cout << c;
        c = ifs.get();
    }
    ifs.close();
}

void try_mkdir(const std::string& fname) {
    if (mkdir(fname.c_str(), ACCESSPERMS) == -1) {
        if (errno != EEXIST)
            std::cerr<<"mkdir error attempting to create "<< fname << ":" << errno;
    }
}

int get_file_names(const char* infile_, std::string &tag, std::string &logfile, std::string&outdir) {
    // Takes infile of the form "DIR/Lysozyme.mol"
    // Stores "Lysozyme" in tag, "output/log/run_Lysozyme" in logfile
    std::string infile = string(infile_);
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
    outdir = "output/__Molecular/"+tag+"/";
    try_mkdir(outdir);
    std::string extension = infile.substr(tagend);
    if (extension != ".mol") {
        std::cerr<<"This file is for coupled calculations. Please provide a .mol file similar to Lysozyme.mol"<<"\n";
        return 1;
    }
    return 0;
}

int main(int argc, const char *argv[]) {
    CmdParser runsettings(argc, argv);
    if (!runsettings.valid_input) {
        return 1;
    }

    std::string name, logname, outdir;

    std::cout<<"Copyright (C) 2020  Alaric Sanders and Alexander Kozlov"<<"\n";
    std::cout<<"Welcome to ee_dynamics, Spencer will be servicing you today with high quality electron content."<<"\n";
    std::cout<<"This is free software, and you are welcome to redistribute it"<<"\n";

    if (get_file_names(argv[1], name, logname, outdir) == 1)
        return 1;

    std::cout<<"Running simulation for target "<<name<<"\n";
    std::cout << "\033[1;32mInitialising... \033[0m" <<"\n";
    ElectronSolver S(argv[1], log); // Contains all of the collision parameters.
    std::cout << "\033[1;32mComputing cross sections... \033[0m" <<"\n";
    S.compute_cross_sections(log, runsettings.recalc);
    if (runsettings.solve_rate_eq) {
        std::cout << "\033[1;32mSolving rate equations... \033[0m" <<"\n";
        S.solve();
        std::cout << "\033[1;32mDone! \033[0m" <<"\n";
        S.save(outdir);
    } else {
        std::cout << "\033[1;32mDone! \033[0m" <<"\n";
    }

    return 0;
}
