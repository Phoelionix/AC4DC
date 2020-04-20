/*
Solver executable. Reads output of the bin/ratecalculator executable.



*/

#include "RateEquationSolver.h"
#include "SystemSolver.h"
#include "Input.h"
#include "CustomDataType.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

namespace filesys  = boost::filesystem;
namespace popts = boost::program_options;

using namespace std;

struct rate_store {
    vector<Rate> auger;
    vector<Rate> fluor;
    vector<Rate> photo;
    vector<CustomDataType::EIIdata> eii;
};


bool read_data( const string& stem, rate_store& store ){
    stringstream stemstream(stem)
    if (!RateIO::ReadRates(stem.native() + "Auger.txt", store.auger)){
        cerr<<"Could not find Auger rates in "<<stemstream<<"Auger.txt"<<endl;
        return false;
    }

    if (!RateIO::ReadRates(stem.native() + "Fluor.txt", store.fluor)){
        cerr<<"Could not find fluorescence rates in "<<stemstream<<"Fluor.txt"<<endl;
        return false;
    }

    if (!RateIO::ReadRates(stem.native() + "Photo.txt", store.photo)){
        cerr<<"Could not find photoionisation rates in "<<stemstream<<"Photo.txt"<<endl;
        return false;
    }

    if (!RateIO::ReadEIIParams(stem.native() + "EII.json", store.eii)){
        cerr<<"Could not find EII parameters in "<<stemstream<<"EII.json"<<endl;
        return false;
    }
    
    return true;
}


int main(int argc, const char *argv[]) {
    try {
        popts::options_description desc{"Options"};
        desc.add_options()
        ("help,h", "Help screen")
        ("file,f", value<string>(), "Input file name");
        // ("age", value<int>()->notifier(on_age), "Age");

        popts::variables_map vm;
        store(popts::parse_command_line(argc, argv, desc), vm);
        popts::notify(vm);

        if (vm.count("help")) {
            cout << desc << '\n';
            return 1;
        } else if (vm.count("file") == 0) {
            cout << "No input file supplied. Exiting..."
        }
    } catch (const error &ex) {
        cerr << ex.what() << '\n';
    }
    filesys::path pathObj(vm["file"].as<filesys::path>())
    if (!pathObj.has_extension()) {
        cerr<<"File does not have extention. Exiting..."
        return 1;
    }
    string extension = pathObj.extension().native();

    filesys::path logfile = pathObj.parent_path() + "/output/log__" + pathObj.filename() + ".log";
    ofstream log(logfile.native());

    filesys::path ratestem = pathObj.parent_path() + "/output/" + pathObj.filename() + "/Xsections/";

    if (extension == ".txt") {
        // Molecular Input.
        cout<<"Attemting to open "<<pathObj.native()<<"..."<<endl;
        MolInp Init(pathObj.native(), log);

        rate_store store;
        read_data(ratestem.native(), store);



    } else if (extension == ".inp") {
        // Atomic Input.
    }

    SystemSolver calc(T, delta_t);



    return 0;
}
