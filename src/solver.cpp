/*
Solver executable. Reads output of the bin/ratecalculator executable.



*/

#include "RateEquationSolver.h"
#include "SystemSolver.h"
#include "Input.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

namespace filesys  = boost::filesystem;
namespace popts = boost::program_options;

int main(int argc, const char *argv[]) {
    try {
        popts::options_description desc{"Options"};
        desc.add_options()
        ("help,h", "Help screen")
        ("file,f", value<std::string>(), "Input file name");
        // ("age", value<int>()->notifier(on_age), "Age");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("help")) {
            std::cout << desc << '\n';
            return 1;
        } else if (vm.count("file") == 0) {
            std::cout << "No input file supplied. Exiting..."
        }
    } catch (const error &ex) {
        std::cerr << ex.what() << '\n';
    }
    filesys::path pathObj(vm["file"].as<filesys::path>())
    if (!pathObj.has_extension()) {
        std::cerr<<"File does not have extention. Exiting..."
        return 1;
    }
    std::string extension = pathObj.extension().string();
    if (extension == ".txt") {
        // Molecular Input.
        

    } else if (extension == ".inp") {
        // Atomic Input.
    }

    SystemSolver calc(T, delta_t);



    return 0;
}
