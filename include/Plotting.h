// https://stackoverflow.com/questions/8293775/live-graph-for-a-c-application
// https://stackoverflow.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value
// https://stackoverflow.com/questions/28269157/plotting-in-a-non-blocking-way-with-matplotlib


#pragma once

#define PY_SSIZE_T_CLEAN
//#include <Python.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
//#include <dlfcn.h>


namespace py = pybind11;


struct Plotting{    
    Plotting();
    void plot_frame(std::vector<double> knot_to_plot, std::vector<double> density);
    // Initialize the Python interpreter
    py::scoped_interpreter guard{};
};

