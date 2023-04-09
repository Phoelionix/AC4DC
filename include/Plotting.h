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



    // Initialize the Python interpreter
    py::scoped_interpreter guard{}; 

    // Import the plot.py module
    //py::module plot = py::module::import("script_generate_frame");
    
    // Import the functions
    // py::object py_plot_frame = py::module_::import("script_generate_frame").attr("plot_frame_from_c");      
    // py::object py_clear_frame = py::module_::import("script_generate_frame").attr("clear_frame_c");   
        
    Plotting() {
        string python_path = std::getenv("WORKSPACE_DIR"); + "/scripts";
        std::wstring wide_py_path = std::wstring(python_path.begin(), python_path.end());
        PySys_SetPath(wide_py_path.c_str());
    // Py_Initialize();
    // // load matplotlib for plotting
    // PyRun_SimpleString(
    //     "from matplotlib import pyplot as plt\n"
    //     "plt.ion()\n"
    //     "plt.show()\n"//"plt.show(block=False)\n"
    //     );

    
    }
    // static void uninitializePlotting() {
    //     PyRun_SimpleString("plt.ioff()\nplt.show()");
    //     Py_Finalize();
    // }
    // static void plot_frame_old(std::vector<double> knot, std::vector<double> density) {
    //     initializePlotting();
    //     PyObject* tmp_string = PyString_FromString((char*)"scripts/script_generate_frame");
    //     PyObject* frame_plot_module = PyImport_Import(tmp_string);

    //     PyObject* plot_frame = PyObject_GetAttrString(frame_plot_module,(char*)"plot_frame");
    //     // chatgpt moment
    //     PyObject* plot_args = PyTuple_New(2);
    //     PyTuple_SetItem(plot_args, 0, PyList_New(1));
    //     PyTuple_SetItem(plot_args, 1, PyList_New(1));
    //     for (int i = 0; i < 100; i++) {
    //         PyObject* x_val = PyFloat_FromDouble(knot[i]);
    //         PyObject* y_val = PyFloat_FromDouble(density[i]);
    //         PyList_SetItem(PyTuple_GetItem(plot_args, 0), i, x_val);
    //         PyList_SetItem(PyTuple_GetItem(plot_args, 1), i, y_val);
    //     }
    //     PyObject_CallObject(myFunction,plot_args);        
    // }
    void plot_frame(std::vector<double> knot, std::vector<double> density){
        std::vector<double> test_vec{1, 2, 3, 4, 5};
        py::list knot_list = py::cast(knot);
        py::list density_list = py::cast(density);     
        py::object py_clear_frame = py::module_::import("script_generate_frame").attr("clear_frame_c")();    
        py::object py_plot_frame = py::module_::import("script_generate_frame").attr("plot_frame_from_c")(knot,density);       
        //py_clear_frame();
        //py_plot_frame(knot,density);
    }
};