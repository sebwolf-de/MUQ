#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::Utilities::PythonBindings;
namespace py = pybind11;


PYBIND11_PLUGIN(pymuqUtilities) {
    py::module m("pymuqUtilities", 
                 "Python bindings for the muqUtilities library.");

    //LinearAlgebraWrapper(m);
    //MultiIndicesWrapper(m);
    GeneralUtilitiesWrapper(m);
    
    return m.ptr();
}