#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::Approximation::PythonBindings;
namespace py = pybind11;


PYBIND11_PLUGIN(pymuqApproximation) {
    py::module m("pymuqApproximation", 
                 "Python bindings for the muqApproximation library.");

    KernelWrapper(m);
    GaussianWrapper(m);
    KLWrapper(m);
    
    return m.ptr();
}
