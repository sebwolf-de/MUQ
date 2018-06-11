#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;


PYBIND11_PLUGIN(pymuqSamplingAlgorithms) {
    py::module m("pymuqSamplingAlgorithms", 
                 "Python bindings for the muqSamplingAlgorithms library.");

    KernelWrapper(m);
    ProposalWrapper(m);
    SampleWrapper(m);
    MCMCWrapper(m);
    
    return m.ptr();
}
