#include "AllClassWrappers.h"

#include "MUQ/Approximation/GaussianProcesses/CoregionalKernel.h"
#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
namespace py = pybind11;

void muq::Approximation::PythonBindings::KernelWrapper(py::module &m)
{
    
}