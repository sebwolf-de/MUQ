#include "AllClassWrappers.h"

#include ".h"
#include "KarhunenLoeveFactory.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
namespace py = pybind11;

void muq::Approximation::PythonBindings::KLWrapper(py::module &m)
{
    
}