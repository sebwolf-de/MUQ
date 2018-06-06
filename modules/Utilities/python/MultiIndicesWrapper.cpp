#include "AllClassWrappers.h"

#include "MultiIndex.h"
#include "MultiIndexFactory.h"
#include "MultiIndexLimiter.h"
#include "MultiIndexSet.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::MultiIndicesWrapper(py::module &m)
{
    
}